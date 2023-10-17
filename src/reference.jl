@lazy mutable struct Reference
    const genomes::Set{Genome}
    # Sequence name => (sequence, targets)
    const targets_by_name::Dict{String, Tuple{Sequence, Vector{Target}}}
    const clades::Vector{Vector{Clade{Genome}}}
    @lazy shortest_seq_len::Int
    @lazy fraction_assembled::Float64
end

"""
    Reference(::Union{IO, AbstractString}; [min_seq_length=1])

A `Reference` contains the ground truth to benchmark against.
Conceptually, it consists of the following parts:
* A list of genomes, each with sources
* The full taxonomic tree, as lists of clades
* A list of sequences, each with a list of (source, span) to where it maps.

Normally, the types `FlagSet` `Genome`, `Source`, `Clade` and `Sequence` do not
need to be constructed manually, but are constructed when the `Reference` is loaded
from a JSON file.

# Examples
```jldoctest
julia> ref = Reference(path_to_ref_file; min_seq_length=3);

julia> ref isa Reference
true

julia> length(genomes(ref))
3

julia> nseqs(ref)
11

julia> first(ref.genomes) isa Genome
true
```

See also: [`subset`](@ref), [`Genome`](@ref), [`Clade`](@ref)
"""
Reference

function Reference(::Unsafe)
    Reference(
        Set{Genome}(),
        Dict{String, Tuple{Sequence, Vector{Target}}}(),
        Vector{Clade{Genome}}[],
        uninit,
        uninit,
    )
end

Base.show(io::IO, ::Reference) = print(io, "Reference()")
function Base.show(io::IO, ::MIME"text/plain", x::Reference)
    if get(io, :compact, false)
        show(io, x)
    else
        print(
            io,
            "Reference",
            "\n  Genomes:    ",
            length(genomes(x)),
            "\n  Sequences:  ",
            nseqs(x),
            "\n  Ranks:      ",
            nranks(x),
            "\n  Seq length: ",
            x.shortest_seq_len,
            "\n  Assembled:  ",
            round(x.fraction_assembled * 100; digits=1),
            " %",
        )
    end
end

top_clade(x::Reference) = only(last(x.clades))
genomes(x::Reference) = x.genomes
nseqs(x::Reference) = length(x.targets_by_name)
function nranks(x::Reference)
    isempty(x.genomes) && return 0
    length(x.clades) + 1
end

function finish!(ref::Reference)
    @isinit(ref.fraction_assembled) && return ref
    foreach(finish!, ref.genomes)
    assembly_size = genome_size = 0
    for genome in ref.genomes
        assembly_size += genome.assembly_size
        genome_size += genome.genome_size
    end
    shortest_seq_len =
        minimum(i -> length(first(i)), values(ref.targets_by_name); init=typemax(Int))
    shortest_seq_len == typemax(Int) &&
        error("Cannot initialize a Reference with no sequences")
    @init! ref.shortest_seq_len = shortest_seq_len
    @init! ref.fraction_assembled = assembly_size / genome_size
    ref
end

"""
    subset!(
            ref::Reference;
            sequences::Function=Returns(true),
            genomes::Function=Returns(true)
    )::Reference

Mutate `ref` in place, removing genomes and sequences.
Keep only sequences S where `sequences(S)` returns `true` and genomes G for which
`genomes(G)` returns `true`.

See also: [`subset`](@ref), [`Reference`](@ref)

# Examples
```jldoctest
julia> ref
Reference
  Genomes:    3
  Sequences:  11
  Ranks:      3
  Seq length: 10
  Assembled:  66.8 %

julia> subset(ref; genomes=g -> Flags.organism in flags(g))
Reference
  Genomes:    2
  Sequences:  11
  Ranks:      3
  Seq length: 10
  Assembled:  91.3 %

julia> VambBenchmarks.subset(ref; sequences=s -> length(s) ≥ 25)
Reference
  Genomes:    3
  Sequences:  9
  Ranks:      3
  Seq length: 25
  Assembled:  61.1 %
```
"""
function subset!(
    ref::Reference;
    sequences::Function=Returns(true),
    genomes::Function=Returns(true),
)::Reference
    ref = uninit!(ref)
    genomes_to_remove = Genome[]
    sources_to_remove = Set{Source}()
    filter!(ref.targets_by_name) do (_, v)
        sequences(first(v))::Bool
    end
    filter!(ref.genomes) do g
        keep = genomes(g)::Bool
        keep || push!(genomes_to_remove, g)
        keep
    end
    for genome in genomes_to_remove
        union!(sources_to_remove, genome.sources)
    end
    for (_, targets) in values(ref.targets_by_name)
        filter!(targets) do (source, _)
            !in(source, sources_to_remove)
        end
    end
    for genome in genomes_to_remove
        recursively_delete_child!(genome)
    end
    for i in (length(ref.clades) - 1):-1:1
        empty!(ref.clades[i])
        for parent in ref.clades[i + 1]
            union!(ref.clades[i], parent.children::Vector{Clade{Genome}})
        end
    end
    for genome in ref.genomes, source in genome.sources
        filter!(i -> sequences(first(i))::Bool, source.sequences)
    end
    finish!(ref)
end

# TODO: This deepcopy is quite slow, and it would be nice to optimise it.
# However, manual deepcopying of references is quite error prone, and having
# an incomplete deepcopy could lead to nasty bugs, so I'll just eat it for now.
"""
    subset(ref::Reference; kwargs...)

Non-mutating copying version of `subset`.

See also: [`subset!`](@ref)
"""
subset(ref::Reference; kwargs...) = subset!(deepcopy(ref); kwargs...)

function uninit!(ref::Reference)
    @uninit! ref.fraction_assembled
    @uninit! ref.shortest_seq_len
    for genome in ref.genomes
        @uninit! genome.genome_size
        @uninit! genome.assembly_size
        for source in genome.sources
            @uninit! source.assembly_size
        end
    end
    ref
end

function add_genome!(ref::Reference, genome::Genome)
    in(genome, ref.genomes) && error(lazy"Genome $(genome.name) already in reference")
    push!(ref.genomes, genome)
    ref
end

function add_sequence!(ref::Reference, seq::Sequence, targets::Vector{Target})
    for (source, span) in targets
        add_sequence!(source, seq, span)
    end
    if last(get!(ref.targets_by_name, seq.name, (seq, targets))) !== targets
        error(lazy"Duplicate sequence in reference: $(seq.name)")
    end
    ref
end

function parse_bins(
    io::IO,
    ref::Reference,
    binsplit_sep::Union{Nothing, AbstractString, Char}=nothing,
)::Vector{Bin}
    [
        Bin(binname, seqs, ref.targets_by_name) for
        (binname, seqs) in parse_bins(io, Dict, ref, binsplit_sep)
    ]
end

function parse_bins(
    io::IO,
    ::Type{Dict},
    ref::Reference,
    binsplit_sep::Union{Nothing, AbstractString, Char}=nothing,
)::Dict{<:AbstractString, Vector{Sequence}}
    itr = tab_pairs(eachline(io))
    if binsplit_sep !== nothing
        itr = binsplit_tab_pairs(itr, binsplit_sep)
    end
    seqs_by_binname = Dict{SubString{String}, Vector{Sequence}}()
    for (binname, seqname) in itr
        (seq, _) = ref.targets_by_name[seqname]
        push!(get!(valtype(seqs_by_binname), seqs_by_binname, binname), seq)
    end
    seqs_by_binname
end

const JSON_VERSION = 1
struct ReferenceJSON
    version::Int
    # [(name, flags, [(sourcename, length)])]
    genomes::Vector{Tuple{String, Int, Vector{Tuple{String, Int}}}}
    # [Sequence => sequence_length, [(subject, from, to)]]
    sequences::Vector{Tuple{String, Int, Vector{Tuple{String, Int, Int}}}}
    # [[(child, parent)], [(parent, grandparent)] ...]
    taxmaps::Vector{Vector{Tuple{String, Union{String, Nothing}}}}
end
StructTypes.StructType(::Type{ReferenceJSON}) = StructTypes.Struct()

function Reference(path::AbstractString; min_seq_length::Integer=1)
    open_perhaps_gzipped(i -> Reference(i; min_seq_length), String(path))
end

function Reference(io::IO; min_seq_length::Integer=1)
    Reference(JSON3.read(io, ReferenceJSON), Int(min_seq_length))
end

function Reference(json_struct::ReferenceJSON, min_seq_length::Int)
    if json_struct.version != JSON_VERSION
        @warn (
            "Deserializing reference JSON of version $(json_struct.version), " *
            "but the supported version of the currently loaded version of VambBenchmarks is $(JSON_VERSION)."
        )
    end
    ref = Reference(unsafe)

    # Parse genomes
    for (genomename, flags, sourcesdict) in json_struct.genomes
        genome = Genome(genomename, FlagSet(UInt(flags)))
        add_genome!(ref, genome)
        for (source_name, source_length) in sourcesdict
            add_source!(genome, source_name, source_length)
        end
    end

    # Check for unique sources
    source_by_name = Dict{String, Source{Genome}}()
    for genome in ref.genomes, source in genome.sources
        if haskey(source_by_name, source.name)
            error(
                lazy"Duplicate source: \"$(source.name)\" belongs to both genome ",
                lazy"\"$(source_by_name[source.name].genome.name)\" and \"$(genome.name)\".",
            )
        end
        source_by_name[source.name] = source
    end

    # Parse sequences
    for (seq_name, seq_length, targs) in json_struct.sequences
        seq_length ≥ min_seq_length || continue
        targets = map(targs) do (source_name, from, to)
            source = get(source_by_name, source_name, nothing)
            if source === nothing
                error(
                    lazy"Sequence \"$(seq_name)\" maps to source \"$(source_name)\", but no such source in reference",
                )
            end
            (source, Int(from):Int(to))
        end
        seq = Sequence(seq_name, seq_length)
        add_sequence!(ref, seq, targets)
    end

    # Add taxonomy
    copy!(ref.clades, parse_taxonomy(ref.genomes, json_struct.taxmaps))

    # Finalize the reference
    finish!(ref)
end

function save(io::IO, ref::Reference)
    json_dict = Dict{Symbol, Any}()
    json_dict[:version] = JSON_VERSION
    # Genomes
    json_dict[:genomes] = [
        (genome.name, Int(genome.flags.x), [(s.name, s.length) for s in genome.sources]) for genome in ref.genomes
    ]

    # Sequences
    json_dict[:sequences] = [
        (
            seq.name,
            seq.length,
            [(source.name, first(span), last(span)) for (source, span) in targets],
        ) for (_, (seq, targets)) in ref.targets_by_name
    ]

    # Taxmaps
    taxmaps = [
        Tuple{String, Union{String, Nothing}}[
            (genome.name, genome.parent.name) for genome in ref.genomes
        ],
    ]
    json_dict[:taxmaps] = taxmaps
    for clades in ref.clades
        v = eltype(taxmaps)()
        push!(taxmaps, v)
        for clade in clades
            parent = clade.parent
            value = parent === nothing ? nothing : parent.name
            push!(v, (clade.name, value))
        end
    end

    JSON3.write(io, json_dict)
end

function parse_taxonomy(
    genomes::Set{Genome},
    dict::Vector{Vector{Tuple{String, Union{String, Nothing}}}},
)::Vector{Vector{Clade{Genome}}}
    child_by_name = Dict{String, Node}(g.name => g for g in genomes)
    parent_by_name = empty(child_by_name)
    result = Vector{Vector{Clade{Genome}}}()
    for (rank, taxmap) in enumerate(dict)
        cladeset = Clade{Genome}[]
        for (child_name, maybe_parent_name) in taxmap
            # Get child
            child = get(child_by_name, child_name, nothing)
            child === nothing && error(
                "At rank $rank, found child name \"$child_name\", but this does not exist " *
                "on the previous rank.",
            )
            # Create parent if it does not already exist
            parent_name = maybe_parent_name === nothing ? child_name : maybe_parent_name
            parent =
                get(parent_by_name, parent_name, nothing)::Union{Clade{Genome}, Nothing}
            if parent === nothing
                parent = Clade(parent_name, child)
                parent_by_name[parent_name] = parent
                push!(cladeset, parent)
            else
                add_child!(parent, child)
                child.parent = parent::Clade{Genome}
            end
        end
        # Enforce all childrens have parents
        for child in values(child_by_name)
            isdefined(child, :parent) ||
                error("At rank $rank, child $(child.name) has no parent")
        end
        parent_by_name, child_by_name = child_by_name, parent_by_name
        empty!(parent_by_name)
        push!(result, cladeset)
    end
    # If there isn't exactly one remaining clade at the top, make one.
    if length(child_by_name) != 1 || only(values(child_by_name)) isa Genome
        top = Clade("top", first(values(child_by_name)))
        for child in values(child_by_name)
            child === first(top.children) || add_child!(top, child)
        end
        push!(result, [top])
    else
        top = only(values(child_by_name))::Clade{Genome}
        # If multiple redundant top clades, go down to the useful level
        while nchildren(top) == 1 && only(top.children).name == top.name
            top = only(top.children)::Clade{Genome}
            top.parent = nothing
            pop!(result)
        end
    end
    @assert length(last(result)) == 1
    foreach(i -> sort!(i; by=j -> j.name), result)
    return result
end
