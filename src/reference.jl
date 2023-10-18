@lazy mutable struct Reference
    const genomes::Set{Genome}
    # Sequence name => index into targets
    const target_index_by_name::Dict{String, UInt32}
    const targets::Vector{Tuple{Sequence, Vector{Target}}}
    const clades::Vector{Vector{Clade{Genome}}}
    @lazy shortest_seq_len::Int
    @lazy fraction_assembled::Float64
end

"""
    Reference(::Union{IO, AbstractString})

A `Reference` contains the ground truth to benchmark against.
Conceptually, it consists of the following parts:
* A list of genomes, each with sources
* The full taxonomic tree, as lists of clades
* A list of sequences, each with a list of (source, span) to where it maps.

Normally, the types `FlagSet` `Genome`, `Source`, `Clade` and `Sequence` do not
need to be constructed manually, but are constructed when the `Reference` is loaded
from a JSON file.

A `Reference` is loaded from a JSON file, which is specified either as an `IO`,
or its path as an `AbstractString`. If the path ends with `.gz`, automatically
gzip decompress when reading the file.

# Examples
```jldoctest
julia> ref = Reference(path_to_ref_file);

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
        Dict{String, Int}(),
        Vector{Tuple{Sequence, Vector{Target}}}(),
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
nseqs(x::Reference) = length(x.targets)
function nranks(x::Reference)
    isempty(x.genomes) && return 0
    length(x.clades) + 1
end

function finish!(ref::Reference)
    @isinit(ref.fraction_assembled) && return ref
    scratch = Tuple{Int, Int}[]
    foreach(g -> finish!(g, scratch), ref.genomes)
    assembly_size = genome_size = 0
    for genome in ref.genomes
        assembly_size += genome.assembly_size
        genome_size += genome.genome_size
    end
    shortest_seq_len = minimum(i -> length(first(i)), ref.targets; init=typemax(Int))
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
  Assembled:  61.9 %

julia> subset(ref; genomes=g -> Flags.organism in flags(g))
Reference
  Genomes:    2
  Sequences:  11
  Ranks:      3
  Seq length: 10
  Assembled:  91.3 %

julia> BinBencherBackend.subset(ref; sequences=s -> length(s) ≥ 25)
Reference
  Genomes:    3
  Sequences:  9
  Ranks:      3
  Seq length: 25
  Assembled:  56.2 %
```
"""
function subset!(
    ref::Reference;
    sequences::Function=Returns(true),
    genomes::Function=Returns(true),
)::Reference
    ref = uninit!(ref)

    # Cache the sequences and genomes to remove in order to not
    # need to compute the predicates multiple times
    genomes_to_remove = Genome[]
    sources_to_remove = Set{Source}()

    # Update both ref.targets and ref.target_index_by_name
    mask = BitVector(sequences(first(i)) for i in ref.targets)
    new_idx = cumsum(mask)
    keepat!(ref.targets, mask)
    filter!(kv -> mask[last(kv)], ref.target_index_by_name)
    map!(i -> new_idx[i], values(ref.target_index_by_name))

    # Populate genomes_to_remove and remove the genomes from ref.genomes
    filter!(ref.genomes) do g
        keep = genomes(g)::Bool
        keep || push!(genomes_to_remove, g)
        keep
    end

    # Compute sources to remove, and remove the sources from ref.targets
    for genome in genomes_to_remove
        union!(sources_to_remove, genome.sources)
    end
    for (_, targets) in ref.targets
        filter!(targets) do (source, _)
            !in(source, sources_to_remove)
        end
    end

    # Remove references to deleted genomes from parents
    # (and recursively remove now-empty Clades)
    for genome in genomes_to_remove
        recursively_delete_child!(genome)
    end

    # Remove references to Clades we deleted from the Clade tree,
    # but which may still be present in ref.clades
    for i in (length(ref.clades) - 1):-1:1
        empty!(ref.clades[i])
        for parent in ref.clades[i + 1]
            union!(ref.clades[i], parent.children::Vector{Clade{Genome}})
        end
    end

    # Remove sequences from sources
    for genome in ref.genomes, source in genome.sources
        filter!(i -> sequences(first(i))::Bool, source.sequences)
    end

    # Re-initialize the now completely filtered Reference
    finish!(ref)
end

# This deepcopy is quite slow, and it would be nice to optimise it.
# However, manual deepcopying of references is quite error prone, and having
# an incomplete deepcopy could lead to nasty bugs, so I'll just eat it.
"""
    subset(ref::Reference; kwargs...)

Non-mutating copying version of `subset!`.
This is currently much slower than `subset!`.

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
            @uninit! source.total_bp
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
    L = length(ref.target_index_by_name)
    if L == typemax(UInt32)
        error("References can only hold 4294967295 sequences")
    end
    i = get!(ref.target_index_by_name, seq.name, L + 1)
    if i != L + 1
        error(lazy"Duplicate sequence in reference: $(seq.name)")
    end
    push!(ref.targets, (seq, targets))
    ref
end

function parse_bins(
    io::IO,
    ::Type{Dict},
    ref::Reference,
    binsplit_sep::Union{Nothing, AbstractString, Char}=nothing,
    disjoint::Bool=true,
)::Dict{<:AbstractString, <:Vector{<:Integer}}
    lines = eachline(io)
    header = "clustername\tcontigname"
    it = iterate(lines)
    if (isnothing(it) ? nothing : rstrip(first(it))) != header
        error(lazy"Expected following header line in cluster file: $(repr(header))")
    end
    itr = tab_pairs(lines)
    itr = isnothing(binsplit_sep) ? itr : binsplit_tab_pairs(itr, binsplit_sep)
    seen_indices = falses(length(ref.targets))
    idxs_by_binname = Dict{SubString{String}, Vector{UInt32}}()
    @inbounds for (binname, seqname) in itr
        i = ref.target_index_by_name[seqname]
        if seen_indices[i]
            name = first(ref.targets[i]).name
            error(lazy"Sequence \"$(name)\" seen twice in disjoint Binning")
        end
        seen_indices[i] = true
        push!(get!(valtype(idxs_by_binname), idxs_by_binname, binname), i)
    end
    idxs_by_binname
end

const JSON_VERSION = 2
# [(name, flags, [(sourcename, length)])]
const GENOMES_JSON_T = Vector{Tuple{String, Int, Vector{Tuple{String, Int}}}}
# [Sequence => sequence_length, [(subject, from, to)]]
const SEQUENCES_JSON_T = Vector{Tuple{String, Int, Vector{Tuple{String, Int, Int}}}}                   # [[(child, parent)], [(parent, grandparent)] ...]
const TAXMAPS_JSON_T = Vector{Vector{Tuple{String, String}}}

struct ReferenceJSON
    version::Int
    genomes::GENOMES_JSON_T
    sequences::SEQUENCES_JSON_T
    taxmaps::TAXMAPS_JSON_T
end
StructTypes.StructType(::Type{ReferenceJSON}) = StructTypes.Struct()

function Reference(path::AbstractString)
    open_perhaps_gzipped(i -> Reference(i), String(path))
end

function Reference(io::IO)
    Reference(JSON3.read(io, ReferenceJSON))
end

function Reference(json_struct::ReferenceJSON)
    if json_struct.version != JSON_VERSION
        @warn (
            "Deserializing reference JSON of version $(json_struct.version), " *
            "but the supported version of the currently loaded version of BinBencherBackend is $(JSON_VERSION)."
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
        targets = map(targs) do (source_name, from, to)
            source = get(source_by_name, source_name, nothing)
            if source === nothing
                error(
                    lazy"Sequence \"$(seq_name)\" maps to source \"$(source_name)\", but no such source in reference",
                )
            end
            (source, (Int(from), Int(to)))
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
    genomes::GENOMES_JSON_T = [
        (genome.name, Int(genome.flags.x), [(s.name, s.length) for s in genome.sources]) for genome in ref.genomes
    ]
    json_dict[:genomes] = genomes

    # Sequences
    sequences::SEQUENCES_JSON_T = [
        (
            seq.name,
            seq.length,
            [(source.name, first(span), last(span)) for (source, span) in targets],
        ) for (seq, targets) in ref.targets
    ]

    json_dict[:sequences] = sequences

    # Taxmaps
    taxmaps::TAXMAPS_JSON_T =
        [[(genome.name, genome.parent.name) for genome in ref.genomes]]
    json_dict[:taxmaps] = taxmaps

    for clades in ref.clades
        length(clades) == 1 && isnothing(only(clades).parent) && break
        v = eltype(taxmaps)()
        push!(taxmaps, v)
        for clade in clades
            parent = clade.parent::Clade
            push!(v, (clade.name, parent.name))
        end
    end

    JSON3.write(io, json_dict)
end

# Invariants for the input list:
# The input is a vector of ranks, consisting of (child, parent) names
# On every rank, the children is a unique set (each child present once)
# of the parent of the previous rank:
#   - On the first rank, the children is a unique set of genome names
#   - On the last rank, the parent names can be arbitrary.
#     If there is more than a single unique name, then a top node will be added
function parse_taxonomy(
    genomes::Set{Genome},
    taxmaps::Vector{Vector{Tuple{String, String}}},
)::Vector{Vector{Clade{Genome}}}
    child_by_name = Dict{String, Node}(g.name => g for g in genomes)
    assigned_children_names = Set{String}()
    parent_by_name = empty(child_by_name)
    result = Vector{Vector{Clade{Genome}}}()
    for (rank, taxmap) in enumerate(taxmaps)
        cladeset = Clade{Genome}[]
        empty!(assigned_children_names)
        for (child_name, parent_name) in taxmap
            child_name ∈ assigned_children_names && error(
                lazy"At rank $(rank), child \"$(child_name)\" is present more than once.",
            )
            # Get child
            child = get(child_by_name, child_name, nothing)
            child === nothing && error(
                lazy"At rank $(rank), found child name \"$(child_name)\", but this does not exist " *
                "on the previous rank.",
            )
            # Create parent if it does not already exist
            parent =
                get(parent_by_name, parent_name, nothing)::Union{Clade{Genome}, Nothing}
            if parent === nothing
                parent = Clade(parent_name, child)
                parent_by_name[parent_name] = parent
                push!(cladeset, parent)
            else
                add_child!(parent, child)
            end
            push!(assigned_children_names, child_name)
        end
        if length(assigned_children_names) != length(child_by_name)
            missing_child = first(setdiff(keys(child_by_name), assigned_children_names))
            error(lazy"At rank $(rank), child $(missing_child) has no parent")
        end
        parent_by_name, child_by_name = child_by_name, parent_by_name
        empty!(parent_by_name)
        push!(result, cladeset)
    end
    # Create top node, if needed
    if length(last(result)) > 1
        top = Clade("top", first(last(result)))
        for child in @view(result[end][2:end])
            add_child!(top, child)
        end
        push!(result, [top])
    end
    foreach(i -> sort!(i; by=j -> j.name), result)
    return result
end
