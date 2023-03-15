@lazy mutable struct Reference
    const genomes::Set{Genome}
    # Sequence name => (sequence, targets)
    const targets_by_name::Dict{String, Tuple{Sequence, Vector{Target}}}
    const clades::Vector{Vector{Clade{Genome}}}
    @lazy fraction_assembled::Float64
end

function Reference()
    Reference(
        Set{Genome}(),
        Dict{String, Tuple{Sequence, Vector{Target}}}(),
        Vector{Clade{Genome}}[],
        uninit
    )
end

Base.show(io::IO, ::Reference) = print(io, "Reference()")
function Base.show(io::IO, ::MIME"text/plain", x::Reference)
    if get(io, :compact, false)
        show(io, x)
    else
        print(io,
            "Reference",
            "\n  Genomes:   ", ngenomes(x),
            "\n  Sequences: ", nseqs(x),
            "\n  Ranks:     ", nranks(x),
            "\n  Assembled: ", round(x.fraction_assembled * 100; digits=1), " %"
        )
    end
end

top_clade(x::Reference) = only(last(x.clades))
ngenomes(x::Reference) = length(x.genomes)
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
    @init! ref.fraction_assembled = assembly_size / genome_size
    ref
end

"""
    filter_size(ref::Reference, size::Int)::Reference

Create a new, independent (deep copied) `Reference`, where all sequences
with a length smaller than `size` has been removed.

# Examples
```julia
julia> ref
Reference
  Genomes:   1057
  Sequences: 1247324
  Ranks:     8

julia> filter_size(ref, 5000)
Reference
  Genomes:   1057
  Sequences: 30501
  Ranks:     8
```
"""
function filter_size(ref::Reference, size::Int)
    ref = deepcopy(ref)
    @uninit! ref.fraction_assembled
    filter!(ref.targets_by_name) do (_, v)
        first(v).length ≥ size
    end
    for genome in ref.genomes
        @uninit! genome.genome_size
        @uninit! genome.assembly_size
        for source in genome.sources
            @uninit! source.assembly_size
            filter!(((seq, _),) -> seq.length ≥ size, source.sequences)
        end
    end
    finish!(ref)
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
    binsplit_sep::Union{Nothing, AbstractString, Char}=nothing
)
    itr = tab_pairs(eachline(io))
    if binsplit_sep !== nothing
        itr = binsplit_tab_pairs(itr, binsplit_sep)
    end
    seqs_by_binname = Dict{SubString{String}, Vector{Sequence}}()
    for (binname, seqname) in itr
        (seq, _) = ref.targets_by_name[seqname]
        push!(get!(valtype(seqs_by_binname), seqs_by_binname, binname), seq)
    end
    [Bin(binname, seqs, ref.targets_by_name) for (binname, seqs) in seqs_by_binname]
end

struct ReferenceJSON
    # Genome => (subject => length)
    genomes::Dict{String, Dict{String, Int}}
    # Sequence => (sequence_length, [(subject, from, to)])
    sequences::Dict{String, Tuple{Int, Vector{Tuple{String, Int, Int}}}}
    # child => parent
    taxmaps::Vector{Dict{String, Union{String, Nothing}}}
end
StructTypes.StructType(::Type{ReferenceJSON}) = StructTypes.Struct()

function Reference(io::IO)
    json_struct = JSON3.read(io, ReferenceJSON)
    ref = Reference()

    # Parse genomes
    for (genomename, sourcesdict) in json_struct.genomes
        genome = Genome(genomename)
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
                lazy"\"$(source_by_name[source.name].genome.name)\" and \"$(genome.name)\"."
            )
        end
        source_by_name[source.name] = source
    end

    # Parse sequences
    for (seq_name, (seq_length, targs)) in json_struct.sequences
        targets = map(targs) do (source_name, from, to)
            if to < from
                error(lazy"Sequence \"$(seq_name)\" spans $(from)-$(to), must span at least 1 base.")
            end
            source = get(source_by_name, source_name, nothing)
            if source === nothing
                error(lazy"Sequence \"$(seq_name)\" maps to source \"$(source_name)\", but no such source in reference")
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
    json_dict = Dict{Symbol, Any}(:genomes => Dict(), :sequences => Dict(), :taxmaps => [])
    genome_dict = json_dict[:genomes]
    for genome in ref.genomes
        genome_dict[genome.name] = Dict(s.name => s.length for s in genome.sources)
    end
    sequence_dict = json_dict[:sequences]
    for (_, (sequence, targets)) in ref.targets_by_name
        sequence_dict[sequence.name] = (sequence.length, [(source.name, first(span), last(span)) for (source, span) in targets])
    end
    children = Set{Node}(ref.genomes)
    parents = Set{Node}()
    while !isempty(children)
        d = Dict()
        for child in children
            parent = child.parent
            d[child.name] = (isnothing(parent) ? nothing : parent.name)
            parent === nothing || push!(parents, parent)
        end
        push!(json_dict[:taxmaps], d)
        children, parents = parents, children
        empty!(parents)
    end
    JSON3.write(io, json_dict)
end

function parse_taxonomy(
    genomes::Set{Genome},
    dict::Vector{Dict{String, Union{String, Nothing}}}
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
                "on the previous rank."
            )
            # Create parent if it does not already exist
            parent_name = maybe_parent_name === nothing ? child_name : maybe_parent_name
            parent = get(parent_by_name, parent_name, nothing)::Union{Clade{Genome}, Nothing}
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
            isdefined(child, :parent) || error(
                "At rank $rank, child $(child.name) has no parent"
            )
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
    foreach(i -> sort!(by=j -> j.name, i), result)
    return result
end
