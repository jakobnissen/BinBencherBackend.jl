struct Reference
    genomes::Set{Genome}
    genomeof::Dict{Sequence, Genome}
    sequence_by_name::Dict{String, Sequence}
    clades::Vector{Vector{Clade{Genome}}}
end

function Reference()
    Reference(
        Set{Genome}(),
        Dict{Sequence, Genome}(),
        Dict{String, Sequence}(),
        Vector{Clade{Genome}}[]
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
        )
    end
end

top_clade(x::Reference) = only(last(x.clades))

ngenomes(x::Reference) = length(x.genomes)
nseqs(x::Reference) = length(x.genomeof)
function nranks(x::Reference)
    isempty(x.genomes) && return 0
    length(x.clades) + 1
end

function add_genome!(ref::Reference, genome::Genome)
    in(genome, ref.genomes) && error("Genome $(genome.name) already in reference")
    push!(ref.genomes, genome)
    ref
end

function add_sequence!(ref::Reference, seq::Sequence, genome::Genome)
    source = get(genome.sources, seq.subject, nothing)
    source === nothing && error(
        "Attempted to add sequence $(seq.name) with source $(seq.subject) " *
        "to genome $(genome.name), but genome has no such source"
    )
    in(genome, ref.genomes) || error("Genome $(genome.name) not in reference")
    source < last(seq.span) && error(
        "Attempted to add sequence $(seq.name) with span $(first(seq.span)):$(last(seq.span)) " *
        "to subject $(seq.subject), but the subject only has length $source"
    )
    haskey(ref.genomeof, seq) && error("Reference already has sequence $(seq.name)")
    ref.genomeof[seq] = genome
    ref.sequence_by_name[seq.name] = seq
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
        seq = ref.sequence_by_name[seqname]
        push!(get!(valtype(seqs_by_binname), seqs_by_binname, binname), seq)
    end
    [Bin(binname, seqs, ref.genomeof) for (binname, seqs) in seqs_by_binname]
end

struct ReferenceJSON
    genomes::Dict{String, Dict{String, Tuple{Int, Dict{String, Tuple{Int, Int}}}}}
    taxmaps::Vector{Dict{String, Union{String, Nothing}}}
end
StructTypes.StructType(::Type{ReferenceJSON}) = StructTypes.Struct()

function Reference(io::IO)
    json_struct = JSON3.read(io, ReferenceJSON)
    ref = Reference()

    nseqs = sum(values(json_struct.genomes)) do v1
        sum(i -> length(i[2]), values(v1))
    end
    sizehint!(ref.genomeof, nseqs)
    sizehint!(ref.sequence_by_name, nseqs)

    # Parse genomes
    for (genomename, sourcesdict) in json_struct.genomes
        genome = Genome(genomename)
        add_genome!(ref, genome)
        for (sourcename, (sourcelen, contigdict)) in sourcesdict
            add_source!(genome, sourcename, sourcelen)
            for (seqname, (start, stop)) in contigdict
                seq = Sequence(seqname, sourcename, start:stop)
                add_sequence!(ref, seq, genome)
            end
        end
    end
    copy!(ref.clades, parse_taxonomy(ref.genomes, json_struct.taxmaps))
    ref
end

function save(io::IO, ref::Reference)
    json_dict = Dict{Symbol, Any}(:genomes => Dict(), :taxmaps => [])
    genome_dict = json_dict[:genomes]
    for genome in ref.genomes
        source_dict = Dict()
        genome_dict[genome.name] = source_dict
        for (sourcename, len) in genome.sources
            source_dict[sourcename] = (len, Dict())
        end
    end
    for (seq, genome) in ref.genomeof
        genome_dict[genome.name][seq.subject][2][seq.name] = extrema(seq.span)
    end
    children = Set{Node}(ref.genomes)
    parents = Set{Node}()
    while !isempty(children)
        d = Dict()
        for child in children
            parent = child.parent
            d[child] = parent
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
    return result
end