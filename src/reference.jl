struct Reference
    genomes::Set{Genome}
    targets::Dict{Sequence, Target}
    sequence_by_name::Dict{String, Sequence}
    top_clade::Clade{Genome}
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

ngenomes(x::Reference) = length(x.genomes)
nseqs(x::Reference) = length(x.sequence_by_name)
function nranks(x::Reference)
    isempty(x.genomes) && return 0
    x.top_clade.rank + 1
end

function add_genome!(ref::Reference, genome::Genome)
    in(genome, ref.genomes) && error(lazy"Genome $(genome.name) already in reference")
    push!(ref.genomes, genome)
    ref
end

function add_sequence!(ref::Reference, seq::Sequence, target::Target)
    haskey(ref.targets, seq) && error(lazy"Reference already has sequence $(seq.name)")
    if target isa Tuple{Genome, String, UnitRange}
        (genome, subject, span) = target
        source = get(genome.sources, subject, nothing)
        source === nothing && error(
            "Attempted to add sequence $(seq.name) with source $subject " *
            "to genome $(genome.name), but genome has no such source"
        )
        in(genome, ref.genomes) || error(lazy"Genome $(genome.name) not in reference")
        source < last(span) && error(
            "Attempted to add sequence $(seq.name) with span $span " *
            "to subject $subject, but the subject only has length $source"
        )
    end
    ref.targets[seq] = target
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
    [Bin(binname, seqs, ref.targets) for (binname, seqs) in seqs_by_binname]
end

const Pairs{A, B} = Vector{Tuple{A, B}}

 # Imprecise mappings: Sequences mapping to simply a Genome or Clade (i.e Node)
const ImpreciseMappingJSON = Pairs{
    # Node name, node rank
    Tuple{String, Int},
    # Vector of (sequence_name, sequence_length)
    Pairs{String, Int}
}

# Precise mappings. These maps to a specific area in a subject of a Genome
const PreciseMappingJSON = Pairs{
    # Genome name
    String,
    # Subjects
    Pairs{
        # Subject name, length
        Tuple{String, Int},
        # Seqname, seqlength, from, to
        Vector{Tuple{String, Int, Int, Int}}
    }
}
const SequenceJSON = Tuple{ImpreciseMappingJSON, PreciseMappingJSON}

struct ReferenceJSON
    sequences::SequenceJSON
    # A vector of dicts, representing taxonomic levels from Genome and upwards.
    # Each dict is {child_name => parent_name}.
    # The final dict must have only one unique parent name, which is the top clade.
    taxmaps::Vector{Dict{String, String}}
end
StructTypes.StructType(::Type{ReferenceJSON}) = StructTypes.Struct()

# This is currently type unstable (Julia #46557), but its instability is
# shielded by the function barrier of this function.
function n_seqs(x::ReferenceJSON)::Int
    v = sum(values(x.sequences); init=0) do subdict
        Int(sum(length, values(subdict); init=0))::Int
    end
    Int(v)
end

function Reference(io::IO)
    json_struct = JSON3.read(io, ReferenceJSON)
    nseqs = n_seqs(json_struct)
    targets = sizehint!(Dict{Sequence, Target}(), nseqs)

    # Parse taxonomy, update genomes and clades
    (genomes, top_clade) = parse_taxonomy(json_struct.taxmaps)

    clades_by_rank_name = Dict{Tuple{Int, String}, Node}()
    for node in all_children(top_clade)
        clades_by_rank_name[(rank(node), node.name)] = node
    end

    # Parse sequences
    (imprecise, precise) = json_struct.sequences
    for ((clade_name, clade_rank), imprecise_sequences) in imprecise
        clade = clades_by_rank_name[(clade_rank, clade_name)]
        for (sequence_name, sequence_len) in imprecise_sequences
            targets[Sequence(sequence_name, sequence_len)] = clade
        end
    end
    for (genome_name, subject_vector) in precise
        genome = clades_by_rank_name[(0, genome_name)]
        for ((subject_name, subject_len), precise_sequences) in subject_vector
            add_source!(genome, subject_name, subject_len)
            for (sequence_name, sequence_len, map_from, map_to) in precise_sequences
                sequence = Sequence(sequence_name, sequence_len)
                targets[sequence] = (genome, subject_name, map_from:map_to)
            end
        end
    end
    sequence_by_name = Dict(seq.name => seq for seq in keys(targets))
    Reference(genomes, targets, sequence_by_name, top_clade)
end

function save(io::IO, ref::Reference)
    # Make taxmaps
    taxmaps = [Dict{String, String}() for _ in 1:nranks(ref)-1]
    for node in last(Iterators.peel(all_children(ref.top_clade))) # skip top clade
        taxmaps[node.rank + 1][node.name] = node.parent.name
    end

    # Make Dict versions of the SequenceJSON types, which we convert to Vector later.
    # JSON dict keys must be String, so we use Vector in SequenceJSON.
    imprecise = Dict{Tuple{String, Int}, Vector{Tuple{String, Int}}}()
    precise = Dict{String, Dict{Tuple{String, Int}, Vector{Tuple{String, Int, Int, Int}}}}()
    for (sequence, target) in ref.targets
        if target isa Union{Clade, Genome}
            key = (target.name, rank(target))
            push!(get!(valtype(imprecise), imprecise, key), (sequence.name, length(sequence)))
        else
            (genome, subject, span) = target
            subject_dict = get!(valtype(precise), precise, genome.name)
            key = (subject, genome.sources[subject])
            sequences = get!(valtype(subject_dict), subject_dict, key)
            push!(sequences, (sequence.name, length(sequence), first(span), last(span)))
        end
    end
    imprecise_vec = collect(imprecise)
    precise_vec = [k: collect(v) for (k,v) in precise]
    JSON3.write(io, Dict(:taxmaps => taxmaps, :sequences => (imprecise_vec, precise_vec)))
end

function parse_taxonomy(
    dicts::Vector{Dict{String, String}}
)::Tuple{Set{Genome}, Clade{Genome}}
    isempty(dicts) && error("There must be at least one taxmap")
    child_by_name = Dict{String, Node}((name => Genome(name) for name in keys(first(dicts))))
    parent_by_name = empty(child_by_name)
    genomes = Set(values(child_by_name))
    clades = Vector{Clade{Genome}}[]
    for (rank, taxmap) in enumerate(dicts)
        clades_this_rank = Clade{Genome}[]
        for (child_name, maybe_parent_name) in taxmap
            # Get child
            child = get(child_by_name, child_name, nothing)
            child === nothing && error(
                "At rank $rank, found child name \"$child_name\", but this does not exist " *
                "as a parent on the previous rank."
            )

            # Create parent if it does not already exist
            parent_name = maybe_parent_name === nothing ? child_name : maybe_parent_name
            parent = get(parent_by_name, parent_name, nothing)::Union{Clade{Genome}, Nothing}
            if parent === nothing
                parent = Clade(parent_name, child)
                parent_by_name[parent_name] = parent
                push!(clades_this_rank, parent)
            end
            add_child!(parent, child)
            child.parent = parent::Clade{Genome}
        end
        # Enforce all childrens have parents
        for child in values(child_by_name)
            isdefined(child, :parent) || error(
                "At rank $rank, child $(child.name) has no parent"
            )
        end
        parent_by_name, child_by_name = child_by_name, parent_by_name
        empty!(parent_by_name)
        push!(clades, clades_this_rank)
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
            pop!(clades)
        end
    end
    top_clade = only(last(clades))
    return (genomes, top_clade)
end