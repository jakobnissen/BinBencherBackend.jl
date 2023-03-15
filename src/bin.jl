struct Bin
    name::String
    sequences::Vector{Sequence}
    # Asmsize: Bases covered by sequences in this bin.
    # Foreign: Sum of sequences not mapping to the given genome
    genomes::Dict{Genome, @NamedTuple{asmsize::Int, foreign::Int}}
    # Recall is max recall of all children
    # Precision is (sum of length of mapping seqs) / breadth
    clades::Dict{Clade{Genome}, @NamedTuple{asm_recall::Float64, genome_recall::Float64, precision::Float64}}
    # Sum of lengths of sequences, cached for efficiency
    breadth::Int
end

function Bin(name::AbstractString, sequences, targets::Dict{String, Tuple{Sequence, Vector{Target}}})
    # To make it work with arbitrary iterables of sequence
    seqs = vec(collect(sequences))::Vector{Sequence}
    breadth = sum(i -> i.length, seqs; init=0)
    # Which sequences map to the given genome, ints in bitset is indices into `seq`.
    genome_mapping = Dict{Genome, BitSet}()
    source_mapping = Dict{Source{Genome}, Vector{UnitRange{Int}}}()
    for (i, seq) in enumerate(seqs), (source, span) in last(targets[seq.name])
        push!(get!(valtype(genome_mapping), genome_mapping, source.genome), i)
        push!(get!(valtype(source_mapping), source_mapping, source), span)
    end
    genomes = Dict{Genome, @NamedTuple{asmsize::Int, foreign::Int}}()
    # Set `foreign`, which we can compute simply by knowing which sequences map to the genomes
    for (genome, set) in genome_mapping
        genomes[genome] = (;asmsize=0, foreign=breadth - sum(i -> seqs[i].length, set; init=0))
    end
    # Incrementally update `asmsize`; we need to compute this on a per-source level
    for (source, spans) in source_mapping
        (asmsize, foreign) = genomes[source.genome]
        genomes[source.genome] = (;asmsize=asmsize + assembly_size(spans), foreign=foreign)
    end
    # We store sets of mapping sequences - the set mapping to a clade is the union of those
    # mapping to its children.
    clade_mapping = Dict{Clade{Genome}, Tuple{Float64, Float64, BitSet}}() # (asm_recall, genome_recall)
    for (genome, (asmsize, _)) in genomes
        asm_recall = asmsize / genome.assembly_size
        genome_recall = asmsize / genome.genome_size
        (old_asm_recall, old_genome_recall, mapping) = get!(() -> (0.0, 0.0, BitSet()), clade_mapping, genome.parent)
        clade_mapping[genome.parent] = (
            max(old_asm_recall, asm_recall),
            max(old_genome_recall, genome_recall),
            union!(mapping, genome_mapping[genome])
        )
    end
    # Now, iteratively compute clades at a higher and higher level.
    # Begin with parents of genome (this generation), then iteratively look at parents and update them
    this_generation = Set((clade for clade in keys(clade_mapping)))
    next_generation = empty(this_generation)
    while true
        while !isempty(this_generation)
            clade = pop!(this_generation)
            parent = clade.parent
            # If top level clade: Do not continue to next generation
            parent === nothing && continue
            (parent_asm_recall, parent_genome_recall, parent_mapping) = get!(() -> (0.0, 0.0, BitSet()), clade_mapping, parent)
            (child_asm_recall, child_genome_recall, child_mapping) = clade_mapping[clade]
            clade_mapping[parent] = (
                max(parent_asm_recall, child_asm_recall),
                max(parent_genome_recall, child_genome_recall),
                union!(parent_mapping, child_mapping)
            )
            push!(next_generation, parent)
        end
        isempty(next_generation) && break
        # Reuse the sets for next generation by swapping them and emptying the now-used up
        (this_generation, next_generation) = (next_generation, this_generation)
        empty!(next_generation)
    end
    # Now, having computed the sets of mapping contigs, we can compute the actual precision values
    clades = Dict{Clade{Genome}, @NamedTuple{asm_recall::Float64, genome_recall::Float64, precision::Float64}}()
    for (clade, (asm_recall, genome_recall, set)) in clade_mapping
        precision = sum(i -> seqs[i].length, set; init=0) / breadth
        clades[clade] = (;asm_recall=asm_recall, genome_recall=genome_recall, precision=precision)
    end
    Bin(String(name), sequences, genomes, clades, breadth)
end

nseqs(x::Bin) = length(x.sequences)

Base.show(io::IO, x::Bin) = print(io, "Bin(", x.name, ')')
function Base.show(io::IO, ::MIME"text/plain", x::Bin)
    if get(io, :compact, false)
        show(io, x)
    else
        print(io,
            "Bin \"", x.name,
            "\"\n  Sequences: ", nseqs(x),
            "\n  Breadth:   ", x.breadth,
            "\n  Intersecting ", length(x.genomes), " genomes"
        )
    end
end

function confusion_matrix(genome::Genome, bin::Bin; assembly::Bool=true)
    (tp, fp) = get(bin.genomes, genome, (0, bin.breadth))
    fn = (assembly ? genome.assembly_size : genome.genome_size) - tp
    (tp, fp, fn)
end

function recall_precision(genome::Genome, bin::Bin; assembly::Bool=true)
    (tp, fp, fn) = confusion_matrix(genome, bin; assembly=assembly)
    recall = tp / (tp + fn)
    precision = tp / (tp + fp)
    (recall, precision)
end

function recall_precision(clade::Clade{Genome}, bin::Bin; assembly::Bool=true)
    (;asm_recall, genome_recall, precision) = get(bin.clades, clade, (;asm_recall=0.0, genome_recall=0.0, precision=0.0))
    assembly ? (asm_recall, precision) : (genome_recall, precision)
end

function fscore(genome::Genome, bin::Bin, b::AbstractFloat)
    recall, precision = recall_precision(genome, bin)
    # Some people say the Fscore is undefined in this case.
    # We define it to be 0.0
    if iszero(recall + precision)
        return 0.0
    end
    (1 + b^2) * (recall*precision) / ((b^2*precision)+recall)
end
f1(genome::Genome, bin::Bin) = fscore(genome, bin, 1.0)