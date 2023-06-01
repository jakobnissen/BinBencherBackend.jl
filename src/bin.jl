"""
    Bin(name::AbstractString, sequences, targets)

`Bin`s each represent a bin created by the binner. Conceptually, they are simply
a set of `Sequence` with a name attached.
Practically, every `Bin` is benchmarked against all `Genome`s and `Clade`s of a given `Reference`,
so each `Bin` stores data about its intersection with every genome/clade, e.g. its purity and recall.

Like `Source`s, `Bin`s also have an _assembly size_ for a given source. This is the number
of base pairs in the source covered by any sequence in the `Bin`, which is always a subset
of the `Source`'s assembly size.

Benchmark statistics for a `Bin`/`Genome` can be done with either _assemblies
or _genomes_ as the ground truth.
* True positives (TP) are defined as the sum of assembly sizes over all sources in the genome
* False positives (FP) are the sum of length of sequences in the bin not mapping to the `Genome`
* False negatives (FN) is either the genome assembly size or genome size minus TP.

For `Bin`/`Clade` pairs B/C, recall is the maximal recall of B/Ch for all children Ch of C.
Precision is the sum of lengths of sequences mapping to any child of the clade divided
by the sum of lengths of all sequences in the bin.

See also: [`Binning`](@ref), [`Genome`](@ref), [`Clade`](@ref)

# Examples
```jldoctest
julia> bin = first(binning.bins)
Bin "C1"
  Sequences: 2
  Breadth:   65
  Intersecting 1 genome

julia> first(bin.sequences)
Sequence("s1", 25)

julia> f1(first(ref.genomes), bin)
0.625
```
"""
struct Bin
    name::String
    sequences::Vector{Sequence}
    # Asmsize: Bases covered by sequences in this bin.
    # Foreign: Sum of sequences mapping to other sequences, but not this genome
    genomes::Dict{Genome, @NamedTuple{asmsize::Int, foreign::Int}}
    # Recall is max recall of all children
    # Precision is (sum of length of mapping seqs) / breadth
    clades::Dict{
        Clade{Genome},
        @NamedTuple{asm_recall::Float64, genome_recall::Float64, precision::Float64}
    }
    # Sum of lengths of sequences mapping to anything, cached for efficiency
    breadth::Int
end

function Bin(
    name::AbstractString,
    sequences,
    targets::Dict{String, Tuple{Sequence, Vector{Target}}},
)
    # To make it work with arbitrary iterables of sequence
    seqs = vec(collect(sequences))::Vector{Sequence}
    mapping_breadth =
        sum(length(s) for s in seqs if !isempty(last(targets[s.name])); init=0)

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
        genomes[genome] =
            (; asmsize=0, foreign=mapping_breadth - sum(i -> seqs[i].length, set; init=0))
    end
    # Incrementally update `asmsize`; we need to compute this on a per-source level
    for (source, spans) in source_mapping
        (asmsize, foreign) = genomes[source.genome]
        genomes[source.genome] = (; asmsize=asmsize + assembly_size(spans), foreign=foreign)
    end
    # We store sets of mapping sequences - the set mapping to a clade is the union of those
    # mapping to its children.
    clade_mapping = Dict{Clade{Genome}, Tuple{Float64, Float64, BitSet}}() # (asm_recall, genome_recall)
    for (genome, (asmsize, _)) in genomes
        asm_recall = asmsize / genome.assembly_size
        genome_recall = asmsize / genome.genome_size
        (old_asm_recall, old_genome_recall, mapping) =
            get!(() -> (0.0, 0.0, BitSet()), clade_mapping, genome.parent)
        clade_mapping[genome.parent] = (
            max(old_asm_recall, asm_recall),
            max(old_genome_recall, genome_recall),
            union!(mapping, genome_mapping[genome]),
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
            (parent_asm_recall, parent_genome_recall, parent_mapping) =
                get!(() -> (0.0, 0.0, BitSet()), clade_mapping, parent)
            (child_asm_recall, child_genome_recall, child_mapping) = clade_mapping[clade]
            clade_mapping[parent] = (
                max(parent_asm_recall, child_asm_recall),
                max(parent_genome_recall, child_genome_recall),
                union!(parent_mapping, child_mapping),
            )
            push!(next_generation, parent)
        end
        isempty(next_generation) && break
        # Reuse the sets for next generation by swapping them and emptying the now-used up
        (this_generation, next_generation) = (next_generation, this_generation)
        empty!(next_generation)
    end
    # Now, having computed the sets of mapping contigs, we can compute the actual precision values
    clades = Dict{
        Clade{Genome},
        @NamedTuple{asm_recall::Float64, genome_recall::Float64, precision::Float64}
    }()
    for (clade, (asm_recall, genome_recall, set)) in clade_mapping
        precision = sum(i -> seqs[i].length, set; init=0) / mapping_breadth
        clades[clade] =
            (; asm_recall=asm_recall, genome_recall=genome_recall, precision=precision)
    end
    breadth = sum(seqs |> imap(length); init=0)
    Bin(String(name), sequences, genomes, clades, breadth)
end

nseqs(x::Bin) = length(x.sequences)

"""
    intersecting([Genome, Clade]=Genome, x::Bin)

Get an iterator of the `Genome`s or `Clade`s that bin `x` intersects with.
`intersecting(::Bin)` defaults to genomes.

# Example
```jldoctest
julia> collect(intersecting(bin))
1-element Vector{Genome}:
 Genome(gA)

julia> sort!(collect(intersecting(Clade, bin)); by=i -> i.name)
2-element Vector{Clade{Genome}}:
 Species "D", 2 genomes
 Genus "F", 3 genomes
```
"""
intersecting(x::Bin) = intersecting(Genome, x)
intersecting(::Type{Genome}, x::Bin) = keys(x.genomes)
intersecting(::Type{<:Clade}, x::Bin) = keys(x.clades)

Base.show(io::IO, x::Bin) = print(io, "Bin(", x.name, ')')
function Base.show(io::IO, ::MIME"text/plain", x::Bin)
    if get(io, :compact, false)
        show(io, x)
    else
        ngenomes = length(x.genomes)
        suffix = ngenomes > 1 ? "s" : ""
        print(
            io,
            "Bin \"",
            x.name,
            "\"\n  Sequences: ",
            nseqs(x),
            "\n  Breadth:   ",
            x.breadth,
            "\n  Intersecting ",
            length(x.genomes),
            " genome",
            suffix,
        )
    end
end

function confusion_matrix(genome::Genome, bin::Bin; assembly::Bool=true)
    (tp, fp) = get(bin.genomes, genome, (0, bin.breadth))
    fn = (assembly ? genome.assembly_size : genome.genome_size) - tp
    (tp, fp, fn)
end

"""
    recall_precision(x::Union{Genome, Clade}, bin::Bin; assembly::Bool=true)

Get the recall, precision as a 2-tuple of `Float64` for the given genome/bin pair.
See the docstring for `Bin` for how this is computed.

See also: [`Bin`](@ref), [`Binning`](@ref)

# Examples
```jldoctest
julia> bingenome = only(intersecting(bin));

julia> recall_precision(bingenome, bin)
(recall = 0.45454545454545453, precision = 1.0)

julia> recall_precision(bingenome, bin; assembly=false)
(recall = 0.4, precision = 1.0)

julia> recall_precision(bingenome.parent, bin; assembly=false)
(recall = 0.4, precision = 1.0)
```
"""
function recall_precision(genome::Genome, bin::Bin; assembly::Bool=true)
    (tp, fp, fn) = confusion_matrix(genome, bin; assembly=assembly)
    recall = tp / (tp + fn)
    precision = tp / (tp + fp)
    (; recall, precision)
end

function recall_precision(clade::Clade{Genome}, bin::Bin; assembly::Bool=true)
    (; asm_recall, genome_recall, precision) =
        get(bin.clades, clade, (; asm_recall=0.0, genome_recall=0.0, precision=0.0))
    assembly ? (; recall=asm_recall, precision) : (; recall=genome_recall, precision)
end

function fscore(genome::Genome, bin::Bin, b::Real)
    (; recall, precision) = recall_precision(genome, bin)
    # Some people say the Fscore is undefined in this case.
    # We define it to be 0.0
    if iszero(recall + precision)
        return 0.0
    end
    (1 + b^2) * (recall * precision) / ((b^2 * precision) + recall)
end
f1(genome::Genome, bin::Bin) = fscore(genome, bin, 1)
