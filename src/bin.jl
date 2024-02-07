"""
    Bin(name::AbstractString, ref::Reference, sequences)

`Bin`s each represent a bin created by the binner. Conceptually, they are simply
a set of `Sequence` with a name attached.
Practically, every `Bin` is benchmarked against all `Genome`s and `Clade`s of a given `Reference`,
so each `Bin` stores data about its intersection with every genome/clade, e.g. its purity and recall.

Like `Source`s, `Bin`s also have an _assembly size_ for a given `Genome`. This is the number
of base pairs in the genomes covered by any sequence in the `Bin`, which is always a subset
of the genome's assembly size.

Benchmark statistics for a `Bin`/`Genome` can be done with either _assemblies_
or _genomes_ as the ground truth.
* True positives (TP) are defined as the sum of assembly sizes over all sources in the genome
* False positives (FP) are the sum of length of sequences in the bin not mapping to the genome
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

# Note: This constructor is user-facing. We use the more efficient
# `bin_by_indices` when constructing Bins in practise
function Bin(
    name::AbstractString,
    ref::Reference,
    sequences, # iterator of Sequence
    considered_genomes::Union{Nothing, Set{Genome}}=nothing,
)
    indices = [ref.target_index_by_name[s.name] for s in sequences]
    scratch = Vector{Tuple{Int, Int}}()
    bin_by_indices(name, indices, ref.targets, scratch, considered_genomes)
end

function bin_by_indices(
    name::AbstractString,
    seq_indices::Vector{<:Integer},
    targets::Vector{Tuple{Sequence, Vector{Target}}},
    scratch::Vector{Tuple{Int, Int}},
    considered_genomes::Union{Nothing, Set{Genome}},
)
    seqs = [first(targets[i]) for i in seq_indices]

    # Which sequences map to the given genome, ints in bitset is indices into `seqs`.
    genome_mapping = Dict{Genome, BitSet}()
    source_mapping = Dict{Source{Genome}, Vector{Tuple{Int, Int}}}()
    mapping_breadth = 0
    for (i, (seq, idx)) in enumerate(zip(seqs, seq_indices))
        seq_targets = last(targets[idx])
        isempty(seq_targets) || (mapping_breadth += length(seq))
        for (source, span) in seq_targets
            genome = source.genome
            !isnothing(considered_genomes) && !in(genome, considered_genomes) && continue
            push!(get!(valtype(genome_mapping), genome_mapping, genome), i)
            push!(get!(valtype(source_mapping), source_mapping, source), span)
        end
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
        genomes[source.genome] = (;
            asmsize=asmsize + assembly_size!(identity, scratch, spans, source.length),
            foreign=foreign,
        )
    end
    # We store sets of mapping sequences - the set mapping to a clade is the union of those
    # mapping to its children.
    clade_mapping = Dict{
        Clade{Genome},
        @NamedTuple{asm_recall::Float64, genome_recall::Float64, mapping_seqs::BitSet}
    }()
    for (genome, (; asmsize)) in genomes
        asm_recall = asmsize / genome.assembly_size
        genome_recall = asmsize / genome.genome_size
        (old_asm_recall, old_genome_recall, mapping) = get!(
            () -> (; asm_recall=0.0, genome_recall=0.0, mapping_seqs=BitSet()),
            clade_mapping,
            genome.parent,
        )
        clade_mapping[genome.parent] = (;
            asm_recall=max(old_asm_recall, asm_recall),
            genome_recall=max(old_genome_recall, genome_recall),
            mapping_seqs=union!(mapping, genome_mapping[genome]),
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
            (parent_asm_recall, parent_genome_recall, parent_mapping) = get!(
                () -> (; asm_recall=0.0, genome_recall=0.0, mapping_seqs=BitSet()),
                clade_mapping,
                parent,
            )
            (child_asm_recall, child_genome_recall, child_mapping) = clade_mapping[clade]
            clade_mapping[parent] = (
                asm_recall=max(parent_asm_recall, child_asm_recall),
                genome_recall=max(parent_genome_recall, child_genome_recall),
                mapping_seqs=union!(parent_mapping, child_mapping),
            )
            push!(next_generation, parent)
        end
        isempty(next_generation) && break
        # Reuse the sets for next generation by swapping them and emptying the now-used up
        (this_generation, next_generation) = (next_generation, this_generation)
        @assert isempty(next_generation)
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
    Bin(String(name), seqs, genomes, clades, breadth)
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

# NB: Returns NaN if recall and precision is zero
function fscore(recall::Real, precision::Real, b::Real)
    (1 + b^2) * (recall * precision) / ((b^2 * precision) + recall)
end
f1(recall::Real, precision::Real) = fscore(recall, precision, 1)

function fscore(genome::Genome, bin::Bin, b::Real; assembly::Bool=true)
    (; recall, precision) = recall_precision(genome, bin; assembly)
    fscore(recall, precision, b)
end
f1(genome::Genome, bin::Bin; assembly::Bool=true) = fscore(genome, bin, 1; assembly)

function recalls_precisions(::Type{Genome}, bin::Bin; assembly::Bool=true)
    bin.genomes |> imap() do (genome, (; asmsize, foreign))
        fn = (assembly ? genome.assembly_size : genome.genome_size) - asmsize
        recall = asmsize / (asmsize + fn)
        precision = asmsize / (asmsize + foreign)
        (; genome, recall, precision)
    end
end

function recalls_precisions(::Type{<:Clade}, bin::Bin; assembly::Bool=true)
    bin.clades |> imap() do (clade, (; asm_recall, genome_recall, precision))
        recall = assembly ? asm_recall : genome_recall
        (; clade, recall, precision)
    end
end

function recalls_precisions(bin::Bin; assembly::Bool=true)
    recalls_precisions(Genome, bin; assembly)
end

# TODO: Why does this function allocate?
# Compute recall/precision of the genome with highest F1 for this bin
function recall_prec_max_f1(bin::Bin; assembly::Bool=true)
    (max_recall, max_precision, max_f1) = (0.0, 0.0, 0.0)
    for (; recall, precision) in recalls_precisions(bin; assembly)
        this_f1 = f1(recall, precision)
        if this_f1 > max_f1
            (max_recall, max_precision, max_f1) = (recall, precision, this_f1)
        end
    end
    # This can happen if the bin only contain sequences unassigned to any genome
    # in which case recalls_precisions returns an iterable with zero elements
    iszero(max_f1) ? nothing : (; recall=max_recall, precision=max_precision)
end

"""
    passes_f1(bin::Bin, threshold::Real; assembly::Bool=false)::Bool

Computes if `bin` has an F1 score equal to, or higher than `threshold` for any genome.

# Examples
```jldoctest
julia> obs_f1 = f1(only(intersecting(bin)), bin)
0.625

julia> passes_f1(bin, obs_f1)
true

julia> passes_f1(bin, obs_f1 + 0.001)
false
```
"""
function passes_f1(bin::Bin, threshold::Real; assembly::Bool=true)::Bool
    any(recalls_precisions(Genome, bin; assembly)) do (; recall, precision)
        f1(recall, precision) ≥ threshold
    end
end

"""
    passes_recall_precision(bin::Bin, recall::Real, precision::Real; assembly::Bool=false)::Bool

Computes if `bin` intersects with any `Genome` with at least the given recall and precision thresholds.

# Examples
```jldoctest
julia> (r, p) = recall_precision(only(intersecting(bin)), bin)
(recall = 0.45454545454545453, precision = 1.0)

julia> passes_recall_precision(bin, 0.45, 1.0)
true

julia> passes_recall_precision(bin, 0.46, 1.0)
false
```
"""
function passes_recall_precision(bin::Bin, rec::Real, prec::Real; assembly::Bool=true)::Bool
    any(recalls_precisions(Genome, bin; assembly)) do (; recall, precision)
        recall ≥ rec && precision ≥ prec
    end
end
