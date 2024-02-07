const DEFAULT_RECALLS = (0.6, 0.7, 0.8, 0.9, 0.95, 0.99)
const DEFAULT_PRECISIONS = (0.6, 0.7, 0.8, 0.9, 0.95, 0.99)

struct BinStats
    # Recall and precision are computed by pairing each bin
    # with the Genome that yields the highest F1 score
    mean_bin_recall::Float64
    mean_bin_precision::Float64
    mean_bin_f1::Float64
end

function BinStats(bins::Vector{Bin}; assembly::Bool=true)
    (mean_bin_recall, mean_bin_precision, mean_bin_f1) =
        mean_bin_recall_prec(bins; assembly)
    BinStats(mean_bin_recall, mean_bin_precision, mean_bin_f1)
end

function mean_bin_recall_prec(bins::Vector{Bin}; assembly::Bool=true)
    (recall_sum, prec_sum, f1_sum, nbins) = foldl(
        Iterators.filter(
            !isnothing,
            Iterators.map(b -> recall_prec_max_f1(b; assembly), bins),
        );
        init=(0.0, 0.0, 0.0, 0),
    ) do (recall_sum, prec_sum, f1_sum, nbins), (; recall, precision)
        (recall_sum + recall, prec_sum + precision, f1_sum + f1(recall, precision), nbins + 1)
    end
    (recall_sum / nbins, prec_sum / nbins, f1_sum / nbins)
end

"""
    Binning(::Union{IO, AbstractString}, ::Reference; kwargs...)

A `Binning` represents a set of `Bin`s benchmarked against a `Reference`.
`Binning`s can be created given a set of `Bin`s and a `Reference`, where the
bins may potentially be loaded from a `.tsv` file.
The fields `recovered_asms` and `recovered_genomes` are used for benchmarking,
these are normally output using the `print_matrix` function.

A `Binning` is loaded from a tsv file, which is specified either as an `IO`,
or its path as an `AbstractString`. If the path ends with `.gz`, automatically
gzip decompress when reading the file.

See also: [`print_matrix`](@ref), [`Bin`](@ref), [`Reference`](@ref)

# Examples
```jldoctest
julia> bins = Binning(path_to_bins_file, ref);

julia> bins isa Binning
true

julia> BinBencherBackend.n_nc(binning)
0
```

# Extended help
Create with:
```julia
open(file) do io
    Binning(
        io::Union{IO, AbstractString},
        ref::Reference;
        min_size::Integer=1,
        min_seqs::Integer=1,
        binsplit_separator::Union{AbstractString, Char, Nothing}=nothing,
        disjoint::Bool=true,
        recalls=DEFAULT_RECALLS,
        precisions=DEFAULT_PRECISIONS,
        filter_genomes=Returns(true)
)
```

* `min_size`: Filter away bins with breadth lower than this
* `min_seqs`: Filter away bins with fewer sequences that this
* `binsplit_separator`: Split bins based on this separator
  (`nothing` means no binsplitting)
* `disjoint`: Throw an error if the same sequence is seen in multiple bins
* `recalls` and `precision`: The thresholds to benchmark with
* `filter_genomes`: A function `f(genome)::Bool`. Genomes for which it returns
   `false` are ignored in benchmarking.
"""
struct Binning
    ref::Reference
    bins::Vector{Bin}
    # One matrix per rank - first is genome, then upwards
    recovered_asms::Vector{Matrix{Int}}
    recovered_genomes::Vector{Matrix{Int}}
    recalls::Vector{Float64}
    precisions::Vector{Float64}
    bin_asm_stats::BinStats
    bin_genome_stats::BinStats
end

function Base.show(io::IO, x::Binning)
    nc = n_nc(x)
    print(io, summary(x), '(')
    if nc !== nothing
        print(io, "NC = ", nc)
    end
    print(io, ')')
end

function Base.show(io::IO, ::MIME"text/plain", x::Binning)
    if get(io, :compact, false)
        show(io, x)
    else
        buf = IOBuffer()
        show(buf, MIME"text/plain"(), x.ref)
        seekstart(buf)
        println(io, "Binning")
        for line in eachline(buf)
            println(io, "  ", line)
        end
        print(io, "  Bins:        ", nbins(x))
        nc = n_nc(x)
        if nc !== nothing
            print(io, "\n  NC genomes:  ", nc)
        end
        for (stats, name) in
            [(x.bin_genome_stats, "genome  "), (x.bin_asm_stats, "assembly")]
            print(
                io,
                "\n  Mean bin ",
                name,
                " R/P/F1: ",
                round(stats.mean_bin_recall; digits=3),
                " / ",
                round(stats.mean_bin_precision; digits=3),
                " / ",
                round(stats.mean_bin_f1; digits=3),
            )
        end
        print(io, "\n  Precisions: ", repr([round(i; digits=3) for i in x.precisions]))
        print(io, "\n  Recalls:    ", repr([round(i; digits=3) for i in x.recalls]))
        print(io, "\n  Reconstruction (assemblies):\n")
        seekstart(buf)
        print_matrix(buf, x; level=0, assembly=true)
        seekstart(buf)
        for line in eachline(buf)
            println(io, "    ", line)
        end
    end
end

"""
    n_recovered(::Binning, recall, precision; level=0, assembly=false)::Integer

Return the number of genomes or clades reconstructed in the `Binning` at the given recall and precision levels.
If `assembly` is set, return the number of assemblies reconstructed instead.
The argument `level` sets the taxonomic rank: 0 for `Genome` (or assemblies).

# Examples
```jldoctest
julia> n_recovered(binning, 0.4, 0.71)
1

julia> n_recovered(binning, 0.4, 0.71; assembly=true)
2

julia> n_recovered(binning, 0.4, 0.71; assembly=true, level=2)
1
```
"""
function n_recovered(
    binning::Binning,
    recall::Real,
    precision::Real;
    level::Integer=0,
    assembly::Bool=false,
)::Integer
    ri = searchsortedfirst(binning.recalls, recall)
    ri > length(binning.recalls) && error("Binning did not benchmark at that high recall")
    pi = searchsortedfirst(binning.precisions, precision)
    pi > length(binning.precisions) &&
        error("Binning did not benchmark at that high precision")
    matrices = assembly ? binning.recovered_asms : binning.recovered_genomes
    if level + 1 ∉ eachindex(matrices)
        error(
            lazy"Requested bins at taxonomic level $level but have only level 0:$(lastindex(matrices)-1)",
        )
    end
    m = matrices[level + 1]
    m[pi, ri]
end

"""
    print_matrix(::Binning; level=0, assembly=true)

Print the number of reconstructed assemblies or genomes at the given taxonomic level (rank).
Level 0 corresponds to genomes, level 1 to species, etc.
If `assembly`, print the number of reconstructed assemblies, else print the level
of reconstructed genomes.

See also: [`Binning`](@ref)
"""
print_matrix(x::Binning; kwargs...) = print_matrix(stdout, x; kwargs...)
function print_matrix(io::IO, x::Binning; level::Integer=0, assembly::Bool=true)
    ms = assembly ? x.recovered_asms : x.recovered_genomes
    m = ms[level + 1]
    rnd(x) = string(round(x; digits=3))
    digitwidth(x) = sizeof(rnd(x))
    width = max(
        max(4, ndigits(maximum(m; init=0) + 1)),
        maximum(digitwidth, x.recalls; init=0) + 1,
    )
    col1_width = max(3, maximum(digitwidth, x.precisions))
    println(io, rpad("P\\R", col1_width), join([lpad(i, width) for i in x.recalls]))
    for (prec_index, prec) in enumerate(x.precisions)
        println(
            io,
            rpad(rnd(prec), col1_width),
            join([lpad(i, width) for i in m[prec_index, :]]),
        )
    end
end

function n_nc(x::Binning)::Union{Int, Nothing}
    rec_i = findfirst(isequal(0.90), x.recalls)
    rec_i === nothing && return nothing
    prec_i = findfirst(isequal(0.95), x.precisions)
    prec_i === nothing && return nothing
    x.recovered_genomes[1][prec_i, rec_i]
end

nbins(x::Binning) = length(x.bins)

function Binning(path::AbstractString, ref::Reference; kwargs...)
    open_perhaps_gzipped(io -> Binning(io, ref; kwargs...), String(path))
end

# This is the most common constructor in practice, because we usually load binnings from files,
# and also the most efficient. I immediately convert the sequences to integers for faster
# processing, then convert back to seqs before instantiating the bins.
function Binning(
    io::IO,
    ref::Reference;
    min_size::Integer=1,
    min_seqs::Integer=1,
    binsplit_separator::Union{AbstractString, Char, Nothing}=nothing,
    disjoint::Bool=true,
    recalls=DEFAULT_RECALLS,
    precisions=DEFAULT_PRECISIONS,
    filter_genomes::Function=Returns(true),
)
    idxs_by_binname = parse_bins(io, Dict, ref, binsplit_separator, disjoint)
    filter!(idxs_by_binname) do (_, idxs)
        length(idxs) ≥ min_seqs &&
            sum((length(first(ref.targets[i])) for i in idxs); init=0) ≥ min_size
    end
    scratch = Tuple{Int, Int}[]
    considered_genomes = if filter_genomes === Returns(true)
        nothing
    else
        Set(g for g in genomes(ref) if filter_genomes(g))
    end
    bins = [
        bin_by_indices(binname, seq_idxs, ref.targets, scratch, considered_genomes) for
        (binname, seq_idxs) in idxs_by_binname
    ]
    sort!(bins; by=i -> i.name)
    # We already checked for disjointedness when parsing bins, so we skip it here
    Binning(bins, ref; recalls, precisions, disjoint=false)
end

function Binning(
    bins_,
    ref::Reference;
    recalls=DEFAULT_RECALLS,
    precisions=DEFAULT_PRECISIONS,
    disjoint::Bool=true,
)
    checked_recalls = validate_recall_precision(recalls)
    checked_precisions = validate_recall_precision(precisions)
    bins = vector(bins_)
    disjoint && check_disjoint(bins)
    bin_asm_stats = BinStats(bins; assembly=true)
    bin_genome_stats = BinStats(bins; assembly=false)
    (asm_matrices, genome_matrices) =
        benchmark(ref, bins, checked_recalls, checked_precisions)
    Binning(
        ref,
        bins,
        asm_matrices,
        genome_matrices,
        checked_recalls,
        checked_precisions,
        bin_asm_stats,
        bin_genome_stats,
    )
end

"""
    gold_standard(
        ref::Reference
        [sequences, an iterable of bins or Binning];
        disjoint=true,
        recalls=DEFAULT_RECALLS,
        precisions=DEFAULT_PRECISIONS
    )::Binning

Create the optimal `Binning` object given a `Reference`, by the optimal binning of
`sequences`.
If `disjoint`, assign each sequence to only a single genome.

If `sequences` is not passed, use all sequences in `ref`. If a `Binning` is passed,
use all sequences in any of its bins. Else, pass an iterable of `Sequence`.
The elements of `Sequence` must be unique.

## Extended help
Currently, the `disjoint` option uses a simple greedy algorithm to assign
sequences to genomes.
"""
function gold_standard(
    ref::Reference;
    disjoint::Bool=true,
    recalls=DEFAULT_RECALLS,
    precisions=DEFAULT_PRECISIONS,
)::Binning
    gold_standard(ref, (first(v) for v in ref.targets); disjoint, recalls, precisions)
end

function gold_standard(
    ref::Reference,
    binning::Binning;
    disjoint::Bool=true,
    recalls=DEFAULT_RECALLS,
    precisions=DEFAULT_PRECISIONS,
)::Binning
    seqs = reduce((s, bin) -> union!(s, bin.sequences), binning.bins; init=Set{Sequence}())
    gold_standard(ref, seqs; disjoint, recalls, precisions)
end

function gold_standard(
    ref::Reference,
    sequences; # Iterator of Sequence
    disjoint::Bool=true,
    recalls=DEFAULT_RECALLS,
    precisions=DEFAULT_PRECISIONS,
)::Binning
    sequences_of_genome = Dict{Genome, Set{Sequence}}()
    for sequence::Sequence in sequences
        targets = last(ref.targets[ref.target_index_by_name[sequence.name]])
        isempty(targets) && continue
        best_index = disjoint ? last(findmax(i -> length(last(i)), targets)) : 0
        for (i, (source, _)) in enumerate(targets)
            disjoint && i != best_index && continue
            push!(
                get!(valtype(sequences_of_genome), sequences_of_genome, source.genome),
                sequence,
            )
        end
    end
    scratch = Tuple{Int, Int}[]
    bins = [
        bin_by_indices(
            "bin_" * genome.name,
            [ref.target_index_by_name[i.name] for i in seqs],
            ref.targets,
            scratch,
            nothing,
        ) for (genome, seqs) in sequences_of_genome
    ]
    sort!(bins; by=i -> i.name)
    Binning(bins, ref; recalls=recalls, precisions=precisions, disjoint=false)
end

function check_disjoint(bins)
    nseq = sum(i -> length(i.sequences), bins; init=0)
    seen_seqs = sizehint!(Set{Sequence}(), nseq)
    for bin in bins, seq in bin.sequences
        in!(seen_seqs, seq) &&
            error(lazy"Sequence \"$(seq.name)\" seen twice in disjoint Binning")
    end
    nothing
end

function validate_recall_precision(xs)::Vector{Float64}
    s = Set{Float64}()
    for x_ in xs
        x = Float64(x_)
        x in s && error(lazy"Recall/precision value $x present multiple times")
        if !isfinite(x) || x <= 0.0 || x > 1.0
            error(lazy"Recall precision value $x is not finite in (0,1]")
        end
        push!(s, x)
    end
    isempty(s) && error("Must provide at least 1 recall/precision value")
    sort!(collect(s))
end

function benchmark(
    ref::Reference,
    bins::Vector{Bin},
    recalls::Vector{Float64},
    precisions::Vector{Float64};
)::NTuple{2, Vector{<:Matrix{<:Integer}}}
    # For each genome/clade, we compute the maximal recall at the given precision levels.
    # i.e. if 3rd element of vector is 0.5, it means that at precision precisions[3], this genome/clade
    # is found with a maximal recall of 0.5.
    # We keep two vectors: First for recovered assemblies, second for recovered genomes
    max_genome_recall_at_precision = Dict{Genome, Tuple{Vector{Float64}, Vector{Float64}}}()
    max_clade_recall_at_precision =
        Dict{Clade{Genome}, Tuple{Vector{Float64}, Vector{Float64}}}()
    # Initialize with zeros for all known genomes/clades
    for genome in ref.genomes
        max_genome_recall_at_precision[genome] =
            (zeros(Float64, length(precisions)), zeros(Float64, length(precisions)))
    end
    for level in ref.clades, clade in level
        max_clade_recall_at_precision[clade] =
            (zeros(Float64, length(precisions)), zeros(Float64, length(precisions)))
    end
    for bin in bins
        for (genome, (asmsize, foreign)) in bin.genomes
            (v_asm, v_genome) = max_genome_recall_at_precision[genome]
            precision = asmsize / (asmsize + foreign)
            asm_recall = asmsize / genome.assembly_size
            genome_recall = asmsize / genome.genome_size
            # Find the index corresponding to the given precision. If the precision is lower than
            # the smallest precision in `precisions`, don't do anything
            precision_index = searchsortedlast(precisions, precision)
            iszero(precision_index) && continue
            v_asm[precision_index] = max(asm_recall, v_asm[precision_index])
            v_genome[precision_index] = max(genome_recall, v_genome[precision_index])
        end
        # Same as above, but for clades
        for (clade, (asm_recall, genome_recall, precision)) in bin.clades
            precision_index = searchsortedlast(precisions, precision)
            iszero(precision_index) && continue
            (v_asm, v_genome) = max_clade_recall_at_precision[clade]
            for (v, recall) in ((v_asm, asm_recall), (v_genome, genome_recall))
                v[precision_index] = max(recall, v[precision_index])
            end
        end
    end
    # Make all the vectors cumulative. I.e. if a genome is seen at precisions[3] with recall=0.6,
    # then it's also seen at precisions[2] and precisions[1] with at least recall=0.6.
    for (v1, v2) in values(max_genome_recall_at_precision)
        make_reverse_cumulative!(v1)
        make_reverse_cumulative!(v2)
    end
    for (v1, v2) in values(max_clade_recall_at_precision)
        make_reverse_cumulative!(v1)
        make_reverse_cumulative!(v2)
    end
    # Now make the matrices. If at precision `precisions[3]`, we find a recall of 0.6,
    # then we find the index in recalls that correspond to 0.6, say 4. Then, add 1 to all entries
    # in the matrix[1:4, 3].
    asm_matrices = [zeros(Int, length(recalls), length(precisions)) for i in 1:nranks(ref)]
    genome_matrices = [copy(i) for i in asm_matrices]
    for (v_asm, v_genome) in values(max_genome_recall_at_precision)
        for (v, m) in ((v_asm, asm_matrices[1]), (v_genome, genome_matrices[1]))
            update_matrix!(m, v, recalls)
        end
    end
    for (clade, (v_asm, v_genome)) in max_clade_recall_at_precision
        for (v, matrices) in ((v_asm, asm_matrices), (v_genome, genome_matrices))
            update_matrix!(matrices[clade.rank + 1], v, recalls)
        end
    end
    (map(permutedims, asm_matrices), map(permutedims, genome_matrices))
end

function update_matrix!(
    matrix::Matrix{<:Integer},
    v::Vector{<:AbstractFloat},
    recalls::Vector{Float64},
)
    # Since we iterate over increasing precisions, the recall_index must shrink
    # or stay the same per iteration. So, we can reduce the search space
    # at each iteration by modifying imax
    imax = lastindex(recalls)
    for (precision_index, recall) in enumerate(v)
        recall_index = searchsortedlast(view(recalls, 1:imax), recall)
        imax = min(recall_index, imax)
        matrix[1:recall_index, precision_index] .+= 1
    end
    matrix
end

function make_reverse_cumulative!(v::Vector{<:Real})
    @inbounds for i in (length(v) - 1):-1:1
        v[i] = max(v[i], v[i + 1])
    end
    v
end
