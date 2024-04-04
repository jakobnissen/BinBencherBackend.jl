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
    # TODO: Should these be a dense 3D tensor instead of a vector of matrices?
    recovered_asms::Vector{Matrix{Int}}
    recovered_genomes::Vector{Matrix{Int}}
    bin_asms::Vector{Matrix{Int}}
    bin_genomes::Vector{Matrix{Int}}
    recalls::Vector{Float64}
    precisions::Vector{Float64}
    bin_asm_stats::BinStats
    bin_genome_stats::BinStats
end

function Base.show(io::IO, x::Binning)
    nc = n_nc(x)
    print(io, summary(x), '(')
    if nc isa Integer
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
        if nc isa Integer
            print(io, "\n  NC genomes:  ", nc)
        end
        npb = n_passing_bins(x, 0.9, 0.95)
        if npb isa Integer
            print(io, "\n  HQ bins:     ", npb)
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

struct RecallTooHigh end
struct PrecisionTooHigh end
struct RankOutOfRange end
struct ThresholdError
    x::Union{RecallTooHigh, PrecisionTooHigh, RankOutOfRange}
end

"""
    n_recovered(::Binning, recall, precision; level=0, assembly=false)::Integer

Return the number of genomes or clades reconstructed in the `Binning` at the given recall
and precision levels.
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
)::Union{Integer, ThresholdError}
    matrices = assembly ? binning.recovered_asms : binning.recovered_genomes
    extract_from_matrix(binning, recall, precision, matrices, level)
end

"""
    n_passing_bins(::Binning, recall, precision; level=0, assembly::Bool=false)::Integer

Return the number of bins which correspond to any genome or clade at the given recall
and precision levels.
If `assembly` is set, a recall of 1.0 means a bin corresponds to a whole assembly,
else it corresponds to a whole genome.
The argument `level` sets the taxonomic rank: 0 for `Genome` (or assemblies).

# Examples
```jldoctest
julia> n_passing_bins(binning, 0.4, 0.71)
1

julia> n_passing_bins(binning, 0.65, 0.71)
0
```
"""
function n_passing_bins(
    binning::Binning,
    recall::Real,
    precision::Real;
    level::Integer=0,
    assembly::Bool=false,
)::Union{Integer, ThresholdError}
    matrices = assembly ? binning.bin_asms : binning.bin_genomes
    extract_from_matrix(binning, recall, precision, matrices, level)
end

function extract_from_matrix(
    binning::Binning,
    recall::Real,
    precision::Real,
    matrices::Vector{<:Matrix},
    level::Integer=0,
)::Union{Integer, ThresholdError}
    inds = recall_precision_indices(binning, recall, precision)
    inds isa ThresholdError && return inds
    (ri, pi) = inds
    if level + 1 ∉ eachindex(matrices)
        return ThresholdError(RankOutOfRange())
    end
    m = matrices[level + 1]
    m[pi, ri]
end

function recall_precision_indices(
    binning::Binning,
    recall::Real,
    precision::Real,
)::Union{Tuple{Int, Int}, ThresholdError}
    ri = searchsortedfirst(binning.recalls, recall)
    ri > length(binning.recalls) && return ThresholdError(RecallTooHigh())
    pi = searchsortedfirst(binning.precisions, precision)
    pi > length(binning.precisions) && return ThresholdError(PrecisionTooHigh())
    (ri, pi)
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

function n_nc(x::Binning)::Union{Int, ThresholdError}
    extract_from_matrix(x, 0.9, 0.95, x.recovered_genomes, 0)
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
    (asm_matrices, genome_matrices, bin_asms, bin_genomes) =
        benchmark(ref, bins, checked_recalls, checked_precisions)
    Binning(
        ref,
        bins,
        asm_matrices,
        genome_matrices,
        bin_asms,
        bin_genomes,
        checked_recalls,
        checked_precisions,
        bin_asm_stats,
        bin_genome_stats,
    )
end

"""
    gold_standard(
        ref::Reference
        [sequences, a Binning or an iterable of Sequence];
        disjoint=true,
        recalls=DEFAULT_RECALLS,
        precisions=DEFAULT_PRECISIONS
    )::Binning

Create the optimal `Binning` object given a `Reference`, by the optimal binning of
the `Sequence`s in `sequences`.
If `disjoint`, assign each sequence to only a single genome.

If `sequences` is not passed, use all sequences in `ref`. If a `Binning` is passed,
use all sequences in any of its bins. Else, pass an iterable of `Sequence`.

# Extended help
Currently, the `disjoint` option uses a simple greedy algorithm to assign
sequences to genomes.
"""
function gold_standard(ref::Reference; kwargs...)::Binning
    gold_standard(ref, Set(first(v) for v in ref.targets)::Set{Sequence}; kwargs...)
end

function gold_standard(ref::Reference, binning::Binning; kwargs...)::Binning
    seqs::Set{Sequence} = reduce(binning.bins; init=Set{Sequence}()) do s, bin
        union!(s, bin.sequences)
    end
    gold_standard(ref, seqs; kwargs...)
end

function gold_standard(ref::Reference, sequences; kwargs...)::Binning
    gold_standard(ref, Set(sequences)::Set{Sequence}; kwargs...)
end

function gold_standard(
    ref::Reference,
    sequences::Set{Sequence};
    disjoint::Bool=true,
    recalls=DEFAULT_RECALLS,
    precisions=DEFAULT_PRECISIONS,
)::Binning
    sequences_of_genome = Dict{Genome, Set{Sequence}}()
    for sequence in sequences
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

# TODO: This function is massive. Can I refactor it into smaller functions?
function benchmark(
    ref::Reference,
    bins::Vector{Bin},
    recalls::Vector{Float64},
    precisions::Vector{Float64};
)::NTuple{4, Vector{<:Matrix{<:Integer}}}
    # We have 8 qualitatively different combinations of the three options below:
    # * 1) rank 0 (a Genome), versus 2) another rank (a Clade)
    # * Counting 1) unique genomes/clades recovered at a r/p threshold level, or
    #   2) bins corresponding to any genome/clade at a r/p level
    # * With recall=1.0 defined as 1) The assembly size, or 2) the genome size

    # We loop over bins below. So, to tally up recovered genomes/clades, we need to store
    # them in a dict.
    # Given N distinct precision values, we compute the maximally seen recall value at that
    # given precision value. We could also have done it the other way.
    max_genome_recall_at_precision = Dict{Genome, Tuple{Vector{Float64}, Vector{Float64}}}()
    max_clade_recall_at_precision =
        Dict{Clade{Genome}, Tuple{Vector{Float64}, Vector{Float64}}}()
    for genome in ref.genomes
        max_genome_recall_at_precision[genome] =
            (zeros(Float64, length(precisions)), zeros(Float64, length(precisions)))
    end
    for level in ref.clades, clade in level
        max_clade_recall_at_precision[clade] =
            (zeros(Float64, length(precisions)), zeros(Float64, length(precisions)))
    end

    # Same as above, except that because we loop over bins, we can compute at each iteration
    # whether the bin passes the threshold. So we don't need to store in a dict, but can
    # increment the matrices directly at the end of the loop
    bin_max_genome_recall_at_precision =
        [Vector{Float64}(undef, length(precisions)) for _ in 1:nranks(ref)]
    bin_max_asm_recall_at_precision =
        [Vector{Float64}(undef, length(precisions)) for _ in 1:nranks(ref)]
    bin_genome_matrices = [
        zeros(Int, length(recalls), length(precisions)) for
        _ in bin_max_genome_recall_at_precision
    ]
    bin_asm_matrices = [
        zeros(Int, length(recalls), length(precisions)) for
        _ in bin_max_genome_recall_at_precision
    ]

    for bin in bins
        # Since we re-use these vectors for each bin, we need to zero them out between each iteration
        for vs in (bin_max_genome_recall_at_precision, bin_max_asm_recall_at_precision)
            for v in vs
                fill!(v, zero(eltype(v)))
            end
        end

        # First, handle the genomes (as opposed to the clades)
        bin_genome_recalls_1 = bin_max_genome_recall_at_precision[1]
        bin_asm_recalls_1 = bin_max_asm_recall_at_precision[1]
        for (genome, (; asmsize, foreign)) in bin.genomes
            (v_asm, v_genome) = max_genome_recall_at_precision[genome]
            precision = asmsize / (asmsize + foreign)
            asm_recall = asmsize / genome.assembly_size
            genome_recall = asmsize / genome.genome_size

            # Find the index corresponding to the given precision. If the precision is lower than
            # the smallest precision in `precisions`, don't do anything
            precision_index = searchsortedlast(precisions, precision)
            iszero(precision_index) && continue

            # Get maximum recalls for genomes
            v_asm[precision_index] = max(asm_recall, v_asm[precision_index])
            v_genome[precision_index] = max(genome_recall, v_genome[precision_index])

            # Get maximum recalls for bins
            bin_genome_recalls_1[precision_index] =
                max(bin_genome_recalls_1[precision_index], genome_recall)
            bin_asm_recalls_1[precision_index] =
                max(bin_asm_recalls_1[precision_index], asm_recall)
        end

        # Same as above, but for clades instead of genomes
        for (clade, (asm_recall, clade_recall, precision)) in bin.clades
            precision_index = searchsortedlast(precisions, precision)
            iszero(precision_index) && continue
            (v_asm, v_genome) = max_clade_recall_at_precision[clade]

            # Get maximum recall for clades
            for (v, recall) in ((v_asm, asm_recall), (v_genome, clade_recall))
                v[precision_index] = max(recall, v[precision_index])
            end

            # Get maximum recall for bins
            for (v, recall) in (
                (bin_max_asm_recall_at_precision, asm_recall),
                (bin_max_genome_recall_at_precision, clade_recall),
            )
                vr = v[clade.rank + 1]
                vr[precision_index] = max(vr[precision_index], recall)
            end
        end

        # Now that the maximal recall at given precisions for this bin has been filled out,
        # we use the data to increment the correct row in the matrix.
        for (vs, ms) in (
            (bin_max_genome_recall_at_precision, bin_genome_matrices),
            (bin_max_asm_recall_at_precision, bin_asm_matrices),
        )
            for (v, m) in zip(vs, ms)
                # First, we make the recall vector reverse cumulative.
                # If a bin was seen with recall 0.6 at precision 0.8, then it's also
                # seen with at least recall 0.6 at every lower precision level
                make_reverse_max!(v)

                # Now, increment the correct (recall) row for each (precision) column.
                # E.g. if a bin is seen at precision 0.5 with a recall of 0.6, we increment
                # the corresponding rows/columns.
                # Technically, we need to increment all recall values lower than the given
                # recall value (by the same logic as the comment above),
                # But we don't need to do this in the inner loop here. We can do it at the end
                # where we only need to do it once for each matrix instead of once per bin per matrix.
                for (precision_index, recall) in enumerate(v)
                    recall_index = searchsortedlast(recalls, recall)
                    iszero(recall_index) && continue
                    m[recall_index, precision_index] += 1
                end
            end
        end
    end

    # Just like above with the bin vectors, we need to make the genome/clade vectors reverse max
    for (v1, v2) in values(max_genome_recall_at_precision)
        make_reverse_max!(v1)
        make_reverse_max!(v2)
    end
    for (v1, v2) in values(max_clade_recall_at_precision)
        make_reverse_max!(v1)
        make_reverse_max!(v2)
    end

    # Now make the matrices counting recovered genomes / clades. Similar to above for
    # the bins, we increment the matrix at the correct recall value.
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

    # For all the matrices below, if a bin/genome/whatever is seen at recalls X and
    # precision Y, then it's also seen at every lower value. The `make_reverse_max!` calls
    # earlier updated the lower precision values with the recall values of the higher precision values.
    # Now, we need to do the reverse - update the low recall values with prec. of high rec. values.
    # We can do this once here, in the matrix.
    for mv in (asm_matrices, genome_matrices, bin_genome_matrices, bin_asm_matrices),
        m in mv

        make_columnwise_reverse_cumulative!(m)
    end

    # For historical reasons, the matrices above are transposed.
    # Instead of rewriting this massive function, simply transpose each matrix before returning
    map(
        v -> map(permutedims, v),
        (asm_matrices, genome_matrices, bin_asm_matrices, bin_genome_matrices),
    )
end

"For each precision column in the matrix, add one to the correct row
given by the recall value at the given precision"
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
        iszero(recall_index) && continue
        imax = min(recall_index, imax)
        matrix[recall_index, precision_index] += 1
    end
    matrix
end

function make_reverse_max!(v::Vector{<:Real})
    @inbounds for i in (length(v) - 1):-1:1
        v[i] = max(v[i], v[i + 1])
    end
    v
end

function make_columnwise_reverse_cumulative!(m::Matrix{<:Real})
    for col in eachcol(m)
        for i in (length(col) - 1):-1:1
            col[i] += col[i + 1]
        end
    end
end
