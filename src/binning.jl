const DEFAULT_RECALLS = (0.6, 0.7, 0.8, 0.9, 0.95, 0.99)
const DEFAULT_PRECISIONS = (0.6, 0.7, 0.8, 0.9, 0.95, 0.99)

struct BinStats
    # Recall and precision are computed by pairing each bin
    # with the Genome that yields the highest F1 score
    mean_bin_recall::Float64
    mean_bin_precision::Float64
    mean_bin_f1::Float64
end

function BinStats(bins::Vector{Bin}; assembly::Bool = true)
    (mean_bin_recall, mean_bin_precision, mean_bin_f1) =
        mean_bin_recall_prec(bins; assembly)
    return BinStats(mean_bin_recall, mean_bin_precision, mean_bin_f1)
end

function mean_bin_recall_prec(bins::Vector{Bin}; assembly::Bool = true)
    (recall_sum, prec_sum, f1_sum, nbins) = foldl(
        Iterators.filter(
            !isnothing,
            Iterators.map(b -> recall_prec_max_f1(b; assembly), bins),
        );
        init = (0.0, 0.0, 0.0, 0),
    ) do (recall_sum, prec_sum, f1_sum, nbins), (; recall, precision)
        (recall_sum + recall, prec_sum + precision, f1_sum + f1(recall, precision), nbins + 1)
    end
    return (recall_sum / nbins, prec_sum / nbins, f1_sum / nbins)
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
        filter_genomes=Returns(true),
        filter_bins=Returns(true),
)
```

* `min_size`: Filter away bins with breadth lower than this
* `min_seqs`: Filter away bins with fewer sequences that this
* `binsplit_separator`: Split bins based on this separator
  (`nothing` means no binsplitting)
* `disjoint`: Throw an error if the same sequence is seen in multiple bins.
  Regardless of this option, the same sequence may not be present in the
  same bin more than once.
* `recalls` and `precision`: The thresholds to benchmark with
* `filter_genomes`: A function `f(genome)::Bool`. Genomes for which it returns
   `false` are ignored in benchmarking.
* `filter_bins`: A function `f(bin)::Bool`. Bins for which it returns `false` are
   removed from the `Binning`. This filter is applies after `min_size` and `min_seqs`
   filter.
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
    is_disjoint::Bool
end

function Base.show(io::IO, x::Binning)
    nc = n_nc(x)
    print(io, summary(x), '(')
    if nc isa Integer
        print(io, "NC = ", nc)
    end
    return print(io, ')')
end

function Base.show(io::IO, ::MIME"text/plain", x::Binning)
    return if get(io, :compact, false)
        show(io, x)
    else
        buf = IOBuffer()
        show(buf, MIME"text/plain"(), x.ref)
        seekstart(buf)
        println(io, "Binning")
        for line in eachline(buf)
            println(io, "  ", line)
        end
        print(io, "  Bins:        ", n_bins(x))
        nc = n_nc(x)
        if nc isa Integer
            print(io, "\n  NC genomes:  ", nc)
        end
        npb = extract_if_present(x, x.bin_genomes[1], 0.9, 0.95)
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
                round(stats.mean_bin_recall; digits = 3),
                " / ",
                round(stats.mean_bin_precision; digits = 3),
                " / ",
                round(stats.mean_bin_f1; digits = 3),
            )
        end
        print(io, "\n  Precisions: ", repr([round(i; digits = 3) for i in x.precisions]))
        print(io, "\n  Recalls:    ", repr([round(i; digits = 3) for i in x.recalls]))
        print(io, "\n  Reconstruction (genomes):\n")
        # Use a new buffer, since writing to the old one would not truncate it
        buf = IOBuffer()
        print_matrix(buf, x; level = 0, assembly = false)
        seekstart(buf)
        for line in eachline(buf)
            println(io, "    ", line)
        end
    end
end

const DEFAULT_THRESHOLD_ATOL = 1.0e-6

"""
    n_recovered(::Binning, recall, precision; level=0, assembly=false, atol=1e-6)::Integer

Return the number of genomes or clades reconstructed in the `Binning` at the given recall
and precision thresholds.
If `assembly` is set, return the number of assemblies reconstructed instead.
The argument `level` sets the taxonomic rank: 0 for `Genome` (or assemblies).

`recall` and `precision` must each be within `atol` of one of the thresholds the
`Binning` was benchmarked with (`binning.recalls` and `binning.precisions`),
else an `ArgumentError` is thrown.
To get bounds for thresholds between these, use [`n_recovered_bounds`](@ref).

See also: [`n_passing_bins`](@ref)

# Examples
```jldoctest
julia> n_recovered(binning, 0.6, 0.8)
1

julia> n_recovered(binning, 0.6, 0.8; assembly=true)
2

julia> n_recovered(binning, 0.6, 0.8; assembly=true, level=2)
1
```
"""
function n_recovered(
        binning::Binning,
        recall::Real,
        precision::Real;
        level::Integer = 0,
        assembly::Bool = false,
        atol::Real = DEFAULT_THRESHOLD_ATOL,
    )::Int
    m = get_matrix(assembly ? binning.recovered_asms : binning.recovered_genomes, level)
    return extract_exact(binning, m, recall, precision, atol)
end

"""
    n_recovered_bounds(
        ::Binning, recall, precision; level=0, assembly=false, atol=1e-6
    )::@NamedTuple{lower::Int, upper::Int}

Return lower and upper bounds for the number of genomes or clades reconstructed in the
`Binning` at the given recall and precision thresholds.

The number is only known at the thresholds in `binning.recalls` and `binning.precisions`.
For other values, the lower bound is the number at the nearest benchmarked thresholds
above the given ones, and the upper bound is the number at the nearest thresholds below.
A threshold within `atol` of a benchmarked threshold is taken to be that threshold.
If both thresholds are benchmarked, `lower == upper`.

Throws an `ArgumentError` if `recall` or `precision` is not in [0, 1], or is outside the
range of benchmarked thresholds.
Keyword arguments `level` and `assembly` are as for [`n_recovered`](@ref).

# Examples
```jldoctest
julia> n_recovered_bounds(binning, 0.6, 0.8)
(lower = 1, upper = 1)

julia> n_recovered_bounds(binning, 0.65, 0.75; assembly=true)
(lower = 1, upper = 2)
```
"""
function n_recovered_bounds(
        binning::Binning,
        recall::Real,
        precision::Real;
        level::Integer = 0,
        assembly::Bool = false,
        atol::Real = DEFAULT_THRESHOLD_ATOL,
    )::@NamedTuple{lower::Int, upper::Int}
    m = get_matrix(assembly ? binning.recovered_asms : binning.recovered_genomes, level)
    return extract_bounds(binning, m, recall, precision, atol)
end

"""
    n_passing_bins(::Binning, recall, precision; level=0, assembly=false, atol=1e-6)::Integer

Return the number of bins which correspond to any genome or clade at the given recall
and precision thresholds.
If `assembly` is set, a recall of 1.0 means a bin corresponds to a whole assembly,
else it corresponds to a whole genome.
The argument `level` sets the taxonomic rank: 0 for `Genome` (or assemblies).

`recall` and `precision` must each be within `atol` of one of the thresholds the
`Binning` was benchmarked with (`binning.recalls` and `binning.precisions`),
else an `ArgumentError` is thrown.
To get bounds for thresholds between these, use [`n_passing_bins_bounds`](@ref).

See also: [`n_recovered`](@ref)

# Examples
```jldoctest
julia> n_passing_bins(binning, 0.6, 0.8)
1

julia> n_passing_bins(binning, 0.7, 0.8)
0
```
"""
function n_passing_bins(
        binning::Binning,
        recall::Real,
        precision::Real;
        level::Integer = 0,
        assembly::Bool = false,
        atol::Real = DEFAULT_THRESHOLD_ATOL,
    )::Int
    m = get_matrix(assembly ? binning.bin_asms : binning.bin_genomes, level)
    return extract_exact(binning, m, recall, precision, atol)
end

"""
    n_passing_bins_bounds(
        ::Binning, recall, precision; level=0, assembly=false, atol=1e-6
    )::@NamedTuple{lower::Int, upper::Int}

Return lower and upper bounds for the number of bins which correspond to any genome or
clade at the given recall and precision thresholds, which need not be ones the
`Binning` was benchmarked with.
The bounds are computed as for [`n_recovered_bounds`](@ref), and the keyword
arguments are as for [`n_passing_bins`](@ref).

# Examples
```jldoctest
julia> n_passing_bins_bounds(binning, 0.6, 0.8)
(lower = 1, upper = 1)

julia> n_passing_bins_bounds(binning, 0.65, 0.8)
(lower = 0, upper = 1)
```
"""
function n_passing_bins_bounds(
        binning::Binning,
        recall::Real,
        precision::Real;
        level::Integer = 0,
        assembly::Bool = false,
        atol::Real = DEFAULT_THRESHOLD_ATOL,
    )::@NamedTuple{lower::Int, upper::Int}
    m = get_matrix(assembly ? binning.bin_asms : binning.bin_genomes, level)
    return extract_bounds(binning, m, recall, precision, atol)
end

function get_matrix(matrices::Vector{<:Matrix}, level::Integer)::Matrix{Int}
    if level + 1 ∉ eachindex(matrices)
        throw(
            ArgumentError(
                "Level $level is out of range: must be in 0:$(length(matrices) - 1)"
            )
        )
    end
    return matrices[level + 1]
end

function extract_exact(
        binning::Binning,
        m::Matrix{Int},
        recall::Real,
        precision::Real,
        atol::Real,
    )::Int
    ri = exact_threshold_index(binning.recalls, recall, "Recall", atol)
    pi = exact_threshold_index(binning.precisions, precision, "Precision", atol)
    return m[pi, ri]
end

function extract_bounds(
        binning::Binning,
        m::Matrix{Int},
        recall::Real,
        precision::Real,
        atol::Real,
    )::@NamedTuple{lower::Int, upper::Int}
    (rlo, rhi) = bounding_threshold_indices(binning.recalls, recall, "Recall", atol)
    (plo, phi) = bounding_threshold_indices(binning.precisions, precision, "Precision", atol)
    # Counts can only decrease when either threshold increases
    return (; lower = m[phi, rhi], upper = m[plo, rlo])
end

function check_unit_interval(x::Real, name::String)
    return if !(0 ≤ x ≤ 1)
        throw(ArgumentError("$name threshold $x is not in [0, 1]"))
    end
end

# Index of the threshold within `atol` of `x`, or `nothing` if there is none.
function find_threshold_index(
        thresholds::Vector{Float64},
        x::Real,
        atol::Real,
    )::Union{Int, Nothing}
    i = searchsortedfirst(thresholds, x - atol)
    return if i ≤ lastindex(thresholds) && thresholds[i] ≤ x + atol
        i
    else
        nothing
    end
end

function exact_threshold_index(
        thresholds::Vector{Float64},
        x::Real,
        name::String,
        atol::Real,
    )::Int
    check_unit_interval(x, name)
    i = find_threshold_index(thresholds, x, atol)
    if isnothing(i)
        throw(
            ArgumentError(
                "$name threshold $x is not one of the benchmarked thresholds $thresholds. " *
                    "Use the `_bounds` function variant for other thresholds.",
            )
        )
    end
    return i
end

# Indices of the nearest thresholds at or below, and at or above `x`.
function bounding_threshold_indices(
        thresholds::Vector{Float64},
        x::Real,
        name::String,
        atol::Real,
    )::Tuple{Int, Int}
    check_unit_interval(x, name)
    i = find_threshold_index(thresholds, x, atol)
    isnothing(i) || return (i, i)
    if x < first(thresholds) || x > last(thresholds)
        throw(
            ArgumentError(
                "$name threshold $x is outside the range of benchmarked thresholds " *
                    "[$(first(thresholds)), $(last(thresholds))]",
            )
        )
    end
    hi = searchsortedfirst(thresholds, x)
    return (hi - 1, hi)
end

"""
    print_matrix(::Binning; level=0, assembly=false)

Print the number of reconstructed assemblies or genomes at the given taxonomic level (rank).
Level 0 corresponds to genomes, level 1 to species, etc.
If `assembly`, print the number of reconstructed assemblies, else print the level
of reconstructed genomes.

See also: [`Binning`](@ref)
"""
print_matrix(x::Binning; kwargs...) = print_matrix(stdout, x; kwargs...)
function print_matrix(io::IO, x::Binning; level::Integer = 0, assembly::Bool = false)
    ms = assembly ? x.recovered_asms : x.recovered_genomes
    m = ms[level + 1]
    rnd(x) = string(round(x; digits = 3))
    digitwidth(x) = sizeof(rnd(x))
    width = max(
        maximum(digitwidth, m; init = 0) + 1,
        maximum(digitwidth, x.recalls; init = 0) + 1,
    )
    col1_width = max(3, maximum(digitwidth, x.precisions))
    println(io, rpad("P\\R", col1_width), join([lpad(rnd(i), width) for i in x.recalls]))
    for (prec_index, prec) in enumerate(x.precisions)
        println(
            io,
            rpad(rnd(prec), col1_width),
            join([lpad(i, width) for i in m[prec_index, :]]),
        )
    end
    return
end

# Count at the given thresholds, or `nothing` if they were not benchmarked.
function extract_if_present(
        x::Binning,
        m::Matrix{Int},
        recall::Real,
        precision::Real,
    )::Union{Int, Nothing}
    ri = find_threshold_index(x.recalls, recall, DEFAULT_THRESHOLD_ATOL)
    pi = find_threshold_index(x.precisions, precision, DEFAULT_THRESHOLD_ATOL)
    (isnothing(ri) || isnothing(pi)) && return nothing
    return m[pi, ri]
end

# Number of near-complete genomes, or `nothing` if the thresholds were not benchmarked
n_nc(x::Binning)::Union{Int, Nothing} = extract_if_present(x, x.recovered_genomes[1], 0.9, 0.95)

n_bins(x::Binning) = length(x.bins)

"""
    is_disjoint(::Binning) -> Bool

Check if the `Binning` is disjoint, i.e. if every sequence in the binning
is present in only one bin.

# Examples
```jldoctest
julia> is_disjoint(binning)
true

julia> is_disjoint(gold_standard(ref; disjoint=false))
false
```
"""
is_disjoint(x::Binning) = x.is_disjoint

function Binning(path::AbstractString, ref::Reference; kwargs...)
    return open_perhaps_gzipped(io -> Binning(io, ref; kwargs...), String(path))
end

@noinline function throw_disjoint_error()
    throw(ArgumentError("Binning is not disjoint, but `disjoint` set to `true`"))
end

# This is the most common constructor in practice, because we usually load binnings from files,
# and also the most efficient. I immediately convert the sequences to integers for faster
# processing, then convert back to seqs before instantiating the bins.
function Binning(
        io::IO,
        ref::Reference;
        min_size::Integer = 1,
        min_seqs::Integer = 1,
        binsplit_separator::Union{AbstractString, Char, Nothing} = nothing,
        disjoint::Bool = true,
        recalls = DEFAULT_RECALLS,
        precisions = DEFAULT_PRECISIONS,
        filter_genomes::Function = Returns(true),
        filter_bins::Function = Returns(true),
    )
    (idxs_by_binname, first_non_disjoint) = parse_bins(io, Dict, ref, binsplit_separator)
    is_disjoint = first_non_disjoint === nothing
    if disjoint && !is_disjoint
        (line_number, seqname) = first_non_disjoint
        throw(
            BinningError(
                "Sequence \"$(seqname)\" is present in multiple bins, but `disjoint` is set to `true`",
                line_number,
            ),
        )
    end
    filter!(idxs_by_binname) do (_, idxs)
        length(idxs) ≥ min_seqs &&
            sum((length(first(ref.targets[i])) for i in idxs); init = 0) ≥ min_size
    end
    considered_genomes = _filter_genomes(filter_genomes, ref)
    bins = _getbins(ref, considered_genomes, idxs_by_binname)
    filter_bins === Returns(true) || filter!(filter_bins, bins)
    return _binning(bins, ref, recalls, precisions, is_disjoint)
end

function _filter_genomes(f::Function, ref::Reference)::Union{Nothing, Set{Genome}}
    return if f === Returns(true)
        nothing
    else
        Set(g for g in genomes(ref) if f(g))
    end
end

function _getbins(
        ref::Reference,
        considered_genomes::Union{Nothing, Set{Genome}},
        idxs_by_binname::Dict{SubString{String}, Vector{UInt32}},
    )::Vector{Bin}
    scratch = Tuple{Int, Int}[]
    return [
        bin_by_indices(binname, seq_idxs, ref.targets, scratch, considered_genomes)::Bin for
            (binname, seq_idxs) in idxs_by_binname
    ]
end


"""
    Binning(bins::Vector{Bin}, ref::Reference; kwargs...)

Create a `Binning` from a vector of `Bin`.
The allowed keyword arguments are: `recalls`, `precisions`, `disjoint`,
see the main docstring for `Binning` for their meaning.

```jldoctest
julia> binning = Binning([bin], ref);

julia> n_bins(binning)
1
```
"""
function Binning(
        bins::Vector{Bin},
        ref::Reference;
        recalls = DEFAULT_RECALLS,
        precisions = DEFAULT_PRECISIONS,
        disjoint::Bool = true,
    )
    is_disjoint = check_disjoint(bins)
    disjoint && !is_disjoint && throw_disjoint_error()
    return _binning(bins, ref, recalls, precisions, is_disjoint)
end

function _binning(
        bins::Vector{Bin},
        ref::Reference,
        recalls,
        precisions,
        is_disjoint::Bool,
    )
    bins = sort(bins; by = i -> i.name)
    checked_recalls = validate_recall_precision(recalls)
    checked_precisions = validate_recall_precision(precisions)
    bin_asm_stats = BinStats(bins; assembly = true)
    bin_genome_stats = BinStats(bins; assembly = false)
    (asm_matrices, genome_matrices, bin_asms, bin_genomes) =
        benchmark(ref, bins, checked_recalls, checked_precisions)
    return Binning(
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
        is_disjoint,
    )
end

"""
    gold_standard(
        ref::Reference
        [sequences, a Binning or an iterable of Sequence];
        binsplit_separator::Union{AbstractString, Char, Nothing}=nothing,
        disjoint=true,
        recalls=DEFAULT_RECALLS,
        precisions=DEFAULT_PRECISIONS
    )::Binning

Create the optimal `Binning` object given a `Reference`, by the optimal binning of
the `Sequence`s in `sequences`.
If `disjoint`, assign each sequence to only a single genome.

If `sequences` is not passed, use all sequences in `ref`. If a `Binning` is passed,
use all sequences in any of its bins. Else, pass an iterable of `Sequence`.

# Examples
```jldoctest
julia> gs = gold_standard(ref);

julia> gs isa Binning
true

julia> n_recovered(gs, 0.6, 0.8) >= n_recovered(binning, 0.6, 0.8)
true
```

# Extended help
Currently, the `disjoint` option uses a simple greedy algorithm to assign
sequences to genomes.
"""
function gold_standard(ref::Reference; kwargs...)::Binning
    return gold_standard(ref, [first(v) for v in ref.targets]; kwargs...)
end

function gold_standard(ref::Reference, binning::Binning; kwargs...)::Binning
    seqs = reduce(binning.bins; init = Sequence[]) do v, bin
        append!(v, bin.sequences)
    end
    is_disjoint(binning) || unique!(seqs)
    return gold_standard(ref, seqs; kwargs...)
end

function gold_standard(
        ref::Reference,
        sequences; # iterable of Sequence
        binsplit_separator::Union{AbstractString, Char, Nothing} = nothing,
        disjoint::Bool = true,
        recalls = DEFAULT_RECALLS,
        precisions = DEFAULT_PRECISIONS,
    )::Binning

    # First, assign each sequence to the genome(s) it maps to - only one if disjoint,
    # else all the genomes it maps to.
    sequences_of_genome = Dict{Genome, Vector{Sequence}}()
    for sequence in sequences
        targets = last(ref.targets[sequence_index(ref, sequence)])
        isempty(targets) && continue
        # If disjoint, pick the genome where the sequence aligns to the longest span
        best_index = disjoint ? last(findmax(span_length, targets)) : 0
        for (i, (source, _)) in enumerate(targets)
            disjoint && i != best_index && continue
            push!(
                get!(valtype(sequences_of_genome), sequences_of_genome, source.genome),
                sequence,
            )
        end
    end
    # Now, build the bins. Use the binsplit iterator if given.
    bins = Bin[]
    scratch = Tuple{Int, Int}[]
    binnames_indices = Tuple{String, Vector{UInt32}}[] # indices in ref.targets
    # Used when binsplitting. Values are the indices of the seqs.
    seqs_by_sample = Dict{SubString{String}, Vector{UInt32}}()
    for (genome, seqs) in sequences_of_genome
        # If not disjoint, then we could have the same sequence mapping to two locations
        # in the same genome. However, we don't want the same bin to contain the same sequence
        # multiple times, even when !disjoint, so we remove duplicates here.
        disjoint || unique!(seqs)
        empty!(binnames_indices)
        if isnothing(binsplit_separator)
            push!(binnames_indices, ("bin_" * genome.name, [ref.target_index_by_name[i.name] for i in seqs]))
        else
            empty!(seqs_by_sample)
            for sequence in seqs
                splt = @inline split_once(sequence.name, binsplit_separator)
                (prefix, _) = @something splt throw(
                    ArgumentError(
                        "Binsplit separator $(repr(binsplit_separator)) not found in sequence name $(repr(sequence.name))"
                    )
                )
                index = ref.target_index_by_name[sequence.name]
                push!(get!(valtype(seqs_by_sample), seqs_by_sample, prefix), index)
            end
            for (prefix, indices) in seqs_by_sample
                push!(binnames_indices, ("bin_" * prefix * genome.name, indices))
            end
        end
        for (name, indices) in binnames_indices
            bin = bin_by_indices(name, indices, ref.targets, scratch, nothing)
            push!(bins, bin)
        end
    end
    # Even if we don't force it to be disjoint, it still could be
    bins_are_disjoint = disjoint ? true : check_disjoint(bins)
    return _binning(bins, ref, recalls, precisions, bins_are_disjoint)
end

# Number of positions in the source covered by the span, accounting for
# spans that wrap around circular sources (from > to)
function span_length((source, (from, to))::Target)::Int
    return from ≤ to ? to - from + 1 : source.length - from + to + 1
end

function check_disjoint(bins)
    nseq = sum(i -> length(i.sequences), bins; init = 0)
    seen_seqs = sizehint!(Set{Sequence}(), nseq)
    for bin in bins, seq in bin.sequences
        in!(seen_seqs, seq) && return false
    end
    return true
end

function validate_recall_precision(xs)::Vector{Float64}
    s = Set{Float64}()
    for x_ in xs
        x = Float64(x_)
        x in s && throw(ArgumentError("Recall/precision value $x present multiple times"))
        if !isfinite(x) || x <= 0.0 || x > 1.0
            throw(ArgumentError("Recall/precision value $x is not in (0,1]"))
        end
        push!(s, x)
    end
    isempty(s) && throw(ArgumentError("Must provide at least 1 recall/precision value"))
    return sort!(collect(s))
end

# TODO: This function is massive. Can I refactor it into smaller functions?
function benchmark(
        ref::Reference,
        bins::Vector{Bin},
        recalls::Vector{Float64},
        precisions::Vector{Float64}
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
        Dict{Clade, Tuple{Vector{Float64}, Vector{Float64}}}()
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
            asm_recall = asmsize / assembly_size(genome)
            genome_recall = asmsize / genome_size(genome)

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
    asm_matrices::Vector{Matrix{Int}} = [zeros(Int, length(recalls), length(precisions)) for i in 1:nranks(ref)]
    genome_matrices::Vector{Matrix{Int}} = [copy(i) for i in asm_matrices]
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
    for mv in Vector{Matrix{Int}}[asm_matrices, genome_matrices, bin_genome_matrices, bin_asm_matrices], m in mv
        make_columnwise_reverse_cumulative!(m)
    end

    # For historical reasons, the matrices above are transposed.
    # Instead of rewriting this massive function, simply transpose each matrix before returning
    return map(
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
    return matrix
end

function make_reverse_max!(v::Vector{<:Real})
    @inbounds for i in (length(v) - 1):-1:1
        v[i] = max(v[i], v[i + 1])
    end
    return v
end

function make_columnwise_reverse_cumulative!(m::Matrix{<:Real})
    for col in eachcol(m)
        for i in (length(col) - 1):-1:1
            col[i] += col[i + 1]
        end
    end
    return
end

@noinline function throw_different_binnings_error()
    throw(ArgumentError("Binnings do not have the same sequences"))
end

"""
    adjusted_rand_index(Binning, Binning; intersection::Bool=false) -> Float64

Compute the adjusted rand index (ARI) between two binnings.
The result is a real number between -0.5 and 1.0, inclusive.
The binnings must use the same reference, and must be disjoint.

If `intersection` is `false` (the default), an `ArgumentError` is thrown if the two binnings
do not contain exactly the same sequences. If `true`, ARI will be computed over
only the sequences the two binnings have in common, and the others will be ignored.

To compute the ARI against the ground truth of a reference `ref`, compare to the result
of `gold_standard(ref)`.
Note that if the input binning is binsplit, consider also binsplitting the gold standard.

!!! warning
    ARI is generally not an appropriate metric for evaluating binnings, and is only
    included in this package for comprehensiveness.
    See [the section on ARI in the documentation](@ref ari) for more details.

# Examples
```jldoctest
julia> adjusted_rand_index(binning, binning)
1.0

julia> adjusted_rand_index(binning, gold_standard(ref); intersection=true)
0.4283752695493192

julia> adjusted_rand_index(binning, gold_standard(ref))
ERROR: ArgumentError: Binnings do not have the same sequences
[...]

```
See also: [`gold_standard`](@ref), [`is_disjoint`](@ref)
"""
function adjusted_rand_index(a::Binning, b::Binning; intersection::Bool = false)
    # See the end of the function for the definition of ARI, to see what e.g. a_sums mean.

    if a.ref !== b.ref
        throw(ArgumentError("To compute ARI, the two binnings must have the same reference"))
    end

    if !is_disjoint(a) || !is_disjoint(b)
        throw(ArgumentError("Binnings must be disjoint to compute ARI, but one or more is not."))
    end

    # Two empty binnings are identical
    if isempty(a.bins) && isempty(b.bins)
        return 1.0
    end

    # Test the binnings have the same sequences.
    # ARI is invalid if the set of sequences in the two binnings are not identical.
    # First, we check they both are disjoint. Then, that they have the same number of seqs.
    # Finally, further below, for each seq in `a`, we check it's present in `b`.
    # That is sufficient to check the set of sequences are the same.
    if !intersection
        n_sequences(x::Binning) = sum(i -> length(i.sequences), x.bins; init = 0)
        n_sequences(a) == n_sequences(b) || throw_different_binnings_error()
    end

    ref = b.ref
    target_index_by_name = ref.target_index_by_name

    # This vector is not n_choose_two since it is filled incrementally
    b_sums = zeros(UInt, n_bins(b))
    a_sum_choose_two = 0
    sum_n_choose_two = 0
    n = 0

    # We store both bins and sequences as integers, to make the core algorithm faster.
    # A missing bin is typemax(UInt32)
    b_bin_index_of_seq_index = fill(typemax(UInt32), n_seqs(ref))
    for (bin_index, bin) in enumerate(b.bins), sequence in bin.sequences
        seq_index = target_index_by_name[sequence.name]
        b_bin_index_of_seq_index[seq_index] = bin_index % UInt32
    end

    # B Bin index => number of overlapping basepairs between bin and the current A bin.
    bin_index_lens = Dict{Int, UInt}()
    for bin_a in a.bins
        # This occurs shockingly often, so is worth special casing
        if length(bin_a.sequences) == 1
            sequence = only(bin_a.sequences)
            seq_index = target_index_by_name[sequence.name]
            bin_b_index = b_bin_index_of_seq_index[seq_index]
            if bin_b_index == typemax(UInt32)
                intersection ? continue : throw_different_binnings_error()
            end
            n += length(sequence)
            sum_n_choose_two += n_choose_two(length(sequence) % UInt)
            a_sum_choose_two += n_choose_two(length(sequence) % UInt)
            b_sums[bin_b_index] += length(sequence) % UInt
        else
            # In the general case, we count up the total basepairs that overlap
            # for each bin.
            # That is, for each seq in bin A, we find the corresponding bin B,
            # where the seq is present, then add length(seq) to binA, binB overlap
            empty!(bin_index_lens)
            for sequence in bin_a.sequences
                seq_index = target_index_by_name[sequence.name]
                bin_b_index = b_bin_index_of_seq_index[seq_index]
                if bin_b_index == typemax(UInt32)
                    intersection ? continue : throw_different_binnings_error()
                end
                new_len = get!(bin_index_lens, bin_b_index, UInt(0)) + length(sequence) % UInt
                bin_index_lens[bin_b_index] = new_len
            end
            a_sum = 0
            for (bin_b_index, intersection) in bin_index_lens
                a_sum += intersection
                b_sums[bin_b_index] += intersection
                n += intersection
                sum_n_choose_two += n_choose_two(intersection)
            end
            a_sum_choose_two += n_choose_two(a_sum % UInt)
        end
    end

    b_sum_choose_two = sum(i -> n_choose_two(i), b_sums; init = UInt(0))
    n_2 = n_choose_two(n % UInt)

    # If there is no overlap, either both binnings are empty (we checked that)
    # or else there is no overlap and they are 100% distinct
    iszero(n_2) && return 0.0

    # The denominator below is zero iff both binnings are the same trivial partition,
    # i.e. both are a single bin, or both consist only of single-basepair bins.
    # ARI is undefined (0/0) here, but by convention, identical partitions have ARI 1.
    # We check with exact integer arithmetic to avoid rounding issues.
    if a_sum_choose_two == b_sum_choose_two && (iszero(a_sum_choose_two) || a_sum_choose_two == n_2)
        return 1.0
    end

    # Convert to Float64 before multiplying, since the product of the
    # integers may overflow for large binnings.
    (a_f, b_f, n_2_f) = (Float64(a_sum_choose_two), Float64(b_sum_choose_two), Float64(n_2))
    expected = (a_f * b_f) / n_2_f
    numerator = Float64(sum_n_choose_two) - expected
    denominator = (a_f + b_f) / 2 - expected

    # Handle rounding errors
    return clamp(numerator / denominator, -0.5, 1.0)
end
