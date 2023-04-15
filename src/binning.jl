const DEFAULT_RECALLS = (0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99)
const DEFAULT_PRECISIONS = (0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99)

"""
    Binning

A `Binning` represents a set of `Bin`s benchmarked against a `Reference`.
`Binning`s can be created given a set of `Bin`s and a `Reference`, where the
bins may potentially be loaded from a `.tsv` file.
The field `binning.recoverable_genomes` shows the maximal number of recoverable genomes
at given recall levels given perfect binning.
The fields `recovered_asms` and `recovered_genomes` are used for benchmarking,
these are normally output using the `print_matrix` function.

See also: [`print_matrix`](@ref), [`Bin`](@ref), [`Reference`](@ref)

# Examples
```jldoctest
julia> bins = gold_standard(ref);

julia> bins isa Binning
true

julia> VambBenchmarks.n_nc(binning)
0
```

# Extended help
Create with:
```julia
open(file) do io
    Binning(
        io::IO,
        ref::Reference;
        min_size::Integer=1,
        min_seqs::Integer=1,
        binsplit_separator::Union{AbstractString, Char, Nothing}=nothing,
        disjoint::Bool=true,
        recalls=DEFAULT_RECALLS,
        precisions=DEFAULT_PRECISIONS
)
```

* `min_size`: Filter away bins with breadth lower than this
* `min_seqs`: Filter away bins with fewer sequences that this
* `binsplit_separator`: Split bins based on this separator
  (`nothing` means no binsplitting)
* `disjoint`: Throw an error if the same sequence is seen in multiple bins
* `recalls` and `precision`: The thresholds to benchmark with
"""
struct Binning
    ref::Reference
    bins::Vector{Bin}
    # One matrix per rank - first is genome, then upwards
    recovered_asms::Vector{Matrix{Int}}
    recovered_genomes::Vector{Matrix{Int}}
    # If perfect binning, this number of genomes could be recovered at the
    # various recall levels
    recoverable_genomes::Vector{Int}
    recalls::Vector{Float64}
    precisions::Vector{Float64}
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
        print(io,
            "  Bins:        ", nbins(x)
        )
        nc = n_nc(x)
        if nc !== nothing
            print(io, "\n  NC genomes:  ", nc)
        end
        print(io, "\n  Precisions: ", repr([round(i; digits=3) for i in x.precisions]))
        print(io, "\n  Recalls:    ", repr([round(i; digits=3) for i in x.recalls]))
        print(io, "\n  Recoverable genomes: ", repr(x.recoverable_genomes))
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
    rnd(x) = string(round(x, digits=3))
    digitwidth(x) = sizeof(rnd(x))
    width = max(
        max(4, ndigits(maximum(m) + 1)),
        maximum(digitwidth, x.recalls) + 1
    )
    col1_width = max(3, maximum(digitwidth, x.precisions))
    println(io, rpad("P\\R", col1_width), join([lpad(i, width) for i in x.recalls]))
    for (prec_index, prec) in enumerate(x.precisions)
        println(io, rpad(rnd(prec), col1_width), join([lpad(i, width) for i in m[prec_index, :]]))
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

function Binning(
    io::IO,
    ref::Reference;
    min_size::Integer=1,
    min_seqs::Integer=1,
    binsplit_separator::Union{AbstractString, Char, Nothing}=nothing,
    disjoint::Bool=true,
    recalls=DEFAULT_RECALLS,
    precisions=DEFAULT_PRECISIONS,
)
    bins = sort!(parse_bins(io, ref, binsplit_separator), by=i -> i.name)
    filter!(bins) do bin
        nseqs(bin) >= min_seqs && bin.breadth >= min_size
    end
    Binning(bins, ref; recalls=recalls, precisions=precisions, disjoint)
end

function Binning(bins_,
    ref::Reference;
    recalls=DEFAULT_RECALLS,
    precisions=DEFAULT_PRECISIONS,
    disjoint::Bool = true
)
    checked_recalls = validate_recall_precision(recalls)
    checked_precisions = validate_recall_precision(precisions)
    bins = vector(bins_)
    disjoint && check_disjoint(bins)
    recoverable_genomes = zeros(Float64, length(checked_recalls))
    for genome in ref.genomes
        i = searchsortedlast(checked_recalls, genome.assembly_size / genome.genome_size)
        iszero(i) && continue
        recoverable_genomes[i] += 1
    end
    # If a genome is recoverable at recall 0.5, it's recoverable at recalls < 0.5
    for i in length(recoverable_genomes)-1:-1:1
        recoverable_genomes[i] += recoverable_genomes[i + 1]
    end
    (asm_matrices, genome_matrices) = benchmark(ref, bins, checked_recalls, checked_precisions)
    Binning(ref, bins, asm_matrices, genome_matrices, recoverable_genomes, checked_recalls, checked_precisions)
end

"""
    gold_standard(
        ref::Reference;
        recalls=DEFAULT_RECALLS,
        precisions=DEFAULT_PRECISIONS
    )::Binning

Create the optimal `Binning` object given an assembly. If sequences map to multiple genomes,
the binning is not guaranteed to be disjoint.
"""
function gold_standard(
    ref::Reference;
    recalls=DEFAULT_RECALLS,
    precisions=DEFAULT_PRECISIONS
)::Binning
    sequences_of_genome = Dict{Genome, Set{Sequence}}()
    for (_, (sequence, targets)) in ref.targets_by_name, (source, _) in targets
        push!(get!(valtype(sequences_of_genome), sequences_of_genome, source.genome), sequence)
    end
    bins = [Bin("bin_" * genome.name, collect(seqs), ref.targets_by_name) for (genome, seqs) in sequences_of_genome]
    sort!(bins, by=i -> i.name)
    Binning(bins, ref; recalls=recalls, precisions=precisions, disjoint=false)
end

function check_disjoint(bins)
    nseq = sum(i -> length(i.sequences), bins)
    seen_seqs = sizehint!(Set{Sequence}(), nseq)
    for bin in bins, seq in bin.sequences
        if seq in seen_seqs
            error(lazy"Sequence \"$(seq.name)\" seen twice in disjoint Binning")
        end
        push!(seen_seqs, seq)
    end
    nothing
end

function validate_recall_precision(xs)
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
    precisions::Vector{Float64}
)
    # For each genome/clade, we compute the maximal recall at the given precision levels.
    # i.e. if 3rd element of vector is 0.5, it means that at precision precisions[3], this genome/clade
    # is found with a maximal recall of 0.5.
    # We keep two vectors: First for recovered assemblies, second for recovered genomes
    max_genome_recall_at_precision = Dict{Genome, Tuple{Vector{Float64}, Vector{Float64}}}()
    max_clade_recall_at_precision = Dict{Clade{Genome}, Tuple{Vector{Float64}, Vector{Float64}}}()
    # Initialize with zeros for all known genomes/clades
    for genome in ref.genomes
        max_genome_recall_at_precision[genome] = (zeros(Float64, length(precisions)), zeros(Float64, length(precisions)))
    end
    for level in ref.clades, clade in level
        max_clade_recall_at_precision[clade] = (zeros(Float64, length(precisions)), zeros(Float64, length(precisions)))
    end
    for bin in bins
        for (genome, (asmsize, foreign)) in bin.genomes
            (v_asm, v_genome) = max_genome_recall_at_precision[genome]
            precision = (asmsize) / (asmsize + foreign)
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

function update_matrix!(matrix::Matrix{<:Integer}, v::Vector{<:AbstractFloat}, recalls::Vector{Float64})
    for (precision_index, recall) in enumerate(v)
        # TODO: By keeping track of the maximal recall index (which must shrink)
        # we can reduce search space here, if performance critical
        recall_index = searchsortedlast(recalls, recall)
        matrix[1:recall_index, precision_index] .+= 1
    end
    matrix
end
            
function make_reverse_cumulative!(v::Vector{<:Real})
    for i in length(v)-1:-1:1
        v[i] = max(v[i], v[i+1])
    end
    v
end
