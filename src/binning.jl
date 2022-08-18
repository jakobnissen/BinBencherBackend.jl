const DEFAULT_RECALLS = (0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99)
const DEFAULT_PRECISIONS = (0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99)

struct Binning
    ref::Reference
    bins::Vector{Bin}
    counters::Vector{Matrix{Int}}
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
        seqs = sum(nseqs, x.bins, init=0)
        buf = IOBuffer()
        show(buf, MIME"text/plain"(), x.ref)
        seekstart(buf)
        println(io, "Binning")
        for line in eachline(buf)
            println(io, "  ", line)
        end
        print(io,
            "  Bins:        ", nbins(x),
            "\n  Sequences:   ", seqs
        )
        nc = n_nc(x)
        if nc !== nothing
            print(io, "\n  NC genomes:  ", nc)
        end
        print(io, "\n  Genomes:\n")
        seekstart(buf)
        print_matrix(buf, x, 1)
        seekstart(buf)
        for line in eachline(buf)
            println(io, "    ", line)
        end
    end
end

print_matrix(x::Binning, level::Integer=1) = print_matrix(stdout, x, level)
function print_matrix(io::IO, x::Binning, level::Integer=1)
    rnd(x) = string(round(x, digits=3))
    digitwidth(x) = sizeof(rnd(x))
    m = x.counters[level]
    width = max(4, ndigits(maximum(m) + 1))
    width = max(width, maximum(digitwidth, x.recalls) + 1)
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
    x.counters[1][prec_i, rec_i]
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
    bins = parse_bins(io, ref, binsplit_separator)
    filter!(bins) do bin
        nseqs(bin) >= min_seqs && bin.breadth >= min_size
    end
    Binning(bins, ref; recalls_=recalls, precisions_=precisions, disjoint)
end

function Binning(bins_,
    ref::Reference;
    recalls_=DEFAULT_RECALLS,
    precisions_=DEFAULT_PRECISIONS,
    disjoint::Bool = true
)
    recalls = validate_recall_precision(recalls_)
    precisions = validate_recall_precision(precisions_)
    bins = vector(bins_)
    disjoint && check_disjoint(bins)
    counters = benchmark(ref, bins, recalls, precisions)
    Binning(ref, bins, counters, recalls, precisions)
end

function check_disjoint(bins)
    seen_seqs = Set{Sequence}()
    for bin in bins, seq in bin.sequences
        if seq in seen_seqs
            error("Sequence \"$(seq.name)\" seen twice in disjoint Binning")
        end
        push!(seen_seqs, seq)
    end
    nothing
end

function validate_recall_precision(xs)
    s = Set{Float64}()
    for x_ in xs
        x = Float64(x_)
        x in s && error("Recall/precision value $x present multiple times")
        if !isfinite(x) || x <= 0.0 || x > 1.0
            error("Recall precision value $x is not finite in (0,1]")
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
    counters = Matrix{Int}[]
    rp_by_node = Dict{Node, Dict{Bin, Tuple{Float64, Float64}}}(
        g => Dict{Bin, Tuple{Float64, Float64}}() for g in ref.genomes
    )
    for bin in bins, genome in keys(bin.intersections)
        rp_by_node[genome][bin] = recall_precision(genome, bin)
    end
    push!(counters, get_counts(rp_by_node, recalls, precisions))
    while length(rp_by_node) > 1 || only(keys(rp_by_node)).parent !== nothing
        rp_by_node = uprank_rp_by_node(rp_by_node)
        push!(counters, get_counts(rp_by_node, recalls, precisions))
    end
    counters
end

function get_counts(
    rp_by_node::Dict{Node, Dict{Bin, Tuple{Float64, Float64}}},
    recalls::Vector{Float64},
    precisions::Vector{Float64},
)
    # The current genome has a bin at the given recall with at most precision
    max_precision_at_recalls = zeros(length(recalls))
    counts = zeros(Int, (length(precisions), length(recalls)))
    for dict in values(rp_by_node)
        fill!(max_precision_at_recalls, zero(eltype(max_precision_at_recalls)))
        # For each (rec, prec) pair, find the corresponding recall, and set
        # the max precision in max_precision_at_recalls
        for (_, (recall, precision)) in dict
            i = searchsortedlast(recalls, recall)
            iszero(i) && continue
            max_precision_at_recalls[i] = max(precision, max_precision_at_recalls[i])
        end
        # If max_precision_at_recalls[i] == x, then all previous values in the
        # vector should be at least x. Set it to that if not already.
        last_prec = last(max_precision_at_recalls)
        for i in lastindex(max_precision_at_recalls)-1:-1:1
            last_prec = max(last_prec, max_precision_at_recalls[i])
            max_precision_at_recalls[i] = last_prec
        end
        # Now, for each (recall_index, precision) in pairs(max_precision_at_recalls),
        # find the corresponding precision_index. The bin is seen at that recall,
        # at all precision indices up to and including that index.
        # For the next recall, the precision index can be no larger than the previous
        # one
        max_precision_index = lastindex(precisions)
        for (recall_index, precision) in pairs(max_precision_at_recalls)
            precision_index = searchsortedlast(view(precisions, 1:max_precision_index), precision)
            iszero(precision_index) && break
            max_precision_index = min(precision_index, max_precision_index)
            counts[1:precision_index, recall_index] .+= 1
        end
    end
    counts
end

function uprank_rp_by_node(
    rp_by_node::Dict{Node, Dict{Bin, Tuple{Float64, Float64}}},
)
    result = empty(rp_by_node)
    for (child, child_dict) in rp_by_node
        parent = child.parent::Clade
        parent_dict = get!(valtype(result), result, parent)
        for (bin, (old_recall, old_prec)) in child_dict
            (new_recall, new_prec) = get(parent_dict, bin, (0.0, 0.0))
            parent_dict[bin] = (
                max(old_recall, new_recall),
                old_prec + new_prec
            )
        end
    end
    result
end