struct Bin
    name::String
    sequences::Vector{Sequence}
    # For sequences that can be assigned, we compute the number of bases in the
    # given genome covered
    bases_covered::Dict{Genome, Int}
    # This is simply the sum of lengths of each sequence that is assigned to a Node
    # Only the lowest clade is put here, e.g. is a seq is assigned to E. coli, it is NOT
    # also assigned to Escherischia
    length_sums::Dict{Node, Int}
    # Sum of values of above dict
    length_sum::Int

    function Bin(name::AbstractString, seqs, targets::Dict{Sequence, Target})
        sequences = vector(seqs)::Vector{Sequence}
        bases_covered = Dict{Genome, Int}()
        length_sums = Dict{Node, Int}()
        by_source = Dict{Tuple{Genome, String}, Vector{UnitRange{Int}}}()
        length_sum = 0
        for seq in sequences
            length_sum += length(seq)
            target = targets[seq]
            if target isa Node
                length_sums[target] = get(length_sums, target, 0) + length(seq)
            else
                (genome, subject, span) = target
                length_sums[genome] = get(length_sums, genome, 0) + length(seq)
                push!(get!(valtype(by_source), by_source, (genome, subject)), span)
            end
        end
        for ((genome, _), source_seqs) in by_source
            bases_covered[genome] = get(bases_covered, genome, 0) + overlap_breadth!(source_seqs)
        end
        new(String(name), sequences, bases_covered, length_sums, length_sum)
    end
end

function overlap_breadth!(sequences::Vector{<:UnitRange{<:Integer}})
    sort!(sequences; by=first, alg=QuickSort)
    result = 0
    rightmost_end = 0
    for span in sequences
        result += max(last(span), rightmost_end) - max(first(span), rightmost_end)
        rightmost_end = max(rightmost_end, last(span))
    end
    result
end

nseqs(x::Bin) = length(x.sequences)

Base.show(io::IO, x::Bin) = print(io, "Bin(", x.name, ')')
function Base.show(io::IO, ::MIME"text/plain", x::Bin)
    if get(io, :compact, false)
        show(io, x)
    else
        n_int = length(x.bases_covered)
        print(io,
            "Bin \"", x.name,
            "\"\n  Sequences: ", nseqs(x),
            "\n  Length sum:   ", x.length_sum,
            "\n  Covering ", n_int, " genome", (isone(n_int) ? "" : "s"),
        )
    end
end

get_recall(genome::Genome, bin::Bin) = get(bin.bases_covered, genome, 0) / genome.breadth
get_precision(genome::Genome, bin::Bin) = get(bin.length_sums, genome, 0) / bin.length_sum

function fscore(genome::Genome, bin::Bin, b::AbstractFloat)
    recall = get_recall(genome, bin)
    precision = get_precision(genome, bin)
    # Some people say the Fscore is undefined in this case.
    # We define it to be 0.0
    if iszero(recall + precision)
        return 0.0
    end
    (1 + b^2) * (recall*precision) / ((b^2*precision)+recall)
end
f1(genome::Genome, bin::Bin) = fscore(genome, bin, 1.0)