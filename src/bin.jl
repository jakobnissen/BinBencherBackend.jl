struct Bin
    name::String
    sequences::Vector{Sequence}
    intersections::Dict{Genome, Int}
    breadth::Int

    function Bin(name::AbstractString, seqs, genomeof::Dict{Sequence, Genome})
        sequences = vector(seqs)
        intersections = Dict{Genome, Int}()
        by_source = Dict{Tuple{Genome, String}, Vector{Sequence}}()
        for seq in sequences
            push!(get!(valtype(by_source), by_source, (genomeof[seq], seq.subject)), seq)
        end
        for ((genome, _), source_seqs) in by_source
            intersections[genome] = get(intersections, genome, 0) + overlap_breadth(source_seqs)
        end
        new(String(name), sequences, intersections, sum(values(intersections), init=0))
    end
end

function overlap_breadth(sequences::Vector{Sequence})
    sort!(sequences, by=x -> first(x.span))
    result = 0
    rightmost_end = 0
    for sequence in sequences
        result += max(last(sequence.span), rightmost_end) -
            max(first(sequence.span), rightmost_end)
        rightmost_end = max(rightmost_end, last(sequence.span))
    end
    result
end

nseqs(x::Bin) = length(x.sequences)

Base.show(io::IO, x::Bin) = print(io, "Bin(", x.name, ')')
function Base.show(io::IO, ::MIME"text/plain", x::Bin)
    if get(io, :compact, false)
        show(io, x)
    else
        n_int = length(x.intersections)
        print(io,
            "Bin \"", x.name,
            "\"\n  Sequences: ", nseqs(x),
            "\n  Breadth:   ", x.breadth,
            "\n  Intersecting ", n_int, " genome", (isone(n_int) ? "" : "s"),
        )
    end
end

function confusion_matrix(genome::Genome, bin::Bin)
    breadth = Int(bin.breadth)
    tp = Int(get(bin.intersections, genome, 0))
    fp = breadth - tp
    fn = Int(genome.breadth) - tp
    (tp, fp, fn)
end

function recall_precision(genome::Genome, bin::Bin)
    (tp, fp, fn) = confusion_matrix(genome, bin)
    recall = tp / (tp + fn)
    precision = tp / (tp + fp)
    (recall, precision)
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