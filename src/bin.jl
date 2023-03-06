struct Bin
    name::String
    sequences::Vector{Sequence}
    # True positives (bases of genome covered)
    # False positives (sum of lengths of sequences not mapping to genome)
    intersections::Dict{Genome, Tuple{Int, Int}}
    # Sum of lengths of sequences, cached for efficiency
    breadth::Int
end

function Bin(
    name::AbstractString,
    sequences::Vector{Sequence},
    reference_targets::Dict{String, Tuple{Sequence, Vector{Target}}}
)
    breadth = sum(i -> i.length, sequences, init=0)
    spans = Dict{Source{Genome}, Vector{UnitRange{Int}}}()
    by_genome = Dict{Genome, Set{Sequence}}()
    for seq in sequences, (source, span) in last(reference_targets[seq.name])
        push!(get!(valtype(spans), spans, source), span)
        push!(get!(valtype(by_genome), by_genome, source.genome), seq)
    end
    intersections = Dict{Genome, Tuple{Int, Int}}()
    default = (0, 0)
    for (source, ranges) in spans
        (tp, fp) = get(intersections, source.genome, default)
        intersections[source.genome] = (tp + overlap_breadth(ranges), fp)
    end
    for (genome, seqs) in by_genome
        intersections[genome] = (intersections[genome][1], breadth - sum(i -> i.length, seqs, init=0))
    end
    Bin(String(name), sequences, intersections, breadth)
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
    (tp, fp) = get(bin.intersections, genome, (0, bin.breadth))
    fn = genome.breadth - tp
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