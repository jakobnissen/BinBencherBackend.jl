@lazy mutable struct Source{G}
    const name::String
    const genome::G
    const length::Int
    const sequences::Vector{Tuple{Sequence, UnitRange{Int}}}
    @lazy n_covered::Int

    function Source(genome::G, name::AbstractString, length::Integer) where G
        length ≤ 0 && error("Source length must be at least 1")
        new{G}(String(name), genome, Int(length), Vector{Tuple{Sequence, UnitRange{Int}}}(), uninit)
    end
end

Base.:(==)(s1::Source, s2::Source) = s1.name == s2.name
Base.hash(x::Source, h::UInt) = hash(x.name, h ⊻ UInt(344509130945))

function add_sequence!(source::Source, seq::Sequence, span::UnitRange{Int})
    @assert !isempty(span)
    if !issubset(span, 1:source.length)
        error(
            lazy"Attempted to add sequence \"$(seq.name)\" to source \"$(source.name)\" ",
            lazy"at span $(span), but valid source indices are 1:$(source.length)"
        )
    end
    @isinit(source.n_covered) && error("Can't add sequence to source after calling finish! on it.")
    push!(source.sequences, (seq, span))
    source
end

function finish!(source::Source)
    @isinit(source.n_covered) && return source
    @init! source.n_covered = overlap_breadth(source.sequences, last)
    source
end

function overlap_breadth(v, by=identity)
    breadth = 0
    rightmost_end = 0
    for i in sort!(v; by=i -> first(by(i)), alg=QuickSort)
        span = by(i)::UnitRange
        breadth += max(last(span), rightmost_end) - max(first(span) - 1, rightmost_end)
        rightmost_end = max(rightmost_end, last(span))
    end
    breadth
end

Base.show(io::IO, x::Source) = print(io, "Source(", x.name, ", ", x.length, ')')
function Base.show(io::IO, ::MIME"text/plain", x::Source)
    if get(io, :compact, false)
        show(io, x)
    else
        print(io,
            "Source \"", x.name,
            "\"\ngenome:      ", x.genome,
            "\n  length:      ", x.length,
            "\n  n_covered:   ", x.n_covered,
            "\n  n_sequences: ", length(x.sequences)
        )
    end
end