@lazy mutable struct Source{G}
    const name::String
    const genome::G
    const length::Int
    const sequences::Vector{Tuple{Sequence, UnitRange{Int}}}
    @lazy assembly_size::Int

    function Source(genome::G, name::AbstractString, length::Integer) where {G}
        length ≤ 0 && error("Source length must be at least 1")
        new{G}(
            String(name),
            genome,
            Int(length),
            Vector{Tuple{Sequence, UnitRange{Int}}}(),
            uninit,
        )
    end
end

"""
    Source{Genome}(g::Genome, name::AbstractString, length::Integer)

Sources are the "ground truth" sequences that the binning attempts to recreate.
For example, the assembled contigs of the reference genome (typically full, closed circular
contigs) as found in NCBI or elsewhere are each `Source`s.
Many `Genome`s only contain a single `Source` namely its full assembled genome.
Each `Source` has a single parent `Genome`, and a unique name which identifies it.

`Source`s have zero or more mapping `Sequence`s, that each map to the `Source` at a given
span given by a `UnitRange{Int}`.

`Source`s have an _assembly size_, which is the number of base pairs where any sequence map
to.
"""
Source

Base.:(==)(s1::Source, s2::Source) = s1.name == s2.name
Base.hash(x::Source, h::UInt) = hash(x.name, h ⊻ UInt(344509130945))

function add_sequence!(source::Source, seq::Sequence, span::UnitRange{Int})
    @assert !isempty(span)
    if !issubset(span, 1:(source.length))
        error(
            lazy"Attempted to add sequence \"$(seq.name)\" to source \"$(source.name)\" ",
            lazy"at span $(span), but valid source indices are 1:$(source.length)",
        )
    end
    @isinit(source.assembly_size) &&
        error("Can't add sequence to source after calling finish! on it.")
    push!(source.sequences, (seq, span))
    source
end

function finish!(source::Source)
    @isinit(source.assembly_size) && return source
    asm_size = assembly_size(source.sequences, last)
    @assert asm_size ≤ source.length
    @init! source.assembly_size = asm_size
    source
end

# v: Vector of X, where by(X) isa UnitRange
function assembly_size(v::Vector, by=identity)
    size = 0
    rightmost_end = 0
    for i in sort!(v; by=i -> first(by(i)), alg=QuickSort)
        span = by(i)::UnitRange
        size += max(last(span), rightmost_end) - max(first(span) - 1, rightmost_end)
        rightmost_end = max(rightmost_end, last(span))
    end
    size
end

Base.show(io::IO, x::Source) = print(io, "Source(", x.name, ", ", x.length, ')')
function Base.show(io::IO, ::MIME"text/plain", x::Source)
    if get(io, :compact, false)
        show(io, x)
    else
        print(
            io,
            "Source \"",
            x.name,
            "\"\ngenome:          ",
            x.genome,
            "\n  Length:        ",
            x.length,
            "\n  Assembly size: ",
            x.assembly_size,
            "\n  Sequences:     ",
            length(x.sequences),
        )
    end
end
