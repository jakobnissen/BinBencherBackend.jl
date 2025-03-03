@lazy mutable struct Source{G}
    const name::String
    const genome::G
    const length::Int
    const sequences::Vector{Tuple{Sequence, Tuple{Int, Int}}}
    @lazy assembly_size::Int
    @lazy total_bp::Int

    function Source(genome::G, name::AbstractString, length::Integer) where {G}
        length ≤ 0 && error("Source length must be at least 1")
        new{G}(
            check_valid_identifier(String(name)),
            genome,
            Int(length),
            Vector{Tuple{Sequence, Tuple{Int, Int}}}(),
            uninit,
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
span given by a 2-tuple `Tuple{Int, Int}`.

`Source`s have an _assembly size_, which is the number of base pairs where any sequence map
to.
"""
Source

Base.:(==)(s1::Source, s2::Source) = s1.name == s2.name
Base.hash(x::Source, h::UInt) = hash(x.name, h ⊻ UInt(344509130945))

function add_sequence!(source::Source, seq::Sequence, span::Tuple{Int, Int})
    (small, big) = minmax(span...)
    (small < 1 || big > source.length) && error(
        lazy"Attempted to add sequence \"$(seq.name)\" to source \"$(source.name)\" ",
        lazy"at span $(first(span)):$(last(span)), but valid source indices are 1:$(source.length)",
    )
    @isinit(source.assembly_size) &&
        error("Can't add sequence to source after calling finish! on it.")
    push!(source.sequences, (seq, span))
    return source
end

function finish!(source::Source, scratch::Vector{Tuple{Int, Int}})
    @isinit(source.assembly_size) && return source
    (asm_size, total_bp) = assembly_size!(last, scratch, source.sequences, source.length)
    @assert asm_size ≤ source.length
    @init! source.assembly_size = asm_size
    @init! source.total_bp = total_bp
    return source
end

"""Compute -> (breadth, total_bp), where breadth is the number of positions in `v` covered at least once,
and total_bp the sum of the lengths of the sequences.
`v` must be a `Vector` such that `all(by(i) isa Tuple{Integer, Integer} for i in v)`.
The `scratch` input is mutated.
"""
function assembly_size!(
        by::Function,
        scratch::Vector{Tuple{Int, Int}},
        v::Vector, # v: Vector of X, where by(X) isa Tuple{Integer, Integer}
        source_len::Int,
    )::Tuple{Integer, Integer}
    # First pass: Convert elements into Tuple{Int, Int}, and count the number
    # of circularly mapping spans (i.e. spans mapping from the end of the source
    # to the beginning)
    n_circular_mappings = 0
    resize!(scratch, length(v))
    @inbounds for i in eachindex(scratch, v)
        (start_, stop_) = by(v[i])::Tuple{Integer, Integer}
        (start, stop) = (Int(start_), Int(stop_))
        n_circular_mappings += start > stop
        scratch[i] = (start, stop)
    end
    # If we have any circular mappings, then these needs to be broken up into
    # two non-circular spans.
    # This is probably rare, so this function is written to be significantly
    # faster when this branch is not taken
    if !iszero(n_circular_mappings)
        old_size = length(scratch)
        new_size = old_size + n_circular_mappings
        resize!(scratch, new_size)
        # We write the extra split circular spans from the end of the vector,
        # to avoid overwriting elements that we are currently reading
        written_extras = 0
        @inbounds for i in 1:old_size
            (start, stop) = scratch[i]
            if start > stop
                scratch[i] = (1, stop)
                scratch[new_size - written_extras] = (start, source_len)
                written_extras += 1
            else
                scratch[i] = (start, stop)
            end
        end
        # @assert written_extras == n_circular_mappings # (should always hold)
    end
    # Now we know we have a Vector{Tuple{Int, Int}} with no circular mappings,
    # so we can compute the assembly size
    size = total_bp = rightmost_end = 0
    for span in sort!(scratch; by = first, alg = QuickSort)
        size += max(last(span), rightmost_end) - max(first(span) - 1, rightmost_end)
        total_bp += last(span) - first(span) + 1
        rightmost_end = max(rightmost_end, last(span))
    end
    return (size, total_bp)
end

Base.show(io::IO, x::Source) = print(io, "Source(", x.name, ", ", x.length, ')')
function Base.show(io::IO, ::MIME"text/plain", x::Source)
    return if get(io, :compact, false)
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
