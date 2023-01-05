# Breadth is the total number of covered basepairs.
mutable struct Genome
    name::String
    sources::Dict{String, Int}
    breadth::Int
    # This field is unassigned during construction, and later set.
    parent::Clade{Genome}

    function Genome(name::AbstractString)
        new(String(name), Dict{String, Int}(), 0)
    end
end
rank(::Genome) = 0

const Node = Union{Genome, Clade{Genome}}

# In the reference, a sequence can be assigned to either a genome or a clade,
# or a specific genome, subject and position within that subject
const Target = Union{Genome, Clade{Genome}, Tuple{Genome, String, UnitRange{Int}}}

function add_child!(c::Clade{Genome}, g::Node)
    children = c.children
    if g isa Genome
        push!(children::Vector{Genome}, g)
    else
        push!(children::Vector{Clade{Genome}}, g)
    end
    g.parent = c
    c.ngenomes += (g isa Genome ? 1 : g.ngenomes)
    c
end

Base.:(==)(g1::Genome, g2::Genome) = g1.name == g2.name
Base.hash(x::Genome, h::UInt) = hash(x.name, h ⊻ UInt(21323125590))

function add_source!(genome::Genome, source::String, len::Integer)
    if haskey(genome.sources, source)
        error(lazy"Genome $(genome.name) already have source $(source)")
    end
    genome.sources[source] = Int(len)
    genome.breadth += Int(len)
    genome
end

function remove_source!(genome::Genome, source::String)
    genome.breadth -= pop!(genome.sources, source)
    genome
end

Base.show(io::IO, x::Genome) = print(io, "Genome(", x.name, ')')
function Base.show(io::IO, ::MIME"text/plain", x::Genome)
    if get(io, :compact, false)
        show(io, x)
    else
        print(io,
            "Genome \"", x.name,
            "\"\n  Parent: ", '"', x.parent.name, '"',
            "\n  Breadth: ", x.breadth,
            "\n  Sources: ", length(x.sources)
        )
    end
end

"""
    mrca(a::Node, b::Node)::Node

Compute the most recent common ancestor (MRCA) of `a` and `b`.
"""
function mrca(a::Node, b::Node)
    a === b && return a
    ca = a isa Genome ? a.parent : a
    cb = b isa Genome ? b.parent : b
    (lo, hi) = ca.rank < cb.rank ? (ca, cb) : (cb, ca)
    while lo.rank < hi.rank
        lo = lo.parent::Clade
    end
    while lo !== hi
        lo = lo.parent::Clade
        hi = hi.parent::Clade
    end
    lo
end