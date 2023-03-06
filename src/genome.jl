@lazy mutable struct Genome
    const name::String
    const sources::Set{Source{Genome}}
    @lazy parent::Clade{Genome}
    @lazy breadth::Int

    function Genome(name::AbstractString)
        new(String(name), Set{Source{Genome}}(), uninit, uninit)
    end
end

const Target = Tuple{Source{Genome}, UnitRange{Int}}
const Node = Union{Genome, Clade{Genome}}

function add_child!(c::Clade{Genome}, g::Node)
    children = c.children
    if g isa Genome
        push!(children::Vector{Genome}, g)
        @init! g.parent = c
    else
        push!(children::Vector{Clade{Genome}}, g)
        g.parent = c
    end
    c.ngenomes += (g isa Genome ? 1 : g.ngenomes)
    c
end

Base.:(==)(g1::Genome, g2::Genome) = g1.name == g2.name
Base.hash(x::Genome, h::UInt) = hash(x.name, h ⊻ UInt(21323125590))

function add_source!(genome::Genome, name::AbstractString, length::Integer)
    @isinit(genome.breadth) && error("Can't add source to genome after calling finish! on it.")
    source = Source(genome, name, length)
    in(source, genome.sources) && error(lazy"Genome $(genome.name) already have source $(source.name)")
    push!(genome.sources, source)
    genome
end

function finish!(genome::Genome)
    @isinit(genome.breadth) && return genome
    @isinit(genome.parent) || error(
        lazy"finish! called on genome \"$(genome.name)\" without assigned parent."
    )
    for source in genome.sources
        @isinit(source.n_covered) || finish!(source)
    end
    @init! genome.breadth = sum(i -> i.n_covered, genome.sources)
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