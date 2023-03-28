@lazy mutable struct Genome
    const name::String
    const sources::Set{Source{Genome}}
    const flags::FlagSet
    @lazy parent::Clade{Genome}
    @lazy genome_size::Int
    @lazy assembly_size::Int

    function Genome(name::AbstractString, flags::FlagSet)
        new(String(name), Set{Source{Genome}}(), flags, uninit, uninit, uninit)
    end
end
Genome(name::AbstractString) = Genome(name, FlagSet())

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

function recursively_delete_child!(child::T) where {T <: Node}
    parent = child.parent
    parent === nothing && return nothing
    children = parent.children::Vector{T}
    deleteat!(children, findfirst(i -> i === child, children)::Integer)
    parent.ngenomes -= 1
    isempty(children) && recursively_delete_child!(parent)
end

Base.:(==)(g1::Genome, g2::Genome) = g1.name == g2.name
Base.hash(x::Genome, h::UInt) = hash(x.name, h ⊻ UInt(21323125590))

function add_source!(genome::Genome, name::AbstractString, length::Integer)
    @isinit(genome.genome_size) && error("Can't add source to genome after calling finish! on it.")
    source = Source(genome, name, length)
    in(source, genome.sources) && error(lazy"Genome $(genome.name) already have source $(source.name)")
    push!(genome.sources, source)
    genome
end

function finish!(genome::Genome)
    @isinit(genome.genome_size) && return genome
    @isinit(genome.parent) || error(
        lazy"finish! called on genome \"$(genome.name)\" without assigned parent."
    )
    for source in genome.sources
        @isinit(source.assembly_size) || finish!(source)
    end
    @init! genome.genome_size = sum(i -> i.length, genome.sources; init=0)
    @init! genome.assembly_size = sum(i -> i.assembly_size, genome.sources; init=0)
    genome
end

Base.show(io::IO, x::Genome) = print(io, "Genome(", x.name, ')')
function Base.show(io::IO, ::MIME"text/plain", x::Genome)
    if get(io, :compact, false)
        show(io, x)
    else
        asm = (x.assembly_size / x.genome_size) * 100
        print(io,
            "Genome \"", x.name,
            "\"\n  Parent:        ", '"', x.parent.name, '"',
            "\n  Genome size:   ", x.genome_size,
            "\n  Assembly size: ", x.assembly_size, " (", round(asm, digits=1), " %)",
            "\n  Sources:       ", length(x.sources),
            "\n  Flags:         ", Int(x.flags.x), " (", join(x.flags, ", "), ')'
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