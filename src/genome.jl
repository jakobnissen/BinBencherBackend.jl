@lazy mutable struct Genome
    const name::String
    const sources::Set{Source{Genome}}
    const flags::FlagSet
    @lazy parent::Clade{Genome}
    @lazy genome_size::Int
    @lazy assembly_size::Int

    function Genome(name::AbstractString, flags::FlagSet)
        new(
            check_valid_identifier(String(name)),
            Set{Source{Genome}}(),
            flags,
            uninit,
            uninit,
            uninit,
        )
    end
end
Genome(name::AbstractString) = Genome(name, FlagSet())

"""
    Genome(name::AbstractString [flags::FlagSet])

`Genome`s represent individual target genomes (organisms, plasmids, viruses etc),
and are conceptually the lowest-level clade that can be reconstructed.
`Genome`s contain one or more `Source`s, and belong to a single parent `Clade`.
They are identified uniquely among genomes by their name.

A genome have a _genome size_, which is the sum of the length of all its sources.
We consider this to be the true size of the biological genome (assuming its full
sequence is contained in its sources), as well as an _assembly size_, which represent
the sum of the assembly sizes of each source.

See also: [`Clade`](@ref), [`Source`](@ref), [`mrca`](@ref)

# Examples
```jldoctest
julia> gA, gB, gC = collect(ref.genomes);


julia> flags(gA)
FlagSet with 1 element:
  BinBencherBackend.Flags.organism

julia> mrca(gA, gB)
Species "D", 2 genomes
├─ Genome(gA)
└─ Genome(gB)
```
"""
Genome

const Target = Tuple{Source{Genome}, Tuple{Int, Int}}
const Node = Union{Genome, Clade{Genome}}

"""
    flags(g::Genome)::FlagSet

Returns the `Flag`s of the `Genome` as a `FlagSet`.

See also: [`Flag`](@ref), [`FlagSet`](@ref)

# Example
```jldoctest
julia> flags(genome)
FlagSet with 1 element:
  BinBencherBackend.Flags.organism
```
"""
flags(g::Genome) = g.flags

"""
    is_organism(g::Genome)::Bool

Check if `g` is known to be an organism.

# Example
```jldoctest
julia> is_organism(genome)
true
``` 
"""
is_organism(g::Genome) = Flags.organism ∈ flags(g)

"""
    is_virus(g::Genome)::Bool

Check if `g` is known to be a virus.

# Example
```jldoctest
julia> is_virus(genome)
false
``` 
"""
is_virus(g::Genome) = Flags.virus ∈ flags(g)

"""
    is_plasmid(g::Genome)::Bool

Check if `g` is known to be a plasmid.

# Example
```jldoctest
julia> is_plasmid(genome)
false
``` 
"""
is_plasmid(g::Genome) = Flags.plasmid ∈ flags(g)

function add_child!(c::Clade{Genome}, g::Node)::Clade
    children = c.children
    if g isa Genome
        push!(children::Vector{Genome}, g)
        @init! g.parent = c
    else
        push!(children::Vector{Clade{Genome}}, g)
        g.parent = c
    end
    c.ngenomes += (g isa Genome ? 1 : g.ngenomes)
    return c
end

"Delete a child from the clade tree."
function recursively_delete_child!(child::T) where {T <: Node}
    parent = child.parent
    # If the child is the top level clade, do nothing, as we delete from the bottom up
    parent === nothing && return nothing
    children = parent.children::Vector{T}
    deleteat!(children, findfirst(i -> i === child, children)::Integer)
    # If we delete a genome, we remove that genome from all ancestors count
    if child isa Genome
        p = parent
        while p !== nothing
            p.ngenomes -= 1
            p = p.parent
        end
    end
    # If the clade now has no children, we can delete the clade
    return isempty(children) && recursively_delete_child!(parent)
end

Base.:(==)(g1::Genome, g2::Genome) = g1.name == g2.name
Base.hash(x::Genome, h::UInt) = hash(x.name, h ⊻ UInt(21323125590))

function add_source!(genome::Genome, name::AbstractString, length::Integer)
    @isinit(genome.genome_size) &&
        error("Can't add source to genome after calling finish! on it.")
    source = Source(genome, name, length)
    in(source, genome.sources) &&
        error(lazy"Genome $(genome.name) already have source $(source.name)")
    push!(genome.sources, source)
    return genome
end

function finish!(genome::Genome, scratch::Vector{Tuple{Int, Int}})
    @isinit(genome.genome_size) && return genome
    @isinit(genome.parent) ||
        error(lazy"finish! called on genome \"$(genome.name)\" without assigned parent.")
    for source in genome.sources
        @isinit(source.assembly_size) || finish!(source, scratch)
    end
    @init! genome.genome_size = sum(i -> i.length, genome.sources; init = 0)
    @init! genome.assembly_size = sum(i -> i.assembly_size, genome.sources; init = 0)
    return genome
end

Base.show(io::IO, x::Genome) = print(io, "Genome(", x.name, ')')
function Base.show(io::IO, ::MIME"text/plain", x::Genome)
    return if get(io, :compact, false)
        show(io, x)
    else
        asm = (x.assembly_size / x.genome_size) * 100
        print(
            io,
            "Genome \"",
            x.name,
            "\"\n  Parent:        ",
            '"',
            x.parent.name,
            '"',
            "\n  Genome size:   ",
            x.genome_size,
            "\n  Assembly size: ",
            x.assembly_size,
            " (",
            round(asm; digits = 1),
            " %)",
            "\n  Sources:       ",
            length(x.sources),
            "\n  Flags:         ",
            Int(x.flags.x),
            " (",
            join(x.flags, ", "),
            ')',
        )
    end
end

"""
    mrca(a::Node, b::Node)::Node

Compute the most recent common ancestor (MRCA) of `a` and `b`.
"""
function mrca(a::Node, b::Node)::Node
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
    return lo
end

"""
    descends_from(child, ancestor) -> Bool

Check if `child` descends from `ancestor`.
Also returns `true` if `child` and `ancestor` is identical.
Boths arguments must be a `Union{Genome, Clade}`.

# Examples
```
julia> descends_from(genome, genome.parent)
true

julia> descends_from(genome, genome)
true

julia> descends_from(genome.parent, genome)
false
```
"""
function descends_from(child::Node, ancestor::Node)::Bool
    child === ancestor && return true
    ancestor isa Genome && return false
    for anc in ancestors(child)
        anc.rank > ancestor.rank && return false
        anc === ancestor && return true
    end
    return false
end

# The design here is awkward: Do we store the next clade to emit,
# or the child? If the child, then it's Union{Genome, Clade}, if
# the next clade, it's Union{Clade, Nothing}.
# We instead store first clade to emit, and a bool telling if it's
# empty (i.e. it began at the top clade)
struct AncestorIterator
    start::Clade{Genome}
    empty::Bool
end

# We could actually precompute the length, if we got a reference
# to the Reference object. But we don't have that object at hand,
# and this is unlikely to be performance critical
Base.IteratorSize(::Type{AncestorIterator}) = Base.SizeUnknown()

Base.eltype(::Type{AncestorIterator}) = Clade{Genome}

"""
    ancestors(x::Union{Genome, Clade})

Return an iterator of the ancestors of `x`. The first element will be `x`'s
parent, the next element its grandparent etc.
All elements of the iterator are `Clade{Genome}`.

# Examples
```jldoctest
julia> count(Returns(true), ancestors(genome))
2

julia> println([i.name for i in ancestors(genome)])
["D", "F"]
```
"""
ancestors(child::Genome) = AncestorIterator(child.parent, false)

function ancestors(child::Clade{Genome})
    next = child.parent
    return isnothing(next) ? AncestorIterator(child, true) : AncestorIterator(next, false)
end

# The state is essentially a Union{Clade, Nothing}, but I want the state
# to be a concrete type
function Base.iterate(it::AncestorIterator, state = (it.start, it.empty))
    (next, is_empty) = state
    is_empty && return nothing
    parent = next.parent
    return (next, (isnothing(parent) ? (next, true) : (parent, false)))
end
