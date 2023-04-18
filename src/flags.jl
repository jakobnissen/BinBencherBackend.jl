baremodule Flags
using Base: @enum
@enum Flag::UInt8 organism virus plasmid

export Flag
end

using .Flags

@doc """
        Flag

    A flag is a boolean associated to a `Genome`, stored in a `Flags` object.
    A flag may be e.g. `Flag.organism`, signaling that the genome is known to be
    an organism.

    See also: [`FlagSet`](@ref), [`Genome`](@ref)
    """
Flags.Flag

"""
    FlagSet <: AbstractSet{Flag}

Flags are compact sets of `Flag` associated to a Genome.
You can construct them from an iterable of `Flag`, e.g. a 1-element tuple.
`FlagSet` support most set operations efficiently.

# Examples
```jldoctest
julia> flags = FlagSet((Flags.organism, Flags.virus));

julia> Flags.virus in flags
true

julia> isdisjoint(flags, FlagSet((Flags.organism,)))
false
```

See also: [`Flag`](@ref), [`Genome`](@ref)
"""
struct FlagSet <: AbstractSet{Flag}
    # TODO: Use fewer bits?
    x::UInt64

    FlagSet(x::UInt64) = new(x)
end

FlagSet() = FlagSet(zero(UInt64))

function FlagSet(itr)
    result = FlagSet()
    for i in itr
        result = push(result, convert(Flag, i))
    end
    result
end

function Base.iterate(x::FlagSet, state::UInt64=x.x)
    iszero(state) ? nothing :
    (reinterpret(Flag, trailing_zeros(state) % UInt8), state & (state - 1))
end

push(x::FlagSet, y::Flag) = FlagSet(x.x | (UInt64(1) << (reinterpret(UInt8, y) & 63)))
Base.in(flag::Flag, x::FlagSet) = isodd(x.x >>> reinterpret(UInt8, flag))
Base.length(x::FlagSet) = count_ones(x.x)
Base.isempty(x::FlagSet) = iszero(x.x)
Base.hash(x::FlagSet, h::UInt) = hash(x.x, h ⊻ (0x116133601de11a93 % UInt))
Base.union(x::FlagSet, y::FlagSet) = FlagSet(x.x | y.x)
Base.intersect(x::FlagSet, y::FlagSet) = FlagSet(x.x & y.x)
Base.setdiff(x::FlagSet, y::FlagSet) = FlagSet(x.x & ~y.x)
Base.issubset(x::FlagSet, y::FlagSet) = isempty(setdiff(x, y))
Base.isdisjoint(x::FlagSet, y::FlagSet) = isempty(x ∩ y)
