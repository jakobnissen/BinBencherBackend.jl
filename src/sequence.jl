struct Sequence
    name::String
    subject::String
    span::UnitRange{Int}

    function Sequence(name::AbstractString, subject::AbstractString, span::UnitRange{Int})
        isempty(span) && throw(ArgumentError("Cannot instantiate an empty sequence"))
        new(name, subject, span)
    end
end

function Base.isless(s::Sequence, s2::Sequence)
    isless((s.subject, first(s.span)), (s2.subject, first(s2.span)))
end

Base.length(x::Sequence) = length(x.span)
Base.hash(x::Sequence, h::UInt) = hash(x.name, h ⊻ UInt(24364341))
Base.:(==)(s1::Sequence, s2::Sequence) = s1.name == s2.name

function Sequence(name::AbstractString, len::Integer)
    s = String(name)
    Sequence(s, s, 1:Int(len))
end