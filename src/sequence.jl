struct Sequence
    name::String
    len::Int

    function Sequence(name::AbstractString, len::Integer)
        len < 1 && error("Cannot instantiate empty Sequence")
        new(name, len)
    end
end

Base.length(x::Sequence) = x.len
Base.hash(x::Sequence, h::UInt) = hash(x.name, h ⊻ UInt(24364341))
Base.:(==)(s1::Sequence, s2::Sequence) = s1.name == s2.name