struct Sequence
    name::String
    length::Int

    function Sequence(name::AbstractString, length::Integer)
        str = String(name)
        if isempty(str) || isspace(first(str)) || isspace(last(str))
            error(lazy"Sequence name \"$(str)\" cannot be empty or have leading or trailing whitespace")
        end
        length < 1 && throw(ArgumentError("Cannot instantiate an empty sequence"))
        new(str, Int(length))
    end
end

Base.length(x::Sequence) = x.length
Base.hash(x::Sequence, h::UInt) = hash(x.name, h ⊻ UInt(24364341))
Base.:(==)(s1::Sequence, s2::Sequence) = s1.name == s2.name
