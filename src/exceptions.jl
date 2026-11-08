"""
    BinBencherException <: Exception

Abstract supertype of exceptions thrown by BinBencherBackend when its input files
are invalid.

Exceptions of this type are thrown where the problem is detected, so any
`BinBencherException` signals a problem with the input, and not a bug in
BinBencherBackend.
Misuse of the API (e.g. passing an invalid recall threshold) throws an `ArgumentError`
instead.

See also: [`ReferenceError`](@ref), [`BinningError`](@ref)
"""
abstract type BinBencherException <: Exception end

"""
    ReferenceError(msg::String) <: BinBencherException

Exception thrown when loading a `Reference` from an invalid reference file,
e.g. malformed JSON, duplicate sequences, or an invalid taxonomy.

See also: [`BinBencherException`](@ref), [`Reference`](@ref)
"""
struct ReferenceError <: BinBencherException
    msg::String
end

function Base.showerror(io::IO, e::ReferenceError)
    return print(io, "ReferenceError: ", e.msg)
end

"""
    BinningError(msg::String, line::Union{Int, Nothing}) <: BinBencherException

Exception thrown when loading a `Binning` from an invalid cluster file,
e.g. a missing header, a malformed line, or a sequence not in the reference.
`line` is the 1-based line number in the cluster file where the problem was found,
or `nothing` if the problem is not tied to a single line.

See also: [`BinBencherException`](@ref), [`Binning`](@ref)
"""
struct BinningError <: BinBencherException
    msg::String
    line::Union{Int, Nothing}
end

function Base.showerror(io::IO, e::BinningError)
    print(io, "BinningError: ")
    e.line === nothing || print(io, "line ", e.line, ": ")
    return print(io, e.msg)
end
