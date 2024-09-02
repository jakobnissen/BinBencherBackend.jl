struct Unsafe end
const unsafe = Unsafe()

function tab_pairs(lines)
    enumerate(lines) |> imap() do (line_no, line)
        cu = codeunits(line)
        t1 = findfirst(isequal(UInt8('\t')), cu)
        t1 === nothing && error(lazy"No tabs in line number $line_no")
        t2 = findnext(isequal(UInt8('\t')), cu, t1 + 1)
        t2 === nothing || error(lazy"More than two tab-sep fields in line number $line_no")
        f1 = SubString(line, 1:prevind(line, t1))
        f2 = SubString(line, (t1 + 1):lastindex(line))
        (f1, f2)
    end
end

function binsplit_tab_pairs(t_pairs, sep::Union{Char, AbstractString})
    t_pairs |> imap() do (binname, seqname)
        p = findfirst(sep, seqname)
        p === nothing && error(lazy"Seperator $sep not found in seq name $seqname")
        before = SubString(seqname, 1:last(p))
        new_binname_ = string(before, binname)
        new_binname = SubString(new_binname_, 1:lastindex(new_binname_))
        (new_binname, seqname)
    end
end

function open_perhaps_gzipped(f::Function, path::String)
    if endswith(path, ".gz")
        stream = GzipDecompressorStream(open(path; lock=false))
        try
            f(stream)
        finally
            close(stream)
        end
    else
        open(f, path; lock=false)
    end
end

# This function exists in order to be able to use a Set
# to detect duplicates without hashing and indexing the Set
# twice.
function in!(s::AbstractSet, x)::Bool
    xT = convert(eltype(s), x)
    L = length(s)
    push!(s, xT)
    return length(s) == L
end
