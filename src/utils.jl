function tab_pairs(lines)
    lines |>
    imap(strip) |>
    ifilter(!isempty) |>
    imap() do line
        cu = codeunits(line)
        t1 = findfirst(isequal(UInt8('\t')), cu)
        t1 === nothing && error(lazy"No tabs in line $line")
        t2 = findnext(isequal(UInt8('\t')), cu, t1 + 1)
        t2 === nothing || error(lazy"More than two tab-sep fields in $line")
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
