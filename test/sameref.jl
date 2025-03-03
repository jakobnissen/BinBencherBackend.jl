using Test

macro ret(ex)
    return quote
        local res = $(esc(ex))
        isnothing(res) ? nothing : return res
    end
end

function test_is_same(a::Reference, b::Reference)::Union{Nothing, String}
    # Genomes
    ga = sort!(collect(genomes(a)); by = i -> i.name)
    gb = sort!(collect(genomes(b)); by = i -> i.name)
    length(ga) == length(gb) || return "Different number of genomes"
    for (gga, ggb) in zip(ga, gb)
        @ret test_is_same(gga, ggb)
    end

    # Clades
    length(a.clades) == length(b.clades) || return "Different ranks"
    for (n, (i, j)) in enumerate(zip(a.clades, b.clades))
        length(i) == length(j) || return "At rank $(n), different lengths"
        i = sort!(copy(i); by = i -> i.name)
        j = sort!(copy(j); by = i -> i.name)
        for (q, p) in zip(i, j)
            @ret test_is_same(p, q)
        end
    end

    # Sequences and targets
    length(a.targets) == length(b.targets) || return "Different num seqs"
    @assert length(a.target_index_by_name) == length(b.target_index_by_name)
    for (name, i) in a.target_index_by_name
        j =
            @something get(b.target_index_by_name, name, nothing) return "Different sequences"
        (sa, atargs) = a.targets[i]
        (sb, btargs) = b.targets[j]
        @ret test_is_same(sa, sb)
        atargs = sort!(copy(atargs); by = i -> (i[1].name, i[2][1]))
        btargs = sort!(copy(btargs); by = i -> (i[1].name, i[2][1]))
        length(atargs) == length(btargs) || return "Different targets"
        for ((ta1, (ta2, ta3)), (tb1, (tb2, tb3))) in zip(atargs, btargs)
            ta2 == tb2 || return "Different targets"
            ta3 == tb3 || return "Different targets"
            ta1.name == tb1.name || return "Different targets"
        end
    end
    return nothing
end

function test_is_same(a::Genome, b::Genome)::Union{Nothing, String}
    a.name == b.name || return "Different genome name"
    a.flags == b.flags || return "Different genome flags"
    asrc = sort!(collect(a.sources); by = i -> i.name)
    bsrc = sort!(collect(b.sources); by = i -> i.name)
    length(asrc) == length(bsrc) || return "Different source lengths"
    for (i, j) in zip(asrc, bsrc)
        i.name == j.name || return "Different sources"
        i.length == j.length || return "Different source lengths"
    end
    return nothing
end

function test_is_same(a::Clade{Genome}, b::Clade{Genome})::Union{Nothing, String}
    a.name == b.name || return "Different clade name"
    a.rank == b.rank || return "Different clade orders (should not happen?)"
    agns = sort!(collect(a.children); by = i -> i.name)
    bgns = sort!(collect(b.children); by = i -> i.name)
    length(agns) == length(bgns) || return "Different num children"
    for (i, j) in zip(agns, bgns)
        i.name == j.name || return "Different children"
    end
    return nothing
end

function test_is_same(a::Sequence, b::Sequence)::Union{Nothing, String}
    return if a.name == b.name && length(a) == length(b)
        nothing
    else
        "Different sequences"
    end
end
