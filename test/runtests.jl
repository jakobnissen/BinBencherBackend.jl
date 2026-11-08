using Test
using BinBencherBackend:
    BinBencherBackend,
    FlagSet,
    Flag,
    Flags,
    Reference,
    Binning,
    Bin,
    Genome,
    Clade,
    Source,
    Sequence,
    ancestors,
    descends_from,
    is_organism,
    is_virus,
    mrca,
    gold_standard,
    top_clade,
    genomes,
    f1,
    recall_precision,
    n_recovered,
    n_recovered_bounds,
    n_passing_bins,
    n_passing_bins_bounds,
    n_seqs,
    n_bins,
    subset,
    subset!,
    is_disjoint,
    adjusted_rand_index,
    BinBencherException,
    ReferenceError,
    BinningError

using CodecZlib: GzipCompressor
using JSON: JSON
using JSONSchema: JSONSchema

include("sameref.jl")

const DIR = joinpath(dirname(dirname(pathof(BinBencherBackend))), "files")
REF_PATH = joinpath(DIR, "ref.json")
REF_STR = read(REF_PATH, String)
CLUSTERS_PATH = joinpath(DIR, "clusters.tsv")
CLUSTERS_STR = read(CLUSTERS_PATH, String)

@assert isdir(DIR)
ngenomes(ref) = length(genomes(ref))

@testset "Misc" begin
    @test isnan(f1(0.0, 0.0))
    @test !isnan(f1(1.0e-6, 1.0e-6))
end

@testset "Flags" begin
    empt = FlagSet()
    a = FlagSet([Flags.organism, Flags.plasmid])
    b = FlagSet([Flags.virus])
    c = FlagSet([Flags.plasmid, Flags.virus])
    d = FlagSet([Flags.organism, Flags.plasmid, Flags.plasmid])

    @test only(b) in c
    @test Set(d) == Set([Flags.organism, Flags.plasmid, Flags.plasmid])
    @test only(b) == Flags.virus
    @test_throws Exception only(a)
    @test_throws Exception only(d)

    @test tryparse(Flag, "oRgANisM") == Flags.organism
    @test tryparse(Flag, "banana") === nothing
    @test FlagSet((tryparse(Flag, "virus"), tryparse(Flag, "organism"))) ==
        FlagSet([Flags.organism, Flags.virus])

    flagsets = [empt, a, b, c, d]
    all_pairs = [(i, j) for i in flagsets for j in flagsets]

    @test all(length(i) == length(Set(i)) for i in flagsets)
    @test all(isempty(i) == isempty(Set(i)) for i in flagsets)

    for f in [issubset, isdisjoint]
        @test all(all_pairs) do (a, b)
            f(a, b) == f(Set(a), Set(b))
        end
    end

    for f in [union, intersect, setdiff]
        @test all(all_pairs) do (a, b)
            Set(f(a, b)) == f(Set(a), Set(b))
        end
    end
end

@testset "Construction" begin
    global ref = Reference(IOBuffer(REF_STR))
    global binning = Binning(IOBuffer(CLUSTERS_STR), ref)
    global bins = sort!(collect(binning.bins); by = i -> i.name)

    @test ref isa Reference
    @test binning isa Binning
    @test bins isa Vector{Bin}
end

function test_is_same_reference(a::Reference, b::Reference)
    return @test test_is_same(a, b) === nothing
end

@testset "Reference" begin
    @test length(genomes(ref)) == 3
    @test n_seqs(ref) == 11

    buf = IOBuffer()
    BinBencherBackend.save(buf, ref)
    ref2 = Reference(IOBuffer(take!(buf)))

    test_is_same_reference(ref, ref2)
    ref3 = Reference(REF_PATH)
    test_is_same_reference(ref, ref3)

    # Test versions
    @test_throws "Found reference version 1" Reference(joinpath(DIR, "ref_v1.json"))

    # References with zero genomes error clearly, also with an empty rank,
    # and also when subsetting all genomes away
    seqs_json = """[["s1", 10, []]]"""
    for taxmaps in ["[]", "[[]]"]
        str = """{"version": 2, "taxmaps": $taxmaps, "genomes": [], "sequences": $seqs_json}"""
        @test_throws "zero genomes" Reference(IOBuffer(str))
    end
    @test_throws "zero genomes" subset(ref; genomes = Returns(false))

    # A failing subset! leaves the reference untouched and usable
    ref_copy = deepcopy(ref)
    @test_throws "zero genomes" subset!(ref_copy; genomes = Returns(false))
    @test_throws "no sequences" subset!(ref_copy; sequences = Returns(false))
    test_is_same_reference(ref, ref_copy)
    @test occursin("Genomes:    3", repr(MIME"text/plain"(), ref_copy))

    # Genomes, but no taxonomy: All genomes are placed directly under a top clade
    notax_str = """
    {
        "version": 2,
        "taxmaps": [],
        "genomes": [["A", 1, [["sA", 100]]], ["B", 1, [["sB", 100]]]],
        "sequences": [["s1", 50, [["sA", 1, 50]]], ["s2", 50, [["sB", 1, 50]]]]
    }
    """
    notax = Reference(IOBuffer(notax_str))
    @test BinBencherBackend.nranks(notax) == 2
    @test top_clade(notax).ngenomes == 2
    @test Set(top_clade(notax).children) == genomes(notax)
    notax_bins = Binning(IOBuffer("clustername\tcontigname\nb\ts1\n"), notax)
    @test n_recovered(notax_bins, 0.6, 0.6; level = 1) == 0
    buf = IOBuffer()
    BinBencherBackend.save(buf, notax)
    test_is_same_reference(notax, Reference(IOBuffer(take!(buf))))

    # More ranks than there are named ranks can be displayed
    deep_taxmaps = join(
        ["[[\"$(i == 1 ? "A" : "R$(i - 1)")\", \"R$(i)\"]]" for i in 1:10], ", "
    )
    deep_str = """
    {
        "version": 2,
        "taxmaps": [$deep_taxmaps],
        "genomes": [["A", 1, [["sA", 100]]]],
        "sequences": [["s1", 50, [["sA", 1, 50]]]]
    }
    """
    deep = Reference(IOBuffer(deep_str))
    @test repr(top_clade(deep)) == "<Rank 10> \"R10\", 1 genome"
    @test occursin("<Rank 9> \"R9\"", repr(MIME"text/plain"(), top_clade(deep)))
end

@testset "Sequence" begin
    s1 = Sequence("abc", 5)
    s2 = Sequence("abc abc", 6)
    s3 = Sequence(SubString(" abc", 2:4), 7)
    seqs = [s1, s2, s3]

    @test_throws Exception Sequence("abc", 0)
    @test_throws Exception Sequence("abc", -5)

    # Bad names
    @test_throws Exception Sequence("a\tb", 5)
    @test_throws Exception Sequence("a\rbc", 5)
    @test_throws Exception Sequence("\n\n", 5)

    # Ok names
    @test Sequence("", 1) isa Sequence
    @test Sequence("   ", 1) isa Sequence
    @test Sequence("\0\v\f", 1) isa Sequence
    @test Sequence("\xff\xff", 1) isa Sequence

    @test map(length, seqs) == [5, 6, 7]
    @test s1 == s3 # we might change this behaviour
    @test isequal(s1, s3)
    @test hash(s1) === hash(s3)
end

@testset "Genome" begin
    gens = sort!(collect(genomes(ref)); by = i -> i.name)
    @test is_organism(gens[1])
    @test is_organism(gens[2])
    @test is_virus(gens[3])
    @test !is_organism(gens[3])
    @test !is_virus(gens[1])

    for good_name in ["", "  ", "\xff\0\0"]
        @test Genome(good_name) isa Genome
    end
    for bad_name in ["\r\n", "abc\na", "a\ta"]
        @test_throws Exception Genome(bad_name)
    end
end

@testset "Clade" begin
    # Error message names the child and parents correctly
    g = Genome("myGenome")
    c = Clade("cladeOne", g)
    @test_throws (
        "Attempted to add parent \"cladeTwo\" to child \"myGenome\", " *
            "which already has parent \"cladeOne\""
    ) Clade("cladeTwo", g)
    Clade("cladeThree", c)
    @test_throws "to child \"cladeOne\", which already has parent \"cladeThree\"" Clade(
        "cladeFour", c
    )

    (gA, gB, gC) = sort!(collect(genomes(ref)); by = i -> i.name)
    D = gA.parent
    @test D === gB.parent
    E = gC.parent
    F = D.parent
    @test F === E.parent

    @testset "MRCA" begin
        @test mrca(gA, gB).name == "D"
        @test mrca(gA, gC).name == mrca(gB, gC).name == "F"

        @test mrca(gA, gB) === D
        @test mrca(gA, D) === D

        @test mrca(D, mrca(gA, gC)) === F
        @test mrca(F, F) === F
        @test mrca(gA, F) === F
    end

    @testset "Ancestors" begin
        @test collect(ancestors(gA)) == [D, F]
        @test collect(ancestors(gA)) == collect(ancestors(gB))
        @test only(ancestors(D)) === F
        @test isempty(collect(ancestors(F)))
        @test collect(ancestors(gC)) == [E, F]
        @test only(ancestors(E)) === F
        @test isempty(ancestors(F))
    end

    @testset "Descends from" begin
        lineages = Any[[gA, D, F], [gB, D, F], [gC, E, F]]
        for i in lineages
            for j in eachindex(i)
                for k in j:lastindex(i)
                    @test descends_from(i[j], i[k])
                    if k > j
                        @test !descends_from(i[k], i[j])
                    end
                end
            end
        end
    end
    @test !descends_from(gA, gB)
    @test !descends_from(gA, E)
    @test !descends_from(E, gB)
end

@testset "Subsetting" begin
    seq_pred = s -> length(s) ≥ 25
    genome_pred = !is_virus
    ref = Reference(IOBuffer(REF_STR))
    ref2 = subset(ref; sequences = seq_pred)
    ref3 = subset(ref; genomes = genome_pred)
    ref4 = subset(ref3; sequences = seq_pred)
    ref5 = subset(ref; sequences = seq_pred, genomes = genome_pred)

    refs = (ref, ref2, ref3, ref4, ref5)
    for i in 1:4, j in (i + 1):5
        @test refs[i] !== refs[j]
    end
    @test (ngenomes(ref3) + 1 == ngenomes(ref) == ngenomes(ref4) + 1 == ngenomes(ref5) + 1)
    @test (n_seqs(ref2) == n_seqs(ref4) == n_seqs(ref5))
    @test (n_seqs(ref) == n_seqs(ref3) != n_seqs(ref2))

    @test (
        top_clade(ref3).ngenomes + 1 ==
            top_clade(ref4).ngenomes + 1 ==
            top_clade(ref5).ngenomes + 1 ==
            top_clade(ref).ngenomes
    )

    # Test subsetting preserves relationship between seq names and their targets
    ref6 = deepcopy(ref)
    trgs_of = Dict(name => ref6.targets[i] for (name, i) in ref6.target_index_by_name)
    subset!(ref6; sequences = seq_pred)
    @test all(ref6.target_index_by_name) do (name, i)
        ref6.targets[i] === trgs_of[name]
    end

    # The sequences predicate is only evaluated once per sequence, so even a
    # non-deterministic predicate keeps sources and targets consistent
    ref7 = deepcopy(ref)
    n_calls = Ref(0)
    alternating = s -> isodd(n_calls[] += 1)
    subset!(ref7; sequences = alternating)
    @test n_calls[] == n_seqs(ref)
    mapped_in_targets = Set(seq for (seq, t) in ref7.targets if !isempty(t))
    in_sources = Set(
        first(i) for g in genomes(ref7) for src in g.sources for i in src.sequences
    )
    @test mapped_in_targets == in_sources
end

function test_is_same_binning(a::Binning, b::Binning)
    @test a.ref === b.ref
    @test [i.name for i in a.bins] == [i.name for i in b.bins]
    for field in [:recovered_asms, :recovered_genomes, :recalls, :precisions]
        @test getfield(a, field) == getfield(b, field)
    end
    return
end

@testset "Binning" begin
    ref = Reference(IOBuffer(REF_STR))
    bins = Binning(IOBuffer(CLUSTERS_STR), ref)

    @test bins isa Binning
    @test n_bins(bins) == 6

    @test n_recovered(bins, 0.6, 0.8) == 1
    @test n_recovered(bins, 0.6, 0.8; assembly = true) == 2
    @test n_recovered(bins, 0.6, 0.8; assembly = true, level = 2) == 1

    allgenomes = collect(genomes(ref))
    for (ir, recall) in enumerate(bins.recalls)
        for (ip, precision) in enumerate(bins.precisions)
            for (asm, mats) in
                [(true, bins.recovered_asms), (false, bins.recovered_genomes)]
                found = falses(length(allgenomes))
                for (ig, genome) in enumerate(allgenomes), bin in bins.bins
                    (nr, np) = recall_precision(genome, bin; assembly = asm)
                    found[ig] |= (nr >= recall && np >= precision)
                end
                @test sum(found) == mats[1][ip, ir]

                for (rank, mat) in zip(ref.clades, mats[2:end])
                    found = falses(length(rank))
                    for bin in bins.bins, (ic, clade) in enumerate(rank)
                        (nr, np) = recall_precision(clade, bin; assembly = asm)
                        found[ic] |= (nr >= recall && np >= precision)
                    end
                    @test sum(found) == mat[ip, ir]
                end
            end
        end
    end

    # Test filter_genomes works
    empty_binning = Binning(IOBuffer(CLUSTERS_STR), ref; filter_genomes = Returns(false))
    @test n_recovered(empty_binning, 0.6, 0.6) == 0
    @test all(m -> all(iszero, m), empty_binning.recovered_asms)
    @test all(m -> all(iszero, m), empty_binning.recovered_genomes)

    only_virus = Binning(IOBuffer(CLUSTERS_STR), ref; filter_genomes = is_virus)
    @test BinBencherBackend.n_nc(only_virus) == 0
    @test n_recovered(only_virus, 0.6, 0.6; assembly = true) == 1

    # This test depends on the exact state of the ref and binning used
    @test all(m -> all(isone, m), only_virus.recovered_asms)
    @test all(m -> all(iszero, m), only_virus.recovered_genomes)

    # Sequences mapping only to ignored genomes do not count as contamination
    small_ref = Reference(
        IOBuffer(
            """
            {
                "version": 2,
                "taxmaps": [[["A", "X"], ["B", "Y"]]],
                "genomes": [["A", 1, [["sA", 100]]], ["B", 1, [["sB", 100]]]],
                "sequences": [
                    ["s1", 100, [["sA", 1, 100]]],
                    ["s2", 100, [["sB", 1, 100]]],
                    ["s3", 100, []]
                ]
            }
            """
        )
    )
    A = only(g for g in genomes(small_ref) if g.name == "A")
    small_binning = Binning(
        IOBuffer("clustername\tcontigname\nb\ts1\nb\ts2\nb\ts3\n"),
        small_ref;
        filter_genomes = g -> g.name == "A",
    )
    small_bin = only(small_binning.bins)
    @test recall_precision(A, small_bin) == (recall = 1.0, precision = 1.0)
    @test recall_precision(A.parent, small_bin) == (recall = 1.0, precision = 1.0)
    @test n_recovered(small_binning, 0.99, 0.99) == 1

    # Text display does not contain leftover text from the reused buffer when the
    # matrix is shorter than the display of the reference
    tiny = Binning(IOBuffer(CLUSTERS_STR), ref; recalls = [0.5], precisions = [0.5])
    shown = repr(MIME"text/plain"(), tiny)
    @test count("Genomes:", shown) == 1
    @test endswith(rstrip(shown), r"0\.5\s+\d+")

    bins2 = Binning(CLUSTERS_PATH, ref)
    test_is_same_binning(bins, bins2)

    @test bins.bin_genome_stats.mean_bin_recall ≈ 0.4916363636363636
    @test bins.bin_genome_stats.mean_bin_precision ≈ 1
    @test bins.bin_asm_stats.mean_bin_recall ≈ 0.636734693877551
    @test bins.bin_asm_stats.mean_bin_precision ≈ 1

    # The same sequence may not be present in the same bin twice, even if not disjoint
    dup_str = "clustername\tcontigname\nb\ts1\nb\ts2\nb\ts1\n"
    @test_throws "\"s1\" is present multiple times in bin \"b\"" Binning(
        IOBuffer(dup_str), ref; disjoint = false
    )
    @test_throws "multiple times" Bin("b", ref, [Sequence("s1", 25), Sequence("s1", 25)])
    # But the same sequence may be present in different bins
    nondisjoint = Binning(
        IOBuffer("clustername\tcontigname\na\ts1\nb\ts1\n"), ref; disjoint = false
    )
    @test !is_disjoint(nondisjoint)

    # Test is disjoint
    @test is_disjoint(bins)
end

@testset "Threshold lookup" begin
    ref = Reference(IOBuffer(REF_STR))
    bins = Binning(IOBuffer(CLUSTERS_STR), ref)
    rs = bins.recalls
    ps = bins.precisions
    cases = [
        (n_recovered, n_recovered_bounds, false, bins.recovered_genomes),
        (n_recovered, n_recovered_bounds, true, bins.recovered_asms),
        (n_passing_bins, n_passing_bins_bounds, false, bins.bin_genomes),
        (n_passing_bins, n_passing_bins_bounds, true, bins.bin_asms),
    ]
    for (exact, bounds, assembly, mats) in cases, (i, m) in enumerate(mats)
        level = i - 1
        for (ri, r) in enumerate(rs), (pi, p) in enumerate(ps)
            @test exact(bins, r, p; level, assembly) == m[pi, ri]
            @test bounds(bins, r, p; level, assembly) == (lower = m[pi, ri], upper = m[pi, ri])
        end
        # Between thresholds on both axes
        for ri in 1:(length(rs) - 1), pi in 1:(length(ps) - 1)
            r = (rs[ri] + rs[ri + 1]) / 2
            p = (ps[pi] + ps[pi + 1]) / 2
            @test bounds(bins, r, p; level, assembly) ==
                (lower = m[pi + 1, ri + 1], upper = m[pi, ri])
            # Exact on recall, between on precision
            @test bounds(bins, rs[ri], p; level, assembly) ==
                (lower = m[pi + 1, ri], upper = m[pi, ri])
        end
    end

    # Default thresholds: 0.6, 0.7, 0.8, 0.9, 0.95, 0.99
    @test n_recovered(bins, 0.6, 0.8) == bins.recovered_genomes[1][3, 1]
    # Within tolerance
    @test n_recovered(bins, 0.6 + 1.0e-7, 0.8 - 1.0e-7) == n_recovered(bins, 0.6, 0.8)
    @test n_recovered(bins, 0.6f0, 0.8f0) == n_recovered(bins, 0.6, 0.8)
    @test n_recovered_bounds(bins, 0.6f0, 0.8f0) == n_recovered_bounds(bins, 0.6, 0.8)
    @test n_recovered(bins, 0.6 + 1.0e-4, 0.8; atol = 1.0e-3) == n_recovered(bins, 0.6, 0.8)
    @test_throws ArgumentError n_recovered(bins, 0.6 + 1.0e-4, 0.8)

    for f in [n_recovered, n_passing_bins]
        # Not a benchmarked threshold
        @test_throws ArgumentError f(bins, 0.65, 0.8)
        @test_throws ArgumentError f(bins, 0.6, 0.85)
        @test_throws ArgumentError f(bins, 0.4, 0.8)
    end
    for f in [n_recovered, n_passing_bins, n_recovered_bounds, n_passing_bins_bounds]
        # Outside [0, 1]
        @test_throws ArgumentError f(bins, -0.1, 0.8)
        @test_throws ArgumentError f(bins, 0.6, 1.1)
        @test_throws ArgumentError f(bins, NaN, 0.8)
        # Bad level
        @test_throws ArgumentError f(bins, 0.6, 0.8; level = -1)
        @test_throws ArgumentError f(bins, 0.6, 0.8; level = length(bins.recovered_asms))
    end
    for f in [n_recovered_bounds, n_passing_bins_bounds]
        # Outside the range of benchmarked thresholds
        @test_throws ArgumentError f(bins, 0.5, 0.8)
        @test_throws ArgumentError f(bins, 0.6, 0.5)
        @test_throws ArgumentError f(bins, 0.995, 0.8)
        @test_throws ArgumentError f(bins, 0.6, 1.0)
    end

    # n_nc is `nothing` when 0.9/0.95 were not benchmarked
    @test BinBencherBackend.n_nc(bins) isa Int
    tiny = Binning(IOBuffer(CLUSTERS_STR), ref; recalls = [0.5], precisions = [0.5])
    @test isnothing(BinBencherBackend.n_nc(tiny))
end

@testset "Gold standard" begin
    ref = Reference(IOBuffer(REF_STR))
    gold_standards = [
        gold_standard(ref; disjoint = true)
        gold_standard(ref; disjoint = false)
    ]
    for (disjoint, bins) in zip([true, false], gold_standards)
        @test bins isa Binning
        @test n_bins(bins) == ngenomes(ref)
        @test is_disjoint(bins) == disjoint
    end

    # Disjoint gold standard assigns each sequence to the genome with the longest span,
    # also when the longest span wraps around a circular source
    refstr = """
    {
        "version": 2,
        "taxmaps": [[["A", "X"], ["B", "X"], ["C", "X"]]],
        "genomes": [
            ["A", 1, [["sA", 100]]],
            ["B", 1, [["sB", 100]]],
            ["C", 1, [["sC", 100]]]
        ],
        "sequences": [
            ["s1", 50, [["sA", 1, 5], ["sB", 1, 50]]],
            ["s2", 50, [["sA", 1, 40], ["sC", 71, 20]]]
        ]
    }
    """
    small_ref = Reference(IOBuffer(refstr))
    gs = gold_standard(small_ref; disjoint = true)
    bin_of = Dict(only(s.name for s in b.sequences) => b.name for b in gs.bins)
    @test bin_of == Dict("s1" => "bin_B", "s2" => "bin_C")

    # Using a Binning
    gs = gold_standard(ref, binning)
    @test gs isa Binning
    @test is_disjoint(gs) # the default
    @test n_bins(gs) == ngenomes(ref)
end

n_sequences(x::Binning) = sum(i -> length(i.sequences), x.bins; init = 0)

@testset "Binsplitting" begin
    # Make the unsplit output by merging genomes A and B
    unsplit_clusters_str = let
        lines = IOBuffer(CLUSTERS_STR)
        output = IOBuffer()
        println(output, readline(lines)) # header
        for line in eachline(lines)
            (bin, contig) = split(line, '\t')
            contigno = parse(Int, contig[2:end]) # called e.g. s9
            sample = if contigno in 1:3
                1
            else
                @assert contigno in 4:11
                2
            end
            cname = contigno < 9 ? 1 : 2
            println(output, "$(cname)\tS$(sample)C$(contig)")
        end
        seekstart(output)
        String(read(output))
    end
    unsplit_ref_str = let
        s = REF_STR
        for i in 1:11
            S = i ≤ 3 ? 1 : 2
            s = replace(s, "\"s$(i)\"" => "\"S$(S)Cs$(i)\"")
        end
        s
    end

    # Test basic binsplitting
    ref = Reference(IOBuffer(unsplit_ref_str))
    unsplit = Binning(IOBuffer(unsplit_clusters_str), ref)
    split_bins = Binning(IOBuffer(unsplit_clusters_str), ref; binsplit_separator = 'C')

    @test n_bins(unsplit) < n_bins(split_bins)
    @test n_sequences(unsplit) == n_sequences(split_bins)

    # Test binsplitting with gold_standard
    gs = gold_standard(ref; binsplit_separator = 'C')
    @test n_bins(gs) >= length(genomes(ref))
    target_seqs = count(ref.targets) do (s, t)
        !isempty(t)
    end

    @test n_sequences(gs) == target_seqs
    @test is_disjoint(gs)

    # Test manual bins
    @test all(gs.bins) do bin
        s = Set([split(s.name, 'C')[1] for s in bin.sequences])
        length(s) == 1
    end
end

@testset "From gzipped" begin
    mktempdir() do path
        ref_path = joinpath(path, "ref.json.gz")
        open(io -> write(io, transcode(GzipCompressor, REF_STR)), ref_path, "w")
        ref1 = Reference(REF_PATH)
        ref2 = Reference(ref_path)
        test_is_same_reference(ref1, ref2)

        bins_path = joinpath(path, "bins.tsv.gz")
        open(io -> write(io, transcode(GzipCompressor, CLUSTERS_STR)), bins_path, "w")
        bins1 = Binning(CLUSTERS_PATH, ref1)
        bins2 = Binning(bins_path, ref1)
        test_is_same_binning(bins1, bins2)
    end
end

@testset "ref.json specifics" begin
    ref = Reference(IOBuffer(REF_STR))
    sources_by_name = Dict{String, Source{Genome}}()
    for genome in genomes(ref), source in genome.sources
        sources_by_name[source.name] = source
    end
    # Circular mapping
    C2 = sources_by_name["subjC2"]
    @test C2.length == 100
    @test C2.assembly_size == 21 + 20 + 21

    @test sources_by_name["subjC3"].assembly_size == 0

    # Ref conforms to schema
    schema = JSONSchema.Schema(String(open(read, joinpath(DIR, "schema.json"))))
    ref_data = Dict(JSON.parse(REF_STR))
    @test isnothing(JSONSchema.validate(schema, ref_data))
end

@testset "Adjusted rand index" begin
    # We need to make a Binning and a ref having the same sequences, so we can
    # compare the gold standard of the ref with the binning
    ref = Reference(IOBuffer(REF_STR))
    bins = Binning(IOBuffer(CLUSTERS_STR), ref)
    seqs_in_bin = foldl((i, j) -> union!(i, j.sequences), bins.bins; init = Set{Sequence}())
    setdiff!(seqs_in_bin, [seq for (seq, t) in ref.targets if isempty(t)])
    seq_names = Set([s.name for s in seqs_in_bin])

    @test adjusted_rand_index(bins, bins) === 1.0
    gs = gold_standard(ref)
    @test_throws ArgumentError adjusted_rand_index(bins, gs)
    ari = adjusted_rand_index(bins, gs; intersection = true)
    @test 0.0 ≤ ari ≤ 1.0

    ref2 = subset(ref; sequences = in(seqs_in_bin))
    lines = filter(collect(eachline(IOBuffer(CLUSTERS_STR)))) do line
        name = split(line, '\t')[2]
        name == "contigname" || in(name, seq_names)
    end

    bins = Binning(IOBuffer(join(lines, '\n')), ref2)
    ari = adjusted_rand_index(bins, gold_standard(ref2))
    @test 0.0 ≤ ari ≤ 1.0

    # Compare to a known ARI
    refstr = """
    {
        "version": 2,
        "taxmaps": [[["A", "X"], ["B", "X"], ["C", "X"], ["D", "X"]]],
        "genomes": [
            ["A", 1, [["A", 5]]],
            ["B", 1, [["B", 5]]],
            ["C", 1, [["C", 5]]],
            ["D", 1, [["D", 5]]]
        ],
        "sequences": [
            ["s0", 1, [["A",1,1]]],
            ["s1", 2, [["A",2,3]]],
            ["s2", 1, [["B",1,1]]],
            ["s3", 1, [["B",2,2]]],
            ["s4", 3, [["C",1,3]]],
            ["s5", 1, [["C",4,4]]],
            ["s6", 1, [["D",1,1]]],
            ["s7", 2, [["D",2,3]]],
            ["s8", 1, [["D",4,4]]],
            ["s9", 1, [["D",5,5]]]
        ]
    }
    """
    cstr = """clustername\tcontigname
    a\ts0
    b\ts1
    c\ts2
    c\ts3
    d\ts4
    e\ts5
    f\ts6
    g\ts7
    g\ts8
    g\ts9
    """

    ref = Reference(IOBuffer(refstr))
    binning = Binning(IOBuffer(cstr), ref)

    # Equivalent to this clustering
    # [1,1,1,2,2,3,3,3,3,4,4,4,4,4]
    # [1,2,2,3,3,4,4,4,5,6,7,7,7,7]
    # Hard coded value obtained using Clustering.jl's randindex

    @test isapprox(adjusted_rand_index(binning, gold_standard(ref)), 0.6560268794624108)

    # Two identical binnings each with a single bin
    single = Binning(IOBuffer("clustername\tcontigname\na\ts0\na\ts1\na\ts2\n"), ref)
    @test adjusted_rand_index(single, single) === 1.0

    # Very long sequences, such that the product of the pair counts overflows UInt128.
    # Binning a has both sequences in one bin, b has them separately. Then the
    # number of co-binned pairs equals the expected number, so the ARI is 0.
    big_refstr = """
    {
        "version": 2,
        "taxmaps": [[["A", "X"]]],
        "genomes": [["A", 1, [["A", 5]]]],
        "sequences": [["s1", 10000000000, []], ["s2", 10000000000, []]]
    }
    """
    big_ref = Reference(IOBuffer(big_refstr))
    big_a = Binning(IOBuffer("clustername\tcontigname\na\ts1\na\ts2\n"), big_ref)
    big_b = Binning(IOBuffer("clustername\tcontigname\na\ts1\nb\ts2\n"), big_ref)
    @test isapprox(adjusted_rand_index(big_a, big_b), 0.0; atol = 1.0e-9)
end

# Returns the exception thrown by `f`, or `nothing`
function caught(f)
    return try
        f()
        nothing
    catch e
        e
    end
end

# Check that `f` throws an exception of type `T`, whose message contains `msg`
function throws(f, ::Type{T}, msg::AbstractString) where {T}
    e = caught(f)
    return e isa T && occursin(msg, sprint(showerror, e))
end

@testset "Exceptions" begin
    @test ReferenceError <: BinBencherException
    @test BinningError <: BinBencherException
    @test sprint(showerror, ReferenceError("abc")) == "ReferenceError: abc"
    @test sprint(showerror, BinningError("abc", 5)) == "BinningError: line 5: abc"
    @test sprint(showerror, BinningError("abc", nothing)) == "BinningError: abc"

    @testset "Reference files" begin
        # A small valid reference, with each part overridable
        function ref_json(;
                version = "2",
                genomes = """[["A", 1, [["sA", 100]]], ["B", 1, [["sB", 100]]]]""",
                sequences = """[["s1", 50, [["sA", 1, 50]]], ["s2", 50, [["sB", 1, 50]]]]""",
                taxmaps = """[[["A", "X"], ["B", "X"]]]""",
            )
            return """{"version": $version, "genomes": $genomes, "sequences": $sequences, "taxmaps": $taxmaps}"""
        end
        load(; kwargs...) = Reference(IOBuffer(ref_json(; kwargs...)))
        @test load() isa Reference

        # Malformed JSON
        for str in ["not json", "[1, 2]", "{\"version\": 2}", "{\"genomes\": []}"]
            @test throws(() -> Reference(IOBuffer(str)), ReferenceError, "Malformed JSON")
        end
        @test throws(() -> load(; version = "\"2\""), ReferenceError, "integer field \"version\"")
        @test throws(() -> load(; genomes = "5"), ReferenceError, "Malformed JSON")
        @test throws(() -> load(; genomes = """[["A", "x", []]]"""), ReferenceError, "Malformed JSON")
        @test throws(() -> load(; version = "1"), ReferenceError, "Found reference version 1")

        # Empty
        @test throws(() -> load(; genomes = "[]", taxmaps = "[]"), ReferenceError, "zero genomes")
        @test throws(() -> load(; sequences = "[]"), ReferenceError, "zero sequences")

        # Genomes and sources
        @test throws(ReferenceError, "Invalid identifier") do
            load(; genomes = """[["A\\r", 1, [["sA", 100]]], ["B", 1, [["sB", 100]]]]""")
        end
        for flags in [-1, 8]
            @test throws(ReferenceError, "invalid flags $flags") do
                load(; genomes = """[["A", $flags, [["sA", 100]]], ["B", 1, [["sB", 100]]]]""")
            end
        end
        @test throws(ReferenceError, "Genome \"A\" is present more than once") do
            load(; genomes = """[["A", 1, [["sA", 100]]], ["A", 1, [["sB", 100]]]]""")
        end
        @test throws(ReferenceError, "Invalid identifier") do
            load(; genomes = """[["A", 1, [["s\\nA", 100]]], ["B", 1, [["sB", 100]]]]""")
        end
        @test throws(ReferenceError, "has length 0") do
            load(; genomes = """[["A", 1, [["sA", 0]]], ["B", 1, [["sB", 100]]]]""")
        end
        @test throws(ReferenceError, "has source \"sA\" more than once") do
            load(; genomes = """[["A", 1, [["sA", 100], ["sA", 10]]], ["B", 1, [["sB", 100]]]]""")
        end
        @test throws(ReferenceError, "Duplicate source") do
            load(; genomes = """[["A", 1, [["sA", 100]]], ["B", 1, [["sA", 100]]]]""")
        end

        # Sequences
        @test throws(ReferenceError, "Invalid identifier") do
            load(; sequences = """[["s\\t1", 50, []]]""")
        end
        @test throws(ReferenceError, "has length 0") do
            load(; sequences = """[["s1", 0, []]]""")
        end
        @test throws(ReferenceError, "no such source") do
            load(; sequences = """[["s1", 50, [["sQ", 1, 50]]]]""")
        end
        @test throws(ReferenceError, "valid source indices are 1:100") do
            load(; sequences = """[["s1", 50, [["sA", 60, 101]]]]""")
        end
        @test throws(ReferenceError, "Sequence \"s1\" is present more than once") do
            load(; sequences = """[["s1", 50, []], ["s1", 50, []]]""")
        end

        # Taxonomy (parsed in another task)
        @test throws(ReferenceError, "child \"A\" is present more than once") do
            load(; taxmaps = """[[["A", "X"], ["A", "Y"], ["B", "X"]]]""")
        end
        @test throws(ReferenceError, "does not exist on the previous rank") do
            load(; taxmaps = """[[["A", "X"], ["B", "X"], ["Q", "X"]]]""")
        end
        @test throws(ReferenceError, "child \"B\" has no parent") do
            load(; taxmaps = """[[["A", "X"]]]""")
        end
        @test throws(ReferenceError, "Invalid identifier") do
            load(; taxmaps = """[[["A", "X\\n"], ["B", "X\\n"]]]""")
        end
    end

    @testset "Cluster files" begin
        ref = Reference(IOBuffer(REF_STR))
        header = "clustername\tcontigname\n"
        load(str; kwargs...) = Binning(IOBuffer(str), ref; kwargs...)
        # Returns (message, line) of the thrown BinningError, or nothing
        function binning_error(str; kwargs...)
            e = caught(() -> load(str; kwargs...))
            return e isa BinningError ? (e.msg, e.line) : e
        end

        @test load(header * "a\ts1\n\nb\ts2\r\n") isa Binning
        @test binning_error("") == ("Cluster file is empty, expected a header", 1)
        (msg, line) = binning_error("name\tcontig\na\ts1\n")
        @test occursin("Expected header", msg) && line == 1
        @test binning_error(header * "a\ts1\na s2\n") ==
            ("Expected two tab-separated fields, found one", 3)
        @test binning_error(header * "a\ts1\n\na\ts2\tx\n") ==
            ("Expected two tab-separated fields, found more", 4)
        @test binning_error(header * "a\ts1\na\tsQ\n") ==
            ("Sequence \"sQ\" is not in the reference", 3)
        (msg, line) = binning_error(header * "a\ts1\na\rb\ts2\n")
        @test occursin("Invalid identifier", msg) && line == 3
        (msg, line) = binning_error(header * "a\ts1\n"; binsplit_separator = "C")
        @test occursin("Binsplit separator", msg) && line == 2
        @test binning_error(header * "a\ts1\nb\ts2\na\ts1\n"; disjoint = false) ==
            ("Sequence \"s1\" is present multiple times in bin \"a\"", nothing)
        @test binning_error(header * "a\ts1\nb\ts2\nb\ts1\n") == (
            "Sequence \"s1\" is present in multiple bins, but `disjoint` is set to `true`", 4,
        )
    end

    @testset "API misuse" begin
        ref = Reference(IOBuffer(REF_STR))
        g = Genome("g")
        Clade("c", g)
        for f in [
                () -> Sequence("a\tb", 5),
                () -> Sequence("a", 0),
                () -> Genome("a\nb"),
                () -> Source(g, "s", 0),
                () -> Source(g, "s\r", 5),
                () -> Clade("d", g),
                () -> Bin("b", ref, Sequence[]),
                () -> Bin("b", ref, [Sequence("sQ", 5)]),
                () -> Bin("b", ref, [Sequence("s1", 25), Sequence("s1", 25)]),
                () -> Binning(IOBuffer(CLUSTERS_STR), ref; recalls = [0.0]),
                () -> Binning(IOBuffer(CLUSTERS_STR), ref; precisions = [0.5, 0.5]),
                () -> Binning(IOBuffer(CLUSTERS_STR), ref; recalls = Float64[]),
                () -> subset(ref; genomes = Returns(false)),
                () -> subset(ref; sequences = Returns(false)),
                () -> gold_standard(ref; binsplit_separator = "QQ"),
                () -> gold_standard(ref, [Sequence("sQ", 5)]),
            ]
            @test caught(f) isa ArgumentError
        end
        bin = Bin("b", ref, [Sequence("s1", 25)])
        @test caught(() -> Binning([bin, bin], ref)) isa ArgumentError
    end
end
