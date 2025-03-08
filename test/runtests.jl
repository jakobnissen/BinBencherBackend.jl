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
    n_seqs,
    n_bins,
    subset,
    subset!,
    is_disjoint,
    adjusted_rand_index

using CodecZlib: GzipCompressor
using JSONSchema: JSONSchema
using JSON: JSON

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

    @test n_recovered(bins, 0.4, 0.71) == 1
    @test n_recovered(bins, 0.4, 0.71; assembly = true) == 2
    @test n_recovered(bins, 0.4, 0.71; assembly = true, level = 2) == 1

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
    @test n_recovered(empty_binning, 0.1, 0.1) == 0
    @test all(m -> all(iszero, m), empty_binning.recovered_asms)
    @test all(m -> all(iszero, m), empty_binning.recovered_genomes)

    only_virus = Binning(IOBuffer(CLUSTERS_STR), ref; filter_genomes = is_virus)
    @test BinBencherBackend.n_nc(only_virus) == 0
    @test n_recovered(only_virus, 0.1, 0.1; assembly = true) == 1

    # This test depends on the exact state of the ref and binning used
    @test all(m -> all(isone, m), only_virus.recovered_asms)
    @test all(m -> all(iszero, m), only_virus.recovered_genomes)

    bins2 = Binning(CLUSTERS_PATH, ref)
    test_is_same_binning(bins, bins2)

    @test bins.bin_genome_stats.mean_bin_recall ≈ 0.4916363636363636
    @test bins.bin_genome_stats.mean_bin_precision ≈ 1
    @test bins.bin_asm_stats.mean_bin_recall ≈ 0.636734693877551
    @test bins.bin_asm_stats.mean_bin_precision ≈ 1

    # Test is disjoint
    @test is_disjoint(bins)
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
    ref_data = copy(JSON.parse(REF_STR))::Dict
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

    isapprox(adjusted_rand_index(binning, gold_standard(ref)), 0.6560268794624108)
end
