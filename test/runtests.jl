using Test
using VambBenchmarks

const DIR = joinpath(dirname(dirname(pathof(VambBenchmarks))), "files")
REFSTR = read(joinpath(DIR, "ref.json"), String)

@assert isdir(DIR)
ngenomes(ref) = length(genomes(ref))

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
    global ref = Reference(IOBuffer(REFSTR))
    global binning = open(i -> Binning(i, ref), joinpath(DIR, "clusters.tsv"))
    global bins = sort!(collect(binning.bins); by=i -> i.name)

    @test ref isa Reference
    @test binning isa Binning
    @test bins isa Vector{Bin}
end

@testset "Reference" begin
    @test length(genomes(ref)) == 3
    @test nseqs(ref) == 11

    buf = IOBuffer()
    VambBenchmarks.save(buf, ref)
    ref2 = Reference(IOBuffer(take!(buf)))

    @test genomes(ref) == genomes(ref2)
    @test ref.targets_by_name == ref2.targets_by_name
    cladenames(ref) = [[c.name for c in v] for v in ref.clades]
    @test cladenames(ref) == cladenames(ref2)
end

@testset "Sequence" begin
    s1 = Sequence("abc", 5)
    s2 = Sequence("abc abc", 6)
    s3 = Sequence(SubString(" abc", 2:4), 7)
    seqs = [s1, s2, s3]

    @test_throws Exception Sequence("abc", 0)
    @test_throws Exception Sequence("abc", -5)

    # Bad names
    @test_throws Exception Sequence("", 5)
    @test_throws Exception Sequence(" abc", 5)
    @test_throws Exception Sequence("abc ", 5)

    @test map(length, seqs) == [5, 6, 7]
    @test s1 == s3 # we might change this behaviour
    @test isequal(s1, s3)
    @test hash(s1) === hash(s3)
end

@testset "Genome" begin
    gens = sort!(collect(genomes(ref)); by=i -> i.name)
    @test is_organism(gens[1])
    @test is_organism(gens[2])
    @test is_virus(gens[3])
    @test !is_organism(gens[3])
    @test !is_virus(gens[1])
end

@testset "Clade" begin
    (gA, gB, gC) = sort!(collect(genomes(ref)); by=i -> i.name)

    @test mrca(gA, gB).name == "D"
    @test mrca(gA, gC).name == mrca(gB, gC).name == "F"

    D = mrca(gA, gB)
    @test mrca(gA, D) === D

    F = mrca(D, mrca(gA, gC))
    @test mrca(F, F) == F
    @test mrca(gA, F) == F
end

@testset "Subsetting" begin
    seq_pred = s -> length(s) ≥ 25
    genome_pred = !is_virus
    ref = Reference(IOBuffer(REFSTR))
    ref2 = subset(ref; sequences=seq_pred)
    ref3 = subset(ref; genomes=genome_pred)
    ref4 = subset(ref3; sequences=seq_pred)
    ref5 = subset(ref; sequences=seq_pred, genomes=genome_pred)

    refs = (ref, ref2, ref3, ref4, ref5)
    for i in 1:4, j in (i + 1):5
        @test refs[i] !== refs[j]
    end
    @test (ngenomes(ref3) + 1 == ngenomes(ref) == ngenomes(ref4) + 1 == ngenomes(ref5) + 1)
    @test (nseqs(ref2) == nseqs(ref4) == nseqs(ref5))
    @test (nseqs(ref) == nseqs(ref3) != nseqs(ref2))

    @test (
        top_clade(ref3).ngenomes + 1 ==
        top_clade(ref4).ngenomes + 1 ==
        top_clade(ref5).ngenomes + 1 ==
        top_clade(ref).ngenomes
    )
end

@testset "Binning" begin
    ref = Reference(IOBuffer(REFSTR))
    bins = open(i -> Binning(i, ref), joinpath(DIR, "clusters.tsv"))

    @test bins isa Binning
    @test nbins(bins) == 6

    allgenomes = collect(genomes(ref))
    for (ir, recall) in enumerate(bins.recalls)
        for (ip, precision) in enumerate(bins.precisions)
            for (asm, mats) in
                [(true, bins.recovered_asms), (false, bins.recovered_genomes)]
                found = falses(length(allgenomes))
                for (ig, genome) in enumerate(allgenomes), bin in bins.bins
                    (nr, np) = recall_precision(genome, bin; assembly=asm)
                    found[ig] |= (nr >= recall && np >= precision)
                end
                @test sum(found) == mats[1][ip, ir]

                for (rank, mat) in zip(ref.clades, mats[2:end])
                    found = falses(length(rank))
                    for bin in bins.bins, (ic, clade) in enumerate(rank)
                        (nr, np) = recall_precision(clade, bin; assembly=asm)
                        found[ic] |= (nr >= recall && np >= precision)
                    end
                    @test sum(found) == mat[ip, ir]
                end
            end
        end
    end
end

@testset "Gold standard" begin
    ref = Reference(IOBuffer(REFSTR))
    gold_standards = [
        gold_standard(ref; disjoint=true)
        gold_standard(ref; disjoint=false)
    ]
    for bins in gold_standards
        @test bins isa Binning
        @test nbins(bins) == ngenomes(ref)
    end
    non_disjoint = last(gold_standards)
    @test non_disjoint.recoverable_genomes == non_disjoint.recovered_genomes[1][1, :]
end
