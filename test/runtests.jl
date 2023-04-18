using Test
using VambBenchmarks

const DIR = joinpath(dirname(dirname(pathof(VambBenchmarks))), "files")
REFSTR = read(joinpath(DIR, "ref.json"), String)

@assert isdir(DIR)
ngenomes(ref) = length(genomes(ref))

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
