using Test
using VambBenchmarks

const DIR = joinpath(dirname(dirname(pathof(VambBenchmarks))), "files")
REFSTR = read(joinpath(DIR, "ref.json"), String)

@assert isdir(DIR)

@testset "Reference" begin
    ref = Reference(IOBuffer(REFSTR))
    
    @test ref isa Reference
    @test ngenomes(ref) == 3
    @test nseqs(ref) == 11

    buf = IOBuffer()
    VambBenchmarks.save(buf, ref)
    ref2 = Reference(IOBuffer(take!(buf)))

    @test ref.genomes == ref2.genomes
    @test ref.targets_by_name == ref2.targets_by_name
    cladenames(ref) = [[c.name for c in v] for v in ref.clades]
    @test cladenames(ref) == cladenames(ref2)
end

@testset "Filtering sequences" begin
    ref = Reference(IOBuffer(REFSTR))
    ref2 = filter_sequences(s -> s.length ≥ 25, ref)

    @test ref2 isa Reference
    @test ref2 !== ref
    @test all(i -> first(i).length >= 25, values(ref2.targets_by_name))
end

@testset "Filtering genomes" begin
    ref = Reference(IOBuffer(REFSTR))
    ref2 = filter_genomes(g -> g.name != "gC", ref)

    @test ref2 isa Reference
    @test ref2 !== ref
    @test ngenomes(ref2) + 1 == ngenomes(ref)
    @test top_clade(ref2).ngenomes + 1 == top_clade(ref).ngenomes
end

@testset "Binning" begin
    ref = Reference(IOBuffer(REFSTR))
    bins = open(i -> Binning(i, ref), joinpath(DIR, "clusters.tsv"))

    @test bins isa Binning
    @test nbins(bins) == 6
end

@testset "Gold standard" begin
    ref = Reference(IOBuffer(REFSTR))
    bins = gold_standard(ref)

    @test bins isa Binning
    @test nbins(bins) == ngenomes(ref)
    @test bins.recoverable_genomes == bins.recovered_genomes[1][1, :]
end
