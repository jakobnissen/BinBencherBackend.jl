using Documenter
using VambBenchmarks

meta = quote
    using VambBenchmarks

    (ref, binning, genome, bin) = let
        dir = joinpath(Base.pkgdir(VambBenchmarks), "files")
        ref = open(i -> Reference(i), joinpath(dir, "ref.json"))
        genome = first(sort!(collect(genomes(ref)); by=i -> i.name))
        bins = open(i -> Binning(i, ref), joinpath(dir, "clusters.tsv"))
        bin = first(bins.bins)
        (ref, bins, genome, bin)
    end
end

DocMeta.setdocmeta!(VambBenchmarks, :DocTestSetup, meta; recursive=true)

makedocs(
    sitename = "VambBenchmarks.jl",
    modules = [VambBenchmarks],
    pages = [
        "Home" => "index.md",
        "Walkthrough" => "walkthrough.md",
        "Reference" => "reference.md"
        ],
    checkdocs = :all
)

deploydocs(
    repo = "github.com/jakobnissen/VambBenchmarks.jl.git",
    push_preview = true
)
