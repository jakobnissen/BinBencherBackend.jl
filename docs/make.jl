using Documenter
using VambBenchmarks

meta = quote
    using VambBenchmarks

    (path_to_ref_file, ref, binning, genome, bin) = let
        dir = joinpath(Base.pkgdir(VambBenchmarks), "files")
        path_to_ref_file = joinpath(dir, "ref.json")
        path_to_bins_file = joinpath(dir, "clusters.tsv")
        ref = open(i -> Reference(i), path_to_ref_file)
        genome = first(sort!(collect(genomes(ref)); by=i -> i.name))
        bins = open(i -> Binning(i, ref), path_to_bins_file)
        bin = first(bins.bins)
        (path_to_ref_file, ref, bins, genome, bin)
    end
end

DocMeta.setdocmeta!(VambBenchmarks, :DocTestSetup, meta; recursive=true)

makedocs(;
    sitename="VambBenchmarks.jl",
    modules=[VambBenchmarks],
    pages=[
        "Home" => "index.md",
        "Walkthrough" => "walkthrough.md",
        "Reference" => "reference.md",
    ],
    checkdocs=:all,
)

deploydocs(; repo="github.com/jakobnissen/VambBenchmarks.jl.git", push_preview=true)
