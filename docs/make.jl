using Documenter
using BinBencherBackend

meta = quote
    using BinBencherBackend

    (path_to_ref_file, path_to_bins_file, ref, binning, genome, bin) = let
        dir = joinpath(Base.pkgdir(BinBencherBackend), "files")
        path_to_ref_file = joinpath(dir, "ref.json")
        path_to_bins_file = joinpath(dir, "clusters.tsv")
        ref = open(i -> Reference(i), path_to_ref_file)
        genome = first(sort!(collect(genomes(ref)); by=i -> i.name))
        bins = open(i -> Binning(i, ref), path_to_bins_file)
        bin = first(bins.bins)
        (path_to_ref_file, path_to_bins_file, ref, bins, genome, bin)
    end
end

DocMeta.setdocmeta!(BinBencherBackend, :DocTestSetup, meta; recursive=true)

makedocs(;
    sitename="BinBencherBackend.jl",
    modules=[BinBencherBackend],
    pages=[
        "Home" => "index.md",
        "Walkthrough" => "walkthrough.md",
        "Reference" => "reference.md",
    ],
    checkdocs=:all,
)

deploydocs(; repo="github.com/jakobnissen/BinBencherBackend.jl.git", push_preview=true)
