# VambBenchmarks
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://jakobnissen.github.io/VambBenchmarks.jl/dev)

This package is used to benchmark and interactively explore the results of metagenomic binning given a dataset.

## Installation
* Install Julia (ideally using Juliaup: https://github.com/JuliaLang/juliaup) or from the official website `www.julialang.org`
* Install VambBenchmarks.jl: `julia --startup-file=no --project=@vambbench -e 'using Pkg; Pkg.add(url="https://github.com/jakobnissen/VambBenchmarks.jl")'`

## Documentation
Basic usage:
* Launch Julia: `julia --startup-file=no --project=@vambbench`
```julia
using VambBenchmarks
ref = Reference("files/ref.json")
bins = Binning("files/clusters.tsv", ref)
print_matrix(bins)
```

For more, see See the [documentation](https://jakobnissen.github.io/VambBenchmarks.jl/dev)
