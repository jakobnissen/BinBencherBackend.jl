# VambBenchmarks
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://jakobnissen.github.io/VambBenchmarks.jl/dev)
[![Latest Release](https://img.shields.io/github/release/jakobnissen/VambBenchmarks.jl.svg)](https://github.com/jakobnissen/VambBenchmarks.jl/releases/latest)

This package is used to benchmark and interactively explore the results of metagenomic binning given a dataset.

## Installation
* Install Julia (ideally using Juliaup: https://github.com/JuliaLang/juliaup) or from the official website `www.julialang.org`
* Use the Julia package manager to install this package: `] add VambBenchmarks`

## Documentation
Basic usage:
```julia
using VambBenchmarks
ref = Reference("files/ref.json")
bins = Binning("files/clusters.tsv", ref)
print_matrix(bins)
```

For more, see the [documentation](https://jakobnissen.github.io/VambBenchmarks.jl/dev)
