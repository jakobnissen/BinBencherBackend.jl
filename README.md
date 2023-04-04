# VambBenchmarks
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://jakobnissen.github.io/VambBenchmarks.jl/dev)
This package is used to benchmark and interactively explore the results of metagenomic binning given a dataset.

## Installation
* Install Julia (ideally using Juliaup: https://github.com/JuliaLang/juliaup)
* Launch Julia: `julia`
* Press `]` to enter package mode. You can exit package mode with backspace.
* In package mode, type `add https://github.com/jakobnissen/VambBenchmarks.jl` to download and install the benchmarking software

## Documentation
Basic usage:
```julia
using VambBenchmarks
ref =  open(i -> Reference(i), "files/ref.json")
bins = open(i -> Binning(i, ref), "files/clusters.tsv")
print_matrix(bins)
```

For more, see See the [documentation](https://jakobnissen.github.io/VambBenchmarks.jl/dev)
