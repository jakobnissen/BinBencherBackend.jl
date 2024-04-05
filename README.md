# BinBencherBackend
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://jakobnissen.github.io/BinBencherBackend.jl/dev)
[![Latest Release](https://img.shields.io/github/release/jakobnissen/BinBencherBackend.jl.svg)](https://github.com/jakobnissen/BinBencherBackend.jl/releases/latest)

This package is used to benchmark and interactively explore the results of metagenomic binning given a dataset.
This is the Julia backend of the command-line tool [BinBencher [work in progress]](https://github.com/jakobnissen/BinBencher.jl).

## Installation
* Install Julia (ideally using Juliaup: https://github.com/JuliaLang/juliaup) or from the official website `www.julialang.org`
* Use the Julia package manager to install this package: `] add BinBencherBackend`

## Documentation
Basic usage:
```julia
using BinBencherBackend
ref = Reference("files/ref.json")
bins = Binning("files/clusters.tsv", ref)
print_matrix(bins)
```

For more, see the [documentation](https://jakobnissen.github.io/BinBencherBackend.jl/dev)
