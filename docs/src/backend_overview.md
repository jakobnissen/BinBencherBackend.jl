# BinBencherBackend.jl

This is an ordinary Julia package that the CLI tool BinBencher.jl calls into for its analysis.
The backend is useful for users who want to do more complex analysis not provided by the CLI wrapper.

## Installation
* Install Julia - preferably using `juliaup`: https://github.com/JuliaLang/juliaup
* Launch Julia: `julia`
* Press `]` to enter package mode. You can exit package mode with backspace.
* In package mode, type `add BinBencherBackend` to download and install the benchmarking software

## Quickstart
```julia
using BinBencherBackend
ref =  Reference("files/ref.json")
bins = Binning("files/clusters.tsv", ref)
print_matrix(bins)
```