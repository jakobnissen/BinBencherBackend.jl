# BinBencherBackend
BinBencherBackend.jl is a package for efficient benchmarking and interactive exploration of a set of metagenomic assembled genomes (MAGs) against a reference.
This is designed to be used for benchmarking metagenomic binners against a simulated metagenome.

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

## Concepts
* A `Sequence` is a sequence (e.g. contig) clustered by the binner
* A `Genome` is a target genome that should be reconstructed by the binner.
  It can be a virus, organism, plasmid etc. Every `Genome` have several `Source`s, and one parent `Clade`.
* A `Flag` marks the certaincy about a boolean attribute of a genome, like "is this a virus?".
* `Source`s are the sequences that `Genome`s are composed of.
  These are typically the reference genome sequences originally obtained by assembly of a purified genome (e.g. clonal colony).
  `Sequence`s map to zero or more `Source`s at particular _spans_, i.e. locations.
* A `Clade` contain one or more `Genome`s or `Clade`s. Clades containing genomes are rank 1, and clades containing rank N clades are rank N+1 clades.
  All genomes descend from a chain of exactly N ranks of clades, where N > 0.
* A `Bin` is a set of `Sequence`s created by the binner. Every bin is benchmarked against all genomes and clades in the reference.
* A `Reference` is composed of:
    - The _genomes_, a set of `Genome`s, each with a set of `Source`s and `Flags`
    - The _taxmaps_, a full set of Clades that encompasses every `Genome` at N ranks (where N > 0)
    - The _sequences_, a list of `Sequence`s, each with zero or more mappings to `Source`s.
* A `Binning` is a set of `Bin`s benchmarked against a `Reference`

See the Reference in the left sidebar.
