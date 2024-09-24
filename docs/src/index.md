# BinBencher
BinBencher is a tool for analysis and benchmarking of a set of metagenomic assembled genomes (MAGs) against a reference.
It is designed to be used for benchmarking metagenomic binners against a metagenome with a known ground truth,
typically from a simulated metagenome such as the ones from [CAMISIM](https://github.com/CAMI-challenge/CAMISIM),
or alternatively from a laboratory-created mock community. 

BinBencher is [pre-published at bioRxiv](https://www.biorxiv.org/content/10.1101/2024.05.06.592671v1.full.pdf).
Read that paper for background and motivation for this package.

BinBencher consists of two packages:
* `BinBencherBackend.jl`: The Julia library that does the actual analysis.
* `BinBencher.jl`: An application that wraps the backend in a command-line interface (CLI).

Use the side bar to your left to find the documentation for either the library, for the CLI application. Most users probably want to use the CLI application.
