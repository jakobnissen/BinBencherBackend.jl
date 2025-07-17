# BinBencher.jl
BinBencher.jl is the command-line interface (CLI) around the BinBencherBackend, a program used to benchmark the output of a metagenomic binning,
where the ground truth is known.
The CLI program also includes functionality to [create a reference JSON file](@ref references).

You can learn about the input files in the [relevant section of the documentation](@ref input_files).

See more details about the subcommands of BinBencher.jl by navigating the documentation menu on the left.

## Installation of the CLI program
1. Install the [Julia programming language](https://julialang.org/)
2. Make sure `julia` is on your `$PATH`, e.g. test it by running `julia -v`
3. Install BinBencher with the following command:

```shell
julia --startup=no --project=@binbencher -e 'using Pkg; Pkg.add(url="https://github.com/jakobnissen/BinBencher.jl", rev="v0.1.0"); Pkg.precompile(); Pkg.build()'
```

You will now have an executable script called `binbench` in the `bin` subdirectory your Julia home directory.
For example, on my computer the script is at `~/.julia/bin/binbench`.
Add the `bin` folder to your `$PATH` environmental variable to be able run to BinBencher as `binbench` from the shell.

## Quickstart

```shell
$ # Make a reference JSON file
$ binbench makeref ref_outdir --seq-mapping seq_mapping.tsv --seq-fasta contigs.fna \
--genome-directories organism=genomes --tax tax.tsv --tax-ncbi ncbi.tsv

$ # Benchmark using the reference and a binning TSV file.
$ binbench bench bench_outdir reference.json binning.tsv
```
