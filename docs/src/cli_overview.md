# BinBencher.jl
BinBencher.jl is the command-line interface (CLI) around the BinBencherBackend, a program used to benchmark the output of a metagenomic binning,
where the ground truth is known.
The CLI program also includes functionality to [create a reference JSON file](@ref references).

You can learn about the input files in the [relevant section of the documentation](@ref input_files).

See more details about the subcommands of BinBencher.jl by navigating the documentation menu on the left.

## Installation of the CLI program
1. Install the [Julia programming language](https://julialang.org/)
2. Launch the Julia REPL from a terminal by typing `$ julia`
3. Enter Pkg mode in the REPL by typing `]`. You will see the prompt change from `julia>` to e.g. `(@v1.12) pkg>` depending on the Julia version.
4. Create a new globally accessible project (i.e. virtual environment) by typing `activate @binbencher` in the pkg REPL.
5. Add BinBencher to the environment with `add https://github.com/jakobnissen/BinBencher.jl`
6. In the pkg repl, type `build`

You will now have an executable script called `binbench` in the `bin` subdirectory your Julia home directory, e.g. on my computer it's at `/home/jakni/.julia/bin/binbench`.
Add the `bin` folder to your PATH environmental variable to be able run to BinBencher as `binbench` from the shell.

!!! info
    These installation instruction are convoluted, because, as of September 2024, Julia does not have a concept of installable
    applications.
    There is [recent work to improve this](https://github.com/JuliaLang/Pkg.jl/pull/3772), so hopefully BinBencher will be easier
    to install and run by the time the work is done some time in 2025.

## Quickstart

```shell
$ # Make a reference JSON file
$ binbench makeref ref_outdir --seq-mapping seq_mapping.tsv --seq-fasta contigs.fna \
--genome-directories organism=genomes --tax tax.tsv --tax-ncbi ncbi.tsv

$ # Benchmark using the reference and a binning TSV file.
$ binbench bench bench_outdir reference.json binning.tsv
```
