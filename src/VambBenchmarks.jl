module VambBenchmarks

using AbstractTrees: AbstractTrees
using JSON3: JSON3
using StructTypes: StructTypes
using SnoopPrecompile: @precompile_all_calls

include("utils.jl")
include("sequence.jl")
include("clade.jl")
include("genome.jl")
include("bin.jl")
include("reference.jl")
include("binning.jl")

vector(x) = x isa Vector ? x : vec(collect(x))
imap(f) = x -> Iterators.map(f, x)
ifilter(f) = x -> Iterators.filter(f, x)

include("workload.jl")

# @precompile_all_calls exercise()

export Sequence,
    Clade,
    Genome,
    Bin,
    Reference,
    Binning,
    nbins,
    nseqs,
    ngenomes,
    top_clade,
    gold_standard,
    f1,
    fscore,
    mrca,
    print_matrix

end # module
