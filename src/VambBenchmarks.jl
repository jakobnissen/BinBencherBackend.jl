module VambBenchmarks

using AbstractTrees: AbstractTrees
using JSON3: JSON3
using StructTypes: StructTypes

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

precompile(Reference, (IOStream,))
precompile(Binning, (IOStream, Reference))

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
    f1,
    fscore,
    mrca,
    print_matrix

end # module
