module VambBenchmarks

using AbstractTrees: AbstractTrees
using JSON3: JSON3
using StructTypes: StructTypes
using LazilyInitializedFields: @lazy, @isinit, @init!, uninit
using SnoopPrecompile: @precompile_all_calls

include("utils.jl")
include("sequence.jl")
include("source.jl")
include("clade.jl")
include("genome.jl")
include("bin.jl")
include("reference.jl")
include("binning.jl")

vector(x) = x isa Vector ? x : vec(collect(x))
imap(f) = x -> Iterators.map(f, x)
ifilter(f) = x -> Iterators.filter(f, x)

@precompile_all_calls let
    ref = open(joinpath(dirname(dirname(pathof(VambBenchmarks))), "files", "ref.json")) do io
        Reference(io)
    end
    bins = open(joinpath(dirname(dirname(pathof(VambBenchmarks))), "files", "clusters.tsv")) do io
        Binning(io, ref)
    end
    print_matrix(IOBuffer(), bins, 1)
end

export Sequence,
    Source,
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
