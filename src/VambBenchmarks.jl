module VambBenchmarks

using AbstractTrees: AbstractTrees
using JSON3: JSON3
using StructTypes: StructTypes
using LazilyInitializedFields: @lazy, @isinit, @init!, @uninit!, uninit
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
    dir = joinpath(dirname(dirname(pathof(VambBenchmarks))), "files")
    @assert isdir(dir)
    ref = open(joinpath(dir, "ref.json")) do io
        Reference(io)
    end
    bins = open(joinpath(dir, "clusters.tsv")) do io
        Binning(io, ref)
    end
    print_matrix(IOBuffer(), bins)
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
    filter_size,
    f1,
    fscore,
    mrca,
    print_matrix

end # module
