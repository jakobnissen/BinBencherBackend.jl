module VambBenchmarks

using AbstractTrees: AbstractTrees
using JSON3: JSON3
using StructTypes: StructTypes
using LazilyInitializedFields: @lazy, @isinit, @init!, @uninit!, uninit
using SnoopPrecompile: @precompile_all_calls, @precompile_setup

include("utils.jl")
include("flags.jl")
include("sequence.jl")
include("source.jl")
include("clade.jl")
include("genome.jl")
include("bin.jl")
include("reference.jl")
include("binning.jl")

vector(x)::Vector = x isa Vector ? x : vec(collect(x))
imap(f) = x -> Iterators.map(f, x)
ifilter(f) = x -> Iterators.filter(f, x)

@precompile_setup let
    dir = joinpath(dirname(dirname(pathof(VambBenchmarks))), "files")
    refpath = joinpath(dir, "ref.json")
    binpath = joinpath(dir, "clusters.tsv")

    @precompile_all_calls let
        ref = open(refpath) do io
            Reference(io; min_seq_length=10)
        end
        gold_standard(ref)
        bins = open(binpath) do io
            Binning(io, ref)
        end
        print_matrix(IOBuffer(), bins)
    end
end

export Sequence,
    flags,
    Flags,
    FlagSet,
    Source,
    Clade,
    Genome,
    Bin,
    Reference,
    Binning,
    nbins,
    nseqs,
    intersecting,
    recall_precision,
    genomes,
    is_organism,
    is_virus,
    is_plasmid,
    top_clade,
    gold_standard,
    subset,
    subset!,
    f1,
    fscore,
    mrca,
    print_matrix

end # module
