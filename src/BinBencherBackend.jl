module BinBencherBackend

using AbstractTrees: AbstractTrees
using CodecZlib: GzipDecompressorStream
using JSON3: JSON3
using StructTypes: StructTypes
using LazilyInitializedFields: @lazy, @isinit, @init!, @uninit!, uninit
using PrecompileTools: @setup_workload, @compile_workload

include("utils.jl")
include("flags.jl")
include("sequence.jl")
include("source.jl")
include("clade.jl")
include("genome.jl")
include("reference.jl")
include("bin.jl")
include("binning.jl")

vector(x)::Vector = x isa Vector ? x : vec(collect(x))::Vector
imap(f) = x -> Iterators.map(f, x)
ifilter(f) = x -> Iterators.filter(f, x)

@setup_workload begin
    dir = joinpath(dirname(dirname(pathof(BinBencherBackend))), "files")
    refpath = joinpath(dir, "ref.json")
    binpath = joinpath(dir, "clusters.tsv")

    @compile_workload begin
        ref = open(refpath) do io
            Reference(io)
        end
        subset!(ref; sequences = Returns(true), genomes = Returns(true))
        gold_standard(ref)
        bins = open(binpath) do io
            Binning(io, ref)
        end
        print_matrix(IOBuffer(), bins)
        n_recovered(bins, 0.4, 0.2)
        n_recovered(bins, 0.4, 0.2; assembly = false)
        n_recovered(bins, 0.4, 0.2; level = 1)
    end
end

export Sequence,
    flags,
    Flag,
    Flags,
    FlagSet,
    Source,
    Clade,
    Genome,
    Bin,
    Reference,
    Binning,
    n_bins,
    n_seqs,
    intersecting,
    recall_precision,
    genomes,
    is_organism,
    is_virus,
    is_plasmid,
    top_clade,
    gold_standard,
    n_recovered,
    n_passing_bins,
    subset,
    subset!,
    f1,
    passes_f1,
    passes_recall_precision,
    recalls_precisions,
    fscore,
    mrca,
    ancestors,
    descends_from,
    print_matrix,
    is_disjoint

end # module
