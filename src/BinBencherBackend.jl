module BinBencherBackend

using AbstractTrees: AbstractTrees
using CodecZlib: GzipDecompressorStream
using JSON: JSON
using StructTypes: StructTypes
using PrecompileTools: @setup_workload, @compile_workload

include("exceptions.jl")
include("utils.jl")
include("flags.jl")
include("sequence.jl")
include("types.jl")
include("source.jl")
include("clade.jl")
include("genome.jl")
include("reference.jl")
include("bin.jl")
include("binning.jl")

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
        n_recovered(bins, 0.6, 0.6)
        n_recovered(bins, 0.6, 0.6; assembly = false)
        n_recovered(bins, 0.6, 0.6; level = 1)
        n_recovered_bounds(bins, 0.65, 0.65)
    end
end

export BinBencherException,
    ReferenceError,
    BinningError,
    Sequence,
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
    genome_size,
    assembly_size,
    shortest_seq_len,
    fraction_assembled,
    is_organism,
    is_virus,
    is_plasmid,
    top_clade,
    gold_standard,
    n_recovered,
    n_recovered_bounds,
    n_passing_bins,
    n_passing_bins_bounds,
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
    is_disjoint,
    adjusted_rand_index

public total_bp

end # module
