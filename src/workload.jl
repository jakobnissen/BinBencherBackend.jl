function exercise()
    ref = open(joinpath(dirname(dirname(pathof(VambBenchmarks))), "files", "ref.json")) do io
        Reference(io)
    end
    bins = open(joinpath(dirname(dirname(pathof(VambBenchmarks))), "files", "clusters.tsv")) do io
        Binning(io, ref; binsplit_separator='c')
    end
    print_matrix(IOBuffer(), bins, 1)
end
