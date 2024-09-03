# [Benchmarking a binning](@id benchmarking)
Benchmarking is done with by running `binbench bench` from command line.
You will need the following two input files:

1. A [reference JSON file](@ref references).
2. Your bins in a two-column TSV file.

Benchmarking is run with the following command:
```shell
$ binbench bench out_dir reference.json binning.tsv
```

### Command line arguments
* `-s --sep`: Set the binsplit separator to a string (by default, no separator is passed).
  When a separator is passed, every sequence is assumed to be named according to the format
  `<samplename><separator><sequencename>`, where:
  - `<samplename>` uniquely identifies the sample of origin of the sequence
  - `<separator>` is equal to the binsplit separator
  - The separator does not occur in the samplename
  - `<sequencename>` is an identifier that uniquely identifies a sequence within a sample.
  BinBencher will then split every input bin by their sample as given by this format.
  All bin filters will be applied after this binsplitting procedure.
* `--minsize`: Remove bins where the sum of sequence lengths is smaller than this number
* `--minseqs`: Remove bins with fewer sequences than this number
* `--keep-flags`: See `--remove_flags`. Keep only genomes with all these flags set.
* `--remove-flags`: Every genome in the reference has a set of boolean attributes, e.g.
  `virus`, `organism`, `plasmid` etc. This argument takes a comma-separated list of flags,
  then removes every genome that matches the flag.
  E.g. `--remove-flags organism,virus`, all genomes that
  are flagged with _either_ `organism` or `virus` is removed.
* `--recalls`: Comma-separated list of recall thresholds to use when computing the
  matrices of the binning
* `--precisions`: See `--recalls`
* `--intersect`: Allow sequences to be present in multiple bins. Without this flag, BinBencher
  will throw an error if any sequences is detected in multiple bins.

### Output files
* `log.txt`: The log file.
* `bins.json.gz`: This file contains _assembly recall_, _genomic recall_ and _precision_ for every bin,
  with every genome and clade that intersects with the bin (i.e. have nonzero recall).
* `recovery.json`: This JSON file contains the keys `"precisions"` and `"recalls"` with the precision/recall
  thresholds used when binning. It then contains the four arrays `bins_genomic_recall`, `genomes_genomic_recall`,
  `bins_asm_recall` and `genomes_asm_recall`. Each of these arrays contain N matrices, one per taxonomic rank of the reference,
  with the first matrix being the genomes and the last matrix the top taxonomic rank.
  Each matrix has one row per precision threshold and one column per recall threshold.
  It then gives the number of genomes/bins reconstructed at the given precision/recall thresholds.
  The `bins_*` matrices gives the number of bins passing the threshold, whereas the `genomes_` gives the number of genomes
  for which at least 1 bin passes the threshold with that given genome.
  The `*_genomic_*` matrices counts _genomic recall_, whereas the `*_asm_*` matrices count _asm recall_.

### Benchmarking multiple binnings in one go
If you have several binnings, it's faster to benchmark them together.
This can be done like so
```shell
$ binbench bench out_dir reference.json bins1=binning_1.tsv bins2=binning_2.tsv
```
Here, `out_dir` will be populated by the log file and one sub-directory per binning, here called `bin1` and `bins2`.
The sub-directories contain the same content as the output directory when running BinBencher with a single binning,
except it does not contain the log file.