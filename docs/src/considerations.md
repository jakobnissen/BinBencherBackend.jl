# Benchmarking considerations
### [Adjusted Rand Index is not suitable for benchmarking binnings](@id ari)
BinBencherBackend contains the function [`adjusted_rand_index`](@ref) to compute
this metric.
However, we do not recommend using ARI to evaluate binnings, because it comes with several problems:

##### ARI is invalid if the two clusterings do not contain the same elements. 
That implies that when used to compare a binning against a gold standard binning,
ARI is invalid if the binning does not include every sequence assigned to a genome.
Since most binners do not output every bins containing every single input contig,
this implies ARI cannot be used on the output of most binners

##### ARI cannot handle redundant sequences
Redundant sequences are sequences that map to the same (or overlapping) positions of reference genomes.
For example, if a dataset contain multiple sequence mapping to the same genome,
such that the genome has a coverage from the assembled sequences of 2x,
then it is valid to either output the sequences in one single bin with 2x coverage,
or partitioned in two bins with 1x coverage.
And, if partitioned into two bins, there may be many permutations of input sequences that is equally correct.

When measuring recovered genomes, this is not an issue, as either way will count the genome as being recovered.
In contrast ANI measures against one single target binning and do not have a concept of equivalent or partially equivalent sequences.

While this concern may seem esoteric, [our investigations have shown that the presence of redundant sequences can cause large errors in recall/precision estimates](https://www.biorxiv.org/content/10.1101/2024.05.06.592671v1).

Furthermore, this issue compounds with the problem of microdiversity.
If the reference contains genomes A and B that are 99.9% identitcal on a nucleotide level (ANI),
then the gold standard binning will need to choose between:

1. Distinguishing two genomes with a 99.9% ANI, which is unreasonable at our current level
   of technology, where we can't reliably assemble at that level, or
2. Collapse the two genomes to a single genome, in which all sequences of the genomes become fully redundant with the other. 