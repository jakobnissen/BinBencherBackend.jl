## Nomenclariture of BinBencher
### [Basic nomenclariture](@id nomenclariture)
* A __sequence__ is the sequences that are to be clustered by the binner.
  I also sometimes call this a __contig__, because in most applications, binners are used to bin contigs.
* A __genome__ is a target genome that should be reconstructed by the binner.
  It can be a virus, organism, plasmid or any other genetic entity.
  Every genome have one or more sources, zero or more flags, and one parent clade.
* A __flag__ marks the certaincy about a boolean attribute of a genome, like "is this a virus?".
* A __source__ are the DNA sequences that genomes are composed of.
  E.g. for a bacteria, its sources are all its chromosomes. For plasmids, its only source is its full-length sequence. 
  Sources are typically the full genome reference sequences from e.g. NCBI's RefSeq.
  Sequences map to zero or more sources at particular spans, i.e. a contiguous range from a starting to a stopping position.
* A __clade__ contain one or more genomes, or clades, known as its children.
  Clades containing genomes are rank 1, and clades containing rank N clades are rank N+1 clades.
  By extension, genomes have rank 0.
  A clade never has children of mixed ranks.
  All genomes descend from a chain of exactly N ranks of clades, where N > 0.
* A __bin__ is a set of sequences created by the binner. Every bin is benchmarked against all genomes and clades in the reference.
* A __reference__ is composed of:
    - The set of genomes__ to benchmark against,
      which all ultimately descend from a single common clades that links all genomes together.
    - Its set of sequences each with zero or more mappings to its genomes' sources.
* A __binning__ is a set of bins benchmarked against a reference, where every bin is a subset of the reference's sequences.

In the equations below, we denote genomes ``G``.
Genomes can be considered disjoint sets of mapping positions (i.e. genomic positions), where each position is a basepair
in a source of the genome.
The total set of positions in the reference is ``Y = \cup_{G}``.
Let ``X`` be the set of sequences ``S`` to be binned.
Each sequence has a length ``L_S`` and can be considered as a set of mapping positions ``S \subseteq Y``.
Note that the cardinality of this set ``|S|`` need not be equal to the sequence's length ``L_S``.

If we have a set of sequences ``x \subseteq X`` and a set of mapping positions ``y \subseteq Y``, let us define
``x \doublecap y := \{S \in x | S \cap y \neq \emptyset \}``, i.e. the set of sequences in ``x``
that have mapping positions in ``y``.

A genome's _assembly_ ``A_G = G \cap \cup_X`` are all mapping positions in the genome covered by any sequence in the dataset.

A bin ``B`` is a set of sequences ``B \subseteq X``. For any bin-genome pair ``\{B,G\}``, we define:
* The true positives ``TP_{\{B,G\}} = |\cup_{S \in B} S \cap G|``, the number of mapping positions in ``G`` that any sequence of ``B`` maps to.
* The false positives ``FP_{\{B,G\}} = \sum_{S \in G \setminus (B \doublecap G)} L_S`` is the sum of lengths of sequences in ``B`` that does not map to any position in ``G``.
* The false _genomic_ negatives ``FG_{\{B,G\}} = |G| - TP_{\{B,G\}}`` is all positions in ``G`` not covered by any sequence in ``B``.
* The false _assembly_ negatives ``FA_{\{B,G\}} = |A_G| - TP_{\{B,G\}}`` is the assembly of G not covered by any sequence in ``B``.

## Derived metrics
BinBencher computes a bunch of different metrics. The most central ones are given here. More may be added in the future.

* The precision of a bin-genome pair is ``\frac{TP_{\{B,G\}}}{TP_{\{B,G\}} + FP_{\{B,G\}}}``
* The _genomic recall_ is ``\frac{TP_{\{B,G\}}}{TP_{\{B,G\}} + FG_{\{B,G\}}}``
* The _assembly recall_ is ``\frac{TP_{\{B,G\}}}{TP_{\{B,G\}} + FA_{\{B,G\}}}``

Note that
* These metrics are only defined for bin-genome pairs. It is not sensical to talk about the precision of a genome without specifying with respect to which bin.
* The assembly recall is never lower than the genomic recall.

The most useful measures are the number of recovered bins or genomes at a given recall/precision threshold.
This counts the total number of genomes (bins, respectively) which, when paired with any bin (genomes, resp.) has a recall and precision over the threhold.
Both these measures can be computed for assembly or genomic recall.
