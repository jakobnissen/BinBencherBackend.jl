# [Creating a reference](@id references)

For more information on the reference, see [the description of the JSON file format](@ref refjson).
For a definition of the concepts like _source_ or _genome_, see [the section on nomenclariture](@ref nomenclariture).

At a high level, the reference contains:
* A list of the underlying genomes that are included in the simulation / the mock community.
* A tree that gives the phylogenetic relationship between these genomes
* A list of sequences included in the binning, and to which positions in the genomes they map to (are simulated from).

In situations where the sequences have been simulated, you have perfect access to this information.
However, there may be situations where this knowledge is only approximate - for example, if you simulate
reads from the genomes and then assemble the reads, it might not be possible to know with 100% certainty
where each contig ultimately is sourced from.

Although you can manually create the reference JSON file, I recommend using the `binbench makeref` command.
It takes a number of input files to extract all the information, which I will go through here:

## Genomes
Each genome must be stored in a single fasta file, and the genome takes its name from the FASTA file (stripping the file extension).
The entries in this fasta file are the sources.
These fasta files are stored in some directory.
To set the flags of the genome, have one directory with fasta files per set of flags.

#### Example
Given the following directory structure
```
.
├── collection_a
│   ├── genome1.fna
│   └── genome2.fna
└── collection_b
    └── genome3.fna
```
And where the content of e.g. `genome1.fna` is:
```
>seq1
TAGA
>seq2
TAA
GA
```
You can pass `--genome-directories plasmid+virus=collection_a,organism=collection_b`.
Thus, the genomes `genome1` and `genome2` will be have the flags `virus` and `plasmid`, and `genome3` will have the flag `organism`.
The genome `genome1` will have two sources `seq1` with length 4 and `seq2` with length 5.

## [Sequences](@id seqs)
The sequence (contig) names and lengths are obtained from reading a FASTA file with the sequences.

The mapping position(s), if any, of each sequence is given by a 4-column TSV file with the sequence name, the source name, mapping start and mapping end (1-based, inclusive).

#### Example
Given sequences `seq1`, `seq2` and `seq3`, and the following mapping file
```
sequence	source	start	end
seq1	src1	4	9
seq1	src1	11	15
seq1	src2	1	6
seq3	src3	15	7
```
* `seq1` maps to three locations, two of which are in the same source. Note that the length of the spans (mapping end to start) need not be the same
* `seq2` is not present - it maps to nothing (or its mapping is unknown)
* `seq3` maps exclusively to `src3`. Note that its mapping start is greater than its mapping end. This indicates that `src3` is circular,
  and `seq3` maps across the breakpoint (origin) of the sequence

Mapping locations use 1-based indexing and are inclusive.

## Taxonomy
The taxonomy is specified with a 3-column TSV file with the header `rank\tchild\tparent`.
The first column gives the rank, which may be one of: 'strain', 'species', 'genus', 'family', 'order', 'class', 'phylum', and gives the rank of the child.
The second column gives the name of the child, which is either a genome name (if the rank is `strain`), or else a parent specified in this file.
The last column gives the name of the parent.

#### Example
```
rank	child	parent
strain	genome1	Bacillus subtilis
strain	genome2	Bacillus subtilis
strain	genome3	Clostridium tetani
species	Bacillus subtilis	Bacillus
species	Clostridium tetani	Clostridium
genus	Bacillus	Bacillaceae
genus	Clostridium	Clostridiaceae
family	Bacillaceae	Bacillales
family	Clostridiaceae	Eubacteriales
order	Bacillales	Bacilli
order	Eubacteriales	Clostridia
class	Bacilli		Bacillota
class	Clostridia	Bacillota
```
This specifies the full lineage of all three genomes up to a single universal ancestor of all the genomes (Bacillota, in this case).

#### Using NCBI for easier taxonomy
It may be annoying to specify the full phylogenetic tree of all clades in the taxonomy file.
To ease it, the parent may instead be listed as `id=XXXX`, where `XXXX` is an NCBI taxonomy ID on the parent rank.
BinBencher will then fill up all higher ranks automatically.
To enable this, `--tax-ncbi` must be passed on the command line to a file that contains a list of all known NCBI taxonomies.
A link to this file can be found at the BinBencher github repo.
So, the taxonomy file above can more easily be written as:

```
rank	child	parent
strain	genome1	id=1423
strain	genome2	id=1423
strain	genome3	id=1513
```

You do not need to specify the NCBI tax id for all entries, but can specify them for a subset of entries if necessary.

## Using intermediate JSON files
The genomes, the sequences and the taxonomy are each parsed from their input files, and saved as intermediate JSON files.

For example, after [parsing in information about the sequences](@ref seqs), BinBencher will create a `seqs.json` file.
For subsequent runs, instead of passing the sequences fasta file and the mapping locations with `--seq-fasta` and `--seq-mapping`,
you can pass the JSON file with `--seq-json`.

You can also manually create these intermediate output JSON files, such that you don't need to supply the corresponding input files.
Suppose you want to create a reference, but that your genomes have such complex sets of flags that it is impractical to pass `--genome-directories`.
You can then manually create the genomes intermediate JSON file and pass `--genome-json` instead of `--genome-directories`.
