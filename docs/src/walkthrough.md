```@meta
DocTestSetup = quote
    using BinBencherBackend

    (path_to_ref_file, path_to_bins_file) = let
        dir = joinpath(Base.pkgdir(BinBencherBackend), "files")
        (joinpath(dir, "ref.json"), joinpath(dir, "clusters.tsv"))
    end
end
```

# Example walkthrough

## Loading the reference
First, let's load the reference:
```jldoctest walk
julia> using BinBencherBackend

julia> ref = Reference(path_to_ref_file)
Reference
  Genomes:    3
  Sequences:  11
  Ranks:      3
  Seq length: 10
  Assembled:  61.9 %
```

This gives us a few statistics about the reference:
* Number of genomes
* Number of sequences
* Number of taxonomic ranks (strain, species, genus...)
* Length of shortest sequence
* Total length of genomes that are assembled

The genomes here contain both plasmids and organisms.
Let's filter the reference using [`subset!`](@ref) to only retain organisms, and sequences of length 10 or more:

```jldoctest walk
julia> subset!(ref;
           genomes = is_organism,
           sequences = s -> length(s) >= 10
       )
Reference
  Genomes:    2
  Sequences:  11
  Ranks:      3
  Seq length: 10
  Assembled:  91.3 %
```

!!! note
    The function [`subset!`](@ref) will mutate the reference, whereas the function [`subset`](@ref)
    will create a new independent reference. At the moment, the latter is much slower.

We removed a single genome, namely one labeled as virus.

## Genomes
We can explore the genomes contained in the reference with the `genomes` function,
which returns an iterable of `Genome` (in this case, a `Set`):

```jldoctest walk; filter = r"\s+Genome\([A-Za-z0-9\.]+\)"
julia> genomes(ref)
Set{Genome} with 2 elements:
  Genome(gA)
  Genome(gB)
```

Let's look at a `Genome` in more detail:
```jldoctest walk
julia> genome, genome2 = genomes(ref);

julia> genome
Genome "gA"
  Parent:        "D"
  Genome size:   100
  Assembly size: 88 (88.0 %)
  Sources:       1
  Flags:         1 (organism)
```

The _flags_ can be extracted with the `flags(genome)` function - each genome contains zero or more flags:

```jldoctest walk
julia> flags(genome)
FlagSet with 1 element:
  BinBencherBackend.Flags.organism
```

... in this case, this genome is an organism as opposed to a plasmid or virus.
You can see all possible flags with `instances(Flags.Flag)`.

See also the helper functions [`is_organism`](@ref), [`is_virus`](@ref) and [`is_plasmid`](@ref)

## Sources
The `genome` has one source - let's look at that

```jldoctest walk
julia> source = only(genome.sources)
Source "subjA1"
genome:          Genome(gA)
  Length:        100
  Assembly size: 88
  Sequences:     6
```

A `Source` is one of the genomic sequences that genomes are composed of.
This is distinct from the assembled sequences that we will be binning - a `Source` represents the reference sequence,
typically the full genome, assembled from a sequencing run on a clonal colony.
For this genome, we can see it has a length of 100 bp, and that 5 sequences map to this source, covering 88 bp.

We can get the sequences mapping to this source:

```jldoctest walk
julia> source.sequences
6-element Vector{Tuple{Sequence, Tuple{Int64, Int64}}}:
 (Sequence("s1", 25), (5, 29))
 (Sequence("s1", 25), (10, 34))
 (Sequence("s2", 40), (1, 40))
 (Sequence("s3", 50), (51, 98))
 (Sequence("s7", 20), (21, 40))
 (Sequence("s8", 25), (2, 26))
```

Where, e.g. the first entrance tells us that the sequence "s2" with a length of 40 maps to positions 1:40 (both inclusive).

## Clades
Genomes are organised into a taxonomic hierarchy.
We can find the immediate parent of a genome by accessing the field `genome.parent`.
Let's look at another genome:

```jldoctest walk
julia> genome2
Genome "gB"
  Parent:        "D"
  Genome size:   50
  Assembly size: 49 (98.0 %)
  Sources:       2
  Flags:         1 (organism)

julia> clade = genome2.parent
Species "D", 2 genomes
├─ Genome(gA)
└─ Genome(gB)
```

The parent is an instance of a `Clade`.
`Clade`s are at a specific rank: Rank 1 for species, 2 for genus, 3 for family, etc.
Every clade has one or more children: These are the clades one rank lower.
Conceptually, rank zero corresponds to `Genome`s (OTUs, for this reference dataset)

```jldoctest walk
julia> clade.children
2-element Vector{Genome}:
 Genome(gA)
 Genome(gB)
```

We can find the most recent common ancestor (MRCA) of `genome` and `genome2` like this:
```jldoctest walk
julia> mrca(genome, genome2)
Species "D", 2 genomes
├─ Genome(gA)
└─ Genome(gB)
```

They are very distantly related, so the domain "Bacteria", one of the highest ranked `Clade`s, are their most recent common ancestor.

The top clade can be found with the `top_clade(ref)` function, which is the universal ancestor of all clades in the reference.

## Binnings
A `Binning` is a set of bins benchmarked against a reference.
We can load a set of Vamb bins and turn it into a `Binning` object like this:

```jldoctest walk
julia> binning = Binning(path_to_bins_file, ref)
Binning
  Reference
    Genomes:    2
    Sequences:  11
    Ranks:      3
    Seq length: 10
    Assembled:  91.3 %
  Bins:        6
  NC genomes:  0
  Mean bin genome   R/P/F1: 0.51 / 1.0 / 0.672
  Mean bin assembly R/P/F1: 0.546 / 1.0 / 0.704
  Precisions: [0.6, 0.7, 0.8, 0.9, 0.95, 0.99]
  Recalls:    [0.6, 0.7, 0.8, 0.9, 0.95, 0.99]
  Reconstruction (assemblies):
    P\R   0.6  0.7  0.8  0.9 0.95 0.99
    0.6     1    0    0    0    0    0
    0.7     1    0    0    0    0    0
    0.8     1    0    0    0    0    0
    0.9     1    0    0    0    0    0
    0.95    1    0    0    0    0    0
    0.99    1    0    0    0    0    0
```

A wealth of information is readily available:
* `binning.ref` gives the underlying `Reference`
* `binning.recalls` and `binning.precisions` gives the recall/precision thresholds used in benchmarking

The function `print_matrix` will display the number of recovered genomes/assemblies.
It takes two optional keyword: `level`, the taxonomic rank (defaults to 0, meaning strain level),
and `assembly` which defaults to `true`. If set to `false`, benchmark number of recovered genomes,
not number of recovered assemblies.

```jldoctest walk
julia> print_matrix(binning; level=1, assembly=false)
P\R   0.6  0.7  0.8  0.9 0.95 0.99
0.6     1    0    0    0    0    0
0.7     1    0    0    0    0    0
0.8     1    0    0    0    0    0
0.9     1    0    0    0    0    0
0.95    1    0    0    0    0    0
0.99    1    0    0    0    0    0
```

You can also get the number of genomes or assemblies reconstructed at a given
precision/recall level directly with `n_recovered`:

```jldoctest walk
julia> n_recovered(binning, 0.6, 0.7; assembly=true)
1

julia> n_recovered(binning, 0.66, 0.91; level=1)
0
```

## Bins
The `Binning` object obviously contains our bins.
Let's pick a random bin:

```jldoctest walk
julia> bin = binning.bins[4]
Bin "C4"
  Sequences: 3
  Breadth:   55
  Intersecting 2 genomes
```

The "breadth" here is the sum of the length of its sequences.
`bin.sequences` gets an iterable of its sequences:

```jldoctest walk
julia> collect(bin.sequences)
3-element Vector{Sequence}:
 Sequence("s5", 25)
 Sequence("s6", 10)
 Sequence("s7", 20)
```

The "Intersecting 2 genomes" means that the sequences map to 2 different genomes - the only two in the reference.
We can get that with the function [`intersecting`](@ref), then get the precision/recall
with [`recall_precision`](@ref):

```jldoctest walk
julia> Set(intersecting(bin)) == genomes(ref)
true

julia> recall_precision(genome2, bin)
(recall = 0.6122448979591837, precision = 1.0)
```

Or do the same for a higher clade - let's say a genus. In this case, we get the same result.
```jldoctest walk
julia> genus = only(Iterators.filter(i -> i.rank == 2, intersecting(Clade, bin)));

julia> recall_precision(genus, bin)
(recall = 0.6122448979591837, precision = 1.0)
```
