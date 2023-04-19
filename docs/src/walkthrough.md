# Example walkthrough
_Note: This is run 2023-04-05 with commit db18e475_

!!! danger
    The documentation on this page uses a large dataset, and which is stored hosted online.
    Hence, the code in this documentation is not tested, and may become outdated.

## Loading the reference
First, let's load the reference:
```jldoctest walk
julia> using VambBenchmarks

julia> ref_path = "/home/jakni/Downloads/vambdata/newref/ref_ptracker_megahit_Skin.json";

julia> bins_path = "/home/jakni/Downloads/vambdata/newbins/megahit_skin.tsv";

julia> ref = open(i -> Reference(i), ref_path)
Reference
  Genomes:    2341
  Sequences:  959970
  Ranks:      9
  Seq length: 200
  Assembled:  25.0 %
```

This gives us a few statistics about the reference:
* Number of genomes
* Number of sequences
* Number of taxonomic ranks (strain, species, genus...)
* Length of shortest sequence
* Total length of genomes that are assembled

The genomes here contain both plasmids and organisms, and the sequence length of 200 bp is too short.
Let's filter the reference using [`subset!`](@ref) to only retain organisms, and sequences of length 1500 or more:

```jldoctest walk
julia> subset!(ref;
           genomes = is_organism,
           sequences = s -> length(s) >= 1500
       )
Reference
  Genomes:    1394
  Sequences:  118267
  Ranks:      9
  Seq length: 1500
  Assembled:  16.1 %
```

!!! note
    The function [`subset!`](@ref) will mutate the reference, whereas the function [`subset`](@ref)
    will create a new independent reference. At the moment, the latter is much slower.

You can see we removed 7/8th of all sequences, but only 1/3rd of the total assembly length.
We also removed almost half of all genomes (namely, all those that wasn't organisms).

## Genomes
We can explore the genomes contained in the reference with the `genomes` function,
which returns an iterable of `Genome` (in this case, a `Set`):

```jldoctest walk; filter = r"\s+Genome\([A-Za-z0-9\.]+\)"
julia> genomes(ref)
Set{Genome} with 1394 elements:
  Genome(OTU_97.34832.0)
  Genome(OTU_97.33610.0)
  Genome(OTU_97.39616.0)
  Genome(OTU_97.23938.0)
  Genome(OTU_97.24774.0)
  Genome(OTU_97.36286.0)
  Genome(OTU_97.6382.0)
  Genome(OTU_97.33661.1)
  Genome(OTU_97.37247.0)
  Genome(OTU_97.1829.0)
  Genome(OTU_97.6602.0)
  Genome(OTU_97.44820.0)
  Genome(OTU_97.24834.0)
  Genome(OTU_97.44856.1)
  Genome(OTU_97.3377.0)
  Genome(OTU_97.19529.1)
  Genome(OTU_97.32371.0)
  Genome(OTU_97.45083.0)
  Genome(OTU_97.15121.1)
  ⋮
```

Let's look at a `Genome` in more detail:
```jldoctest walk
julia> genome, genome2 = sort!(collect(genomes(ref)); by=i -> i.name);

julia> genome
Genome "OTU_97.10046.0"
  Parent:        "Staphylococcus capitis"
  Genome size:   2474232
  Assembly size: 9239 (0.4 %)
  Sources:       1
  Flags:         1 (organism)
```

The _flags_ can be extracted with the `flags(genome)` function - each genome contains zero or more flags:

```jldoctest walk
julia> flags(genome)
FlagSet with 1 element:
  VambBenchmarks.Flags.organism
```

... in this case, this genome is an organism as opposed to a plasmid.
You can see all possible flags with `instances(Flags.Flag)`.

See also the helper functions [`is_organism`](@ref), [`is_virus`](@ref) and [`is_plasmid`](@ref)

## Sources
The `genome` has one source - let's look at that

```jldoctest walk
julia> source = only(genome.sources)
Source "CP007601.1"
genome:          Genome(OTU_97.10046.0)
  Length:        2474232
  Assembly size: 9239
  Sequences:     2
```

A `Source` is one of the genomic sequences that genomes are composed of. In the reference database where our genome
was gotten from, it was assembled into a single contig called "CP001781.1" of 2.47 Mbp - presumably the complete, circular genome.
We can see that 2 sequences map to this source, covering about 9.2 kbp.

We can get the sequences mapping to this source:

```jldoctest walk
julia> source.sequences
2-element Vector{Tuple{Sequence, UnitRange{Int64}}}:
 (Sequence("S13C8728", 2472), 2069906:2072374)
 (Sequence("S13C26910", 7349), 2085051:2091820)
```

Where, e.g. the first entrance tells us that the sequence "S13C8728" with a length of 2472 maps to positions 2069906:2072374 (both inclusive).

## Clades
Genomes are organised into a taxonomic hierarchy.
We can find the immediate parent of a genome by accessing the field `genome.parent`.
Let's look at another genome:

```jldoctest walk
julia> genome2
Genome "OTU_97.10083.0"
  Parent:        "Moraxella bovoculi"
  Genome size:   2220786
  Assembly size: 0 (0.0 %)
  Sources:       1
  Flags:         1 (organism)

julia> clade = genome2.parent
Species "Moraxella bovoculi", 6 genomes
├─ Genome(OTU_97.36230.1)
├─ Genome(OTU_97.11114.0)
├─ Genome(OTU_97.10083.0)
├─ Genome(OTU_97.9999.0)
├─ Genome(OTU_97.36230.0)
└─ Genome(OTU_97.30020.0)
```

The parent is an instance of a `Clade`.
`Clade`s are at a specific rank: Rank 1 for species, 2 for genus, 3 for family, etc.
Every clade has one or more children: These are the clades one rank lower.
Conceptually, rank zero corresponds to `Genome`s (OTUs, for this reference dataset)

```jldoctest walk
julia> clade.children
6-element Vector{Genome}:
 Genome(OTU_97.36230.1)
 Genome(OTU_97.11114.0)
 Genome(OTU_97.10083.0)
 Genome(OTU_97.9999.0)
 Genome(OTU_97.36230.0)
 Genome(OTU_97.30020.0)
```

We can find the most recent common ancestor (MRCA) of `genome` and `g2` like this:
```jldoctest walk
julia> mrca(genome, genome2)
Domain "Bacteria", 1394 genomes
├─ Phylum "Deinococcus-Thermus", 12 genomes
│  └─ Class "Deinococci", 12 genomes
│     └─ Order "Deinococcales", 12 genomes
│        ⋮
│
├─ Phylum "Tenericutes", 5 genomes
│  └─ Class "Mollicutes", 5 genomes
│     └─ Order "Mycoplasmatales", 5 genomes
│        ⋮
│
├─ Phylum "Fusobacteria", 21 genomes
│  └─ Class "Fusobacteriia", 21 genomes
│     └─ Order "Fusobacteriales", 21 genomes
│        ⋮
│
├─ Phylum "Bacteroidetes", 85 genomes
│  ├─ Class "Sphingobacteriia", 1 genome
│  │  └─ Order "Sphingobacteriales", 1 genome
│  │     ⋮
│  │
│  ├─ Class "Bacteroidia", 47 genomes
│  │  └─ Order "Bacteroidales", 47 genomes
│  │     ⋮
│  │
⋮
```

They are very distantly related, so the domain "Bacteria", one of the highest ranked `Clade`s, are their most recent common ancestor.

The top clade can be found with the `top_clade(ref)` function.
In our case, the dataset only consists of bacteria and plasmids, and we filtered away all plasmids, so the universal ancestor
only has one child, namely the domain Bacteria.

## Binnings
A `Binning` is a set of bins benchmarked against a reference.
We can load a set of Vamb bins and turn it into a `Binning` object like this:

```jldoctest walk
julia> binning = open(bins_path) do io
           Binning(io, ref; binsplit_separator='C')
       end
Binning
  Reference
    Genomes:    1394
    Sequences:  118267
    Ranks:      9
    Seq length: 1500
    Assembled:  16.1 %
  Bins:        14535
  NC genomes:  10
  Precisions: [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99]
  Recalls:    [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99]
  Recoverable genomes: [325, 281, 235, 193, 142, 98, 51, 21, 0]
  Reconstruction (assemblies):
    P\R   0.3  0.4  0.5  0.6  0.7  0.8  0.9 0.95 0.99
    0.3   268  244  229  201  167  111   52   20    4
    0.4   247  223  212  188  157  102   50   20    4
    0.5   230  207  196  175  146   95   45   17    4
    0.6   212  192  181  160  135   87   39   14    4
    0.7   191  173  162  143  122   80   36   13    3
    0.8   161  145  136  118  101   71   33   13    3
    0.9   112   98   90   79   71   55   29   12    3
    0.95   94   84   80   70   63   50   26   12    3
    0.99   69   64   62   57   51   39   23   10    2
```

A wealth of information is readily available:
* `binning.ref` gives the underlying `Reference`
* `binning.recalls` and `binning.precision` gives the recall/precision thresholds used in benchmarking
* You can get the number of genomes that are assembled at the various recall levels with `recoverable_genomes`:
  This sets an upper limit of how many genomes can be reconstructed given a hypothetical perfect binning

```jldoctest walk
julia> println(binning.recoverable_genomes)
[325, 281, 235, 193, 142, 98, 51, 21, 0]
```

The function `print_matrix` will display the number of recovered genomes/assemblies.
It takes two optional keyword: `level`, the taxonomic rank (defaults to 0, meaning strain level),
and `assembly` which defaults to `true`. If set to `false`, benchmark number of recovered genomes,
not number of recovered assemblies.

```jldoctest walk
julia> print_matrix(binning; level=1, assembly=false)
P\R   0.3  0.4  0.5  0.6  0.7  0.8  0.9 0.95 0.99
0.3    78   74   72   66   52   34   10    1    0
0.4    77   73   71   64   50   33   10    1    0
0.5    76   72   70   63   49   32   10    1    0
0.6    75   71   69   62   48   31   10    1    0
0.7    72   68   66   59   46   30   10    1    0
0.8    72   68   66   59   46   30   10    1    0
0.9    72   68   66   59   46   30   10    1    0
0.95   70   66   64   56   44   30   10    1    0
0.99   56   52   51   46   37   25    9    1    0
```

## Bins
The `Binning` object obviously contains our bins.
Let's pick a particularly good bin:

```jldoctest walk
julia> bin = binning.bins[261]
Bin "S13Cvae_479"
  Sequences: 47
  Breadth:   2482216
  Intersecting 1 genome
```

The "breadth" here is the sum of the length of its sequences.
`bin.sequences` gets an iterable of its sequences:

```jldoctest walk
julia> collect(bin.sequences)
47-element Vector{Sequence}:
 Sequence("S13C2748", 73726)
 Sequence("S13C13383", 18077)
 Sequence("S13C59849", 95631)
 Sequence("S13C29443", 13209)
 Sequence("S13C56904", 44356)
 Sequence("S13C14948", 9756)
 Sequence("S13C37848", 19070)
 Sequence("S13C2782", 22054)
 Sequence("S13C40546", 46579)
 Sequence("S13C12978", 15372)
 ⋮
 Sequence("S13C33782", 34963)
 Sequence("S13C27111", 21961)
 Sequence("S13C10862", 55887)
 Sequence("S13C54667", 69970)
 Sequence("S13C19757", 199299)
 Sequence("S13C64535", 73792)
 Sequence("S13C2686", 20376)
 Sequence("S13C62986", 19619)
 Sequence("S13C14459", 39940)
```

The "Intersecting 1 genome" means that the sequences map to 1 genome.
We can get that with the function [`intersecting`](@ref), then get the precision/recall
with [`recall_precision`](@ref):

```jldoctest walk
julia> intersecting_genome = only(intersecting(bin))
Genome "OTU_97.9674.1"
  Parent:        "Corynebacterium humireducens"
  Genome size:   2714910
  Assembly size: 2575299 (94.9 %)
  Sources:       1
  Flags:         1 (organism)

julia> recall_precision(only(intersecting(bin)), bin)
(0.9614992278566489, 1.0)
```

Or do the same for a higher clade - let's say a genus:
```jldoctest walk
julia> genus = only(Iterators.filter(i -> i.rank == 2, intersecting(Clade, bin)));

julia> recall_precision(genus, bin)
(0.9614992278566489, 1.0)
```
