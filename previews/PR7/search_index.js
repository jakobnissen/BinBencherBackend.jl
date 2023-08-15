var documenterSearchIndex = {"docs":
[{"location":"reference/#Reference","page":"Reference","title":"Reference","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"Modules = [VambBenchmarks]\nOrder   = [:type, :function]","category":"page"},{"location":"reference/#VambBenchmarks.Bin","page":"Reference","title":"VambBenchmarks.Bin","text":"Bin(name::AbstractString, sequences, targets)\n\nBins each represent a bin created by the binner. Conceptually, they are simply a set of Sequence with a name attached. Practically, every Bin is benchmarked against all Genomes and Clades of a given Reference, so each Bin stores data about its intersection with every genome/clade, e.g. its purity and recall.\n\nLike Sources, Bins also have an assembly size for a given Genome. This is the number of base pairs in the genomes covered by any sequence in the Bin, which is always a subset of the genome's assembly size.\n\nBenchmark statistics for a Bin/Genome can be done with either assemblies or genomes as the ground truth.\n\nTrue positives (TP) are defined as the sum of assembly sizes over all sources in the genome\nFalse positives (FP) are the sum of length of sequences in the bin not mapping to the genome\nFalse negatives (FN) is either the genome assembly size or genome size minus TP.\n\nFor Bin/Clade pairs B/C, recall is the maximal recall of B/Ch for all children Ch of C. Precision is the sum of lengths of sequences mapping to any child of the clade divided by the sum of lengths of all sequences in the bin.\n\nSee also: Binning, Genome, Clade\n\nExamples\n\njulia> bin = first(binning.bins)\nBin \"C1\"\n  Sequences: 2\n  Breadth:   65\n  Intersecting 1 genome\n\njulia> first(bin.sequences)\nSequence(\"s1\", 25)\n\njulia> f1(first(ref.genomes), bin)\n0.625\n\n\n\n\n\n","category":"type"},{"location":"reference/#VambBenchmarks.Binning","page":"Reference","title":"VambBenchmarks.Binning","text":"Binning(::IO, ::Reference; kwargs...)\n\nA Binning represents a set of Bins benchmarked against a Reference. Binnings can be created given a set of Bins and a Reference, where the bins may potentially be loaded from a .tsv file. The field binning.recoverable_genomes shows the maximal number of recoverable genomes at given recall levels given perfect binning. The fields recovered_asms and recovered_genomes are used for benchmarking, these are normally output using the print_matrix function.\n\nSee also: print_matrix, Bin, Reference\n\nExamples\n\njulia> bins = gold_standard(ref);\n\njulia> bins isa Binning\ntrue\n\njulia> VambBenchmarks.n_nc(binning)\n0\n\nExtended help\n\nCreate with:\n\nopen(file) do io\n    Binning(\n        io::IO,\n        ref::Reference;\n        min_size::Integer=1,\n        min_seqs::Integer=1,\n        binsplit_separator::Union{AbstractString, Char, Nothing}=nothing,\n        disjoint::Bool=true,\n        recalls=DEFAULT_RECALLS,\n        precisions=DEFAULT_PRECISIONS\n)\n\nmin_size: Filter away bins with breadth lower than this\nmin_seqs: Filter away bins with fewer sequences that this\nbinsplit_separator: Split bins based on this separator (nothing means no binsplitting)\ndisjoint: Throw an error if the same sequence is seen in multiple bins\nrecalls and precision: The thresholds to benchmark with\n\n\n\n\n\n","category":"type"},{"location":"reference/#VambBenchmarks.Clade","page":"Reference","title":"VambBenchmarks.Clade","text":"Clade{Genome}(name::AbstractString, child::Union{Clade{Genome}, Genome})\n\nA Clade represents any clade above Genome. Every Genome is expected to belong to the same number of clades, e.g. there may be exactly 7 levels of clades above every Genome. Clades always have at least one child (which is either a Genome or a Clade one rank lower), and a parent, unless it's the unique top clade from which all other clades and genomes descend from. The rank of a Genome is 0, clades that contain genomes have rank 1, and clades containing rank-1 clades have rank 2 etc. By default, zero-indexed ranks correspond to OTU, species, genus, family, order, class, phylum and domain.\n\nExamples\n\njulia> top_clade(ref)\nGenus \"F\", 3 genomes\n├─ Species \"D\", 2 genomes\n│  ├─ Genome(gA)\n│  └─ Genome(gB)\n└─ Species \"E\", 1 genome\n   └─ Genome(gC)\n\njulia> top_clade(ref).children\n2-element Vector{Clade{Genome}}:\n Species \"D\", 2 genomes\n Species \"E\", 1 genome\n\n\n\n\n\n","category":"type"},{"location":"reference/#VambBenchmarks.FlagSet","page":"Reference","title":"VambBenchmarks.FlagSet","text":"FlagSet <: AbstractSet{Flag}\n\nFlags are compact sets of Flag associated to a Genome. You can construct them from an iterable of Flag, e.g. a 1-element tuple. FlagSet support most set operations efficiently.\n\nExamples\n\njulia> flags = FlagSet((Flags.organism, Flags.virus));\n\njulia> Flags.virus in flags\ntrue\n\njulia> isdisjoint(flags, FlagSet((Flags.organism,)))\nfalse\n\nSee also: Flag, Genome\n\n\n\n\n\n","category":"type"},{"location":"reference/#VambBenchmarks.Flags.Flag","page":"Reference","title":"VambBenchmarks.Flags.Flag","text":"Flag\n\nA flag is a boolean associated to a Genome, stored in a Flags object. A flag may be e.g. Flag.organism, signaling that the genome is known to be an organism.\n\nSee also: FlagSet, Genome\n\n\n\n\n\n","category":"type"},{"location":"reference/#VambBenchmarks.Genome","page":"Reference","title":"VambBenchmarks.Genome","text":"Genome(name::AbstractString [flags::FlagSet])\n\nGenomes represent individual target genomes (organisms, plasmids, viruses etc), analogous to lowest-level clade that can be reconstructed. Conceptually, Genomes contain one or more Sources, and to a single parent Clade. They are identified uniquely among genomes by their name.\n\nA genome have a genome size, which is the sum of the length of all its sources. We consider this to be the true size of the biological genome (assuming its full sequence is contained in its sources), as well as an assembly size, which represent the sum of the assembly sizes of each source.\n\nSee also: Clade, Source, mrca\n\nExamples\n\njulia> gA, gB, gC = collect(ref.genomes);\n\njulia> flags(gA)\nFlagSet with 1 element:\n  VambBenchmarks.Flags.organism\n\njulia> mrca(gA, gB)\nSpecies \"D\", 2 genomes\n├─ Genome(gA)\n└─ Genome(gB)\n\n\n\n\n\n","category":"type"},{"location":"reference/#VambBenchmarks.Reference","page":"Reference","title":"VambBenchmarks.Reference","text":"Reference(::IO; [min_seq_length=1])\n\nA Reference contains the ground truth to benchmark against. Conceptually, it consists of the following parts:\n\nA list of genomes, each with sources\nThe full taxonomic tree, as lists of clades\nA list of sequences, each with a list of (source, span) to where it maps.\n\nNormally, the types FlagSet Genome, Source, Clade and Sequence do not need to be constructed manually, but are constructed when the Reference is loaded from a JSON file.\n\nExamples\n\njulia> ref = open(path_to_ref_file) do io\n           Reference(io; min_seq_length=1)\n       end;\n\njulia> ref isa Reference\ntrue\n\njulia> length(genomes(ref))\n3\n\njulia> nseqs(ref)\n11\n\njulia> first(ref.genomes) isa Genome\ntrue\n\nSee also: subset, Genome, Clade\n\n\n\n\n\n","category":"type"},{"location":"reference/#VambBenchmarks.Sequence","page":"Reference","title":"VambBenchmarks.Sequence","text":"Sequence(name::AbstractString, length::Integer)\n\nType that represents a binnable sequence. Sequences do not contain other information than their name and their length, and are identified by their name.\n\nExamples\n\njulia> Sequence(\"abc\", 5)\nSequence(\"abc\", 5)\n\njulia> Sequence(\"abc\", 5) == Sequence(\"abc\", 9)\ntrue\n\njulia> Sequence(\"abc\", 0)\nERROR: ArgumentError: Cannot instantiate an empty sequence\n\n\n\n\n\n","category":"type"},{"location":"reference/#VambBenchmarks.Source","page":"Reference","title":"VambBenchmarks.Source","text":"Source{Genome}(g::Genome, name::AbstractString, length::Integer)\n\nSources are the \"ground truth\" sequences that the binning attempts to recreate. For example, the assembled contigs of the reference genome (typically full, closed circular contigs) as found in NCBI or elsewhere are each Sources. Many Genomes only contain a single Source namely its full assembled genome. Each Source has a single parent Genome, and a unique name which identifies it.\n\nSources have zero or more mapping Sequences, that each map to the Source at a given span given by a UnitRange{Int}.\n\nSources have an assembly size, which is the number of base pairs where any sequence map to.\n\n\n\n\n\n","category":"type"},{"location":"reference/#VambBenchmarks.flags-Tuple{Genome}","page":"Reference","title":"VambBenchmarks.flags","text":"flags(g::Genome)::FlagSet\n\nReturns the Flags of the Genome as a FlagSet.\n\nSee also: Flag, FlagSet\n\nExample\n\njulia> flags(genome)\nFlagSet with 1 element:\n  VambBenchmarks.Flags.organism\n\n\n\n\n\n","category":"method"},{"location":"reference/#VambBenchmarks.gold_standard-Tuple{Reference}","page":"Reference","title":"VambBenchmarks.gold_standard","text":"gold_standard(\n    ref::Reference;\n    disjoint=true,\n    recalls=DEFAULT_RECALLS,\n    precisions=DEFAULT_PRECISIONS\n)::Binning\n\nCreate the optimal Binning object given an assembly. If disjoint, assign each sequence to only a single genome.\n\nExtended help\n\nCurrently, the disjoint option uses a simple greedy algorithm to assign sequences to genomes.\n\n\n\n\n\n","category":"method"},{"location":"reference/#VambBenchmarks.intersecting-Tuple{Bin}","page":"Reference","title":"VambBenchmarks.intersecting","text":"intersecting([Genome, Clade]=Genome, x::Bin)\n\nGet an iterator of the Genomes or Clades that bin x intersects with. intersecting(::Bin) defaults to genomes.\n\nExample\n\njulia> collect(intersecting(bin))\n1-element Vector{Genome}:\n Genome(gA)\n\njulia> sort!(collect(intersecting(Clade, bin)); by=i -> i.name)\n2-element Vector{Clade{Genome}}:\n Species \"D\", 2 genomes\n Genus \"F\", 3 genomes\n\n\n\n\n\n","category":"method"},{"location":"reference/#VambBenchmarks.is_organism-Tuple{Genome}","page":"Reference","title":"VambBenchmarks.is_organism","text":"is_organism(g::Genome)::Bool\n\nCheck if g is known to be an organism.\n\nExample\n\njulia> is_organism(genome)\ntrue\n\n\n\n\n\n","category":"method"},{"location":"reference/#VambBenchmarks.is_plasmid-Tuple{Genome}","page":"Reference","title":"VambBenchmarks.is_plasmid","text":"is_plasmid(g::Genome)::Bool\n\nCheck if g is known to be a plasmid.\n\nExample\n\njulia> is_plasmid(genome)\nfalse\n\n\n\n\n\n","category":"method"},{"location":"reference/#VambBenchmarks.is_virus-Tuple{Genome}","page":"Reference","title":"VambBenchmarks.is_virus","text":"is_virus(g::Genome)::Bool\n\nCheck if g is known to be a virus.\n\nExample\n\njulia> is_virus(genome)\nfalse\n\n\n\n\n\n","category":"method"},{"location":"reference/#VambBenchmarks.mrca-Tuple{Union{Clade{Genome}, Genome}, Union{Clade{Genome}, Genome}}","page":"Reference","title":"VambBenchmarks.mrca","text":"mrca(a::Node, b::Node)::Node\n\nCompute the most recent common ancestor (MRCA) of a and b.\n\n\n\n\n\n","category":"method"},{"location":"reference/#VambBenchmarks.n_recovered-Tuple{Binning, Real, Real}","page":"Reference","title":"VambBenchmarks.n_recovered","text":"n_recovered(::Binning, recall, precision; level=0, assembly=false)::Integer\n\nReturn the number of genomes or clades reconstructed in the Binning at the given recall and precision levels. If assembly is set, return the number of assemblies reconstructed instead. The argument level sets the taxonomic rank: 0 for Genome (or assemblies).\n\nExamples\n\njulia> n_recovered(binning, 0.4, 0.71)\n1\n\njulia> n_recovered(binning, 0.4, 0.71; assembly=true)\n2\n\njulia> n_recovered(binning, 0.4, 0.71; assembly=true, level=2)\n1\n\n\n\n\n\n","category":"method"},{"location":"reference/#VambBenchmarks.print_matrix-Tuple{Binning}","page":"Reference","title":"VambBenchmarks.print_matrix","text":"print_matrix(::Binning; level=0, assembly=true)\n\nPrint the number of reconstructed assemblies or genomes at the given taxonomic level (rank). Level 0 corresponds to genomes, level 1 to species, etc. If assembly, print the number of reconstructed assemblies, else print the level of reconstructed genomes.\n\nSee also: Binning\n\n\n\n\n\n","category":"method"},{"location":"reference/#VambBenchmarks.recall_precision-Tuple{Genome, Bin}","page":"Reference","title":"VambBenchmarks.recall_precision","text":"recall_precision(x::Union{Genome, Clade}, bin::Bin; assembly::Bool=true)\n\nGet the recall, precision as a 2-tuple of Float64 for the given genome/bin pair. See the docstring for Bin for how this is computed.\n\nSee also: Bin, Binning\n\nExamples\n\njulia> bingenome = only(intersecting(bin));\n\njulia> recall_precision(bingenome, bin)\n(recall = 0.45454545454545453, precision = 1.0)\n\njulia> recall_precision(bingenome, bin; assembly=false)\n(recall = 0.4, precision = 1.0)\n\njulia> recall_precision(bingenome.parent, bin; assembly=false)\n(recall = 0.4, precision = 1.0)\n\n\n\n\n\n","category":"method"},{"location":"reference/#VambBenchmarks.subset!-Tuple{Reference}","page":"Reference","title":"VambBenchmarks.subset!","text":"subset!(\n        ref::Reference;\n        sequences::Function=Returns(true),\n        genomes::Function=Returns(true)\n)::Reference\n\nMutate ref in place, removing genomes and sequences. Keep only sequences S where sequences(S) returns true and genomes G for which genomes(G) returns true.\n\nSee also: subset, Reference\n\nExamples\n\njulia> ref\nReference\n  Genomes:    3\n  Sequences:  11\n  Ranks:      3\n  Seq length: 10\n  Assembled:  66.8 %\n\njulia> subset(ref; genomes=g -> Flags.organism in flags(g))\nReference\n  Genomes:    2\n  Sequences:  11\n  Ranks:      3\n  Seq length: 10\n  Assembled:  91.3 %\n\njulia> VambBenchmarks.subset(ref; sequences=s -> length(s) ≥ 25)\nReference\n  Genomes:    3\n  Sequences:  9\n  Ranks:      3\n  Seq length: 25\n  Assembled:  61.1 %\n\n\n\n\n\n","category":"method"},{"location":"reference/#VambBenchmarks.subset-Tuple{Reference}","page":"Reference","title":"VambBenchmarks.subset","text":"subset(ref::Reference; kwargs...)\n\nNon-mutating copying version of subset.\n\nSee also: subset!\n\n\n\n\n\n","category":"method"},{"location":"walkthrough/#Example-walkthrough","page":"Walkthrough","title":"Example walkthrough","text":"","category":"section"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"Note: This is run 2023-04-05 with commit AFTER bc99339","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"danger: Danger\nThe documentation on this page uses a large dataset, and which is stored hosted online. Hence, the code in this documentation is not tested, and may become outdated.","category":"page"},{"location":"walkthrough/#Loading-the-reference","page":"Walkthrough","title":"Loading the reference","text":"","category":"section"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"First, let's load the reference:","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"julia> using VambBenchmarks\n\njulia> ref_path = \"/home/jakni/Downloads/vambdata/newref/ref_ptracker_megahit_Skin.json\";\n\njulia> bins_path = \"/home/jakni/Downloads/vambdata/newbins/megahit_skin.tsv\";\n\njulia> ref = open(i -> Reference(i), ref_path)\nReference\n  Genomes:    2341\n  Sequences:  959970\n  Ranks:      9\n  Seq length: 200\n  Assembled:  25.0 %","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"This gives us a few statistics about the reference:","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"Number of genomes\nNumber of sequences\nNumber of taxonomic ranks (strain, species, genus...)\nLength of shortest sequence\nTotal length of genomes that are assembled","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"The genomes here contain both plasmids and organisms, and the sequence length of 200 bp is too short. Let's filter the reference using subset! to only retain organisms, and sequences of length 1500 or more:","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"julia> subset!(ref;\n           genomes = is_organism,\n           sequences = s -> length(s) >= 1500\n       )\nReference\n  Genomes:    1394\n  Sequences:  118267\n  Ranks:      9\n  Seq length: 1500\n  Assembled:  16.1 %","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"note: Note\nThe function subset! will mutate the reference, whereas the function subset will create a new independent reference. At the moment, the latter is much slower.","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"You can see we removed 7/8th of all sequences, but only 1/3rd of the total assembly length. We also removed almost half of all genomes (namely, all those that wasn't organisms).","category":"page"},{"location":"walkthrough/#Genomes","page":"Walkthrough","title":"Genomes","text":"","category":"section"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"We can explore the genomes contained in the reference with the genomes function, which returns an iterable of Genome (in this case, a Set):","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"julia> genomes(ref)\nSet{Genome} with 1394 elements:\n  Genome(OTU_97.34832.0)\n  Genome(OTU_97.33610.0)\n  Genome(OTU_97.39616.0)\n  Genome(OTU_97.23938.0)\n  Genome(OTU_97.24774.0)\n  Genome(OTU_97.36286.0)\n  Genome(OTU_97.6382.0)\n  Genome(OTU_97.33661.1)\n  Genome(OTU_97.37247.0)\n  Genome(OTU_97.1829.0)\n  Genome(OTU_97.6602.0)\n  Genome(OTU_97.44820.0)\n  Genome(OTU_97.24834.0)\n  Genome(OTU_97.44856.1)\n  Genome(OTU_97.3377.0)\n  Genome(OTU_97.19529.1)\n  Genome(OTU_97.32371.0)\n  Genome(OTU_97.45083.0)\n  Genome(OTU_97.15121.1)\n  ⋮","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"Let's look at a Genome in more detail:","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"julia> genome, genome2 = sort!(collect(genomes(ref)); by=i -> i.name);\n\njulia> genome\nGenome \"OTU_97.10046.0\"\n  Parent:        \"Staphylococcus capitis\"\n  Genome size:   2474232\n  Assembly size: 9239 (0.4 %)\n  Sources:       1\n  Flags:         1 (organism)","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"The flags can be extracted with the flags(genome) function - each genome contains zero or more flags:","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"julia> flags(genome)\nFlagSet with 1 element:\n  VambBenchmarks.Flags.organism","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"... in this case, this genome is an organism as opposed to a plasmid. You can see all possible flags with instances(Flags.Flag).","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"See also the helper functions is_organism, is_virus and is_plasmid","category":"page"},{"location":"walkthrough/#Sources","page":"Walkthrough","title":"Sources","text":"","category":"section"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"The genome has one source - let's look at that","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"julia> source = only(genome.sources)\nSource \"CP007601.1\"\ngenome:          Genome(OTU_97.10046.0)\n  Length:        2474232\n  Assembly size: 9239\n  Sequences:     2","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"A Source is one of the genomic sequences that genomes are composed of. In the reference database where our genome was gotten from, it was assembled into a single contig called \"CP001781.1\" of 2.47 Mbp - presumably the complete, circular genome. We can see that 2 sequences map to this source, covering about 9.2 kbp.","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"We can get the sequences mapping to this source:","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"julia> source.sequences\n2-element Vector{Tuple{Sequence, UnitRange{Int64}}}:\n (Sequence(\"S13C8728\", 2472), 2069906:2072374)\n (Sequence(\"S13C26910\", 7349), 2085051:2091820)","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"Where, e.g. the first entrance tells us that the sequence \"S13C8728\" with a length of 2472 maps to positions 2069906:2072374 (both inclusive).","category":"page"},{"location":"walkthrough/#Clades","page":"Walkthrough","title":"Clades","text":"","category":"section"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"Genomes are organised into a taxonomic hierarchy. We can find the immediate parent of a genome by accessing the field genome.parent. Let's look at another genome:","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"julia> genome2\nGenome \"OTU_97.10083.0\"\n  Parent:        \"Moraxella bovoculi\"\n  Genome size:   2220786\n  Assembly size: 0 (0.0 %)\n  Sources:       1\n  Flags:         1 (organism)\n\njulia> clade = genome2.parent\nSpecies \"Moraxella bovoculi\", 6 genomes\n├─ Genome(OTU_97.36230.1)\n├─ Genome(OTU_97.11114.0)\n├─ Genome(OTU_97.10083.0)\n├─ Genome(OTU_97.9999.0)\n├─ Genome(OTU_97.36230.0)\n└─ Genome(OTU_97.30020.0)","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"The parent is an instance of a Clade. Clades are at a specific rank: Rank 1 for species, 2 for genus, 3 for family, etc. Every clade has one or more children: These are the clades one rank lower. Conceptually, rank zero corresponds to Genomes (OTUs, for this reference dataset)","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"julia> clade.children\n6-element Vector{Genome}:\n Genome(OTU_97.36230.1)\n Genome(OTU_97.11114.0)\n Genome(OTU_97.10083.0)\n Genome(OTU_97.9999.0)\n Genome(OTU_97.36230.0)\n Genome(OTU_97.30020.0)","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"We can find the most recent common ancestor (MRCA) of genome and g2 like this:","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"julia> mrca(genome, genome2)\nDomain \"Bacteria\", 1394 genomes\n├─ Phylum \"Deinococcus-Thermus\", 12 genomes\n│  └─ Class \"Deinococci\", 12 genomes\n│     └─ Order \"Deinococcales\", 12 genomes\n│        ⋮\n│\n├─ Phylum \"Tenericutes\", 5 genomes\n│  └─ Class \"Mollicutes\", 5 genomes\n│     └─ Order \"Mycoplasmatales\", 5 genomes\n│        ⋮\n│\n├─ Phylum \"Fusobacteria\", 21 genomes\n│  └─ Class \"Fusobacteriia\", 21 genomes\n│     └─ Order \"Fusobacteriales\", 21 genomes\n│        ⋮\n│\n├─ Phylum \"Bacteroidetes\", 85 genomes\n│  ├─ Class \"Sphingobacteriia\", 1 genome\n│  │  └─ Order \"Sphingobacteriales\", 1 genome\n│  │     ⋮\n│  │\n│  ├─ Class \"Bacteroidia\", 47 genomes\n│  │  └─ Order \"Bacteroidales\", 47 genomes\n│  │     ⋮\n│  │\n⋮","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"They are very distantly related, so the domain \"Bacteria\", one of the highest ranked Clades, are their most recent common ancestor.","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"The top clade can be found with the top_clade(ref) function. In our case, the dataset only consists of bacteria and plasmids, and we filtered away all plasmids, so the universal ancestor only has one child, namely the domain Bacteria.","category":"page"},{"location":"walkthrough/#Binnings","page":"Walkthrough","title":"Binnings","text":"","category":"section"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"A Binning is a set of bins benchmarked against a reference. We can load a set of Vamb bins and turn it into a Binning object like this:","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"julia> binning = open(bins_path) do io\n           Binning(io, ref; binsplit_separator='C')\n       end\nBinning\n  Reference\n    Genomes:    1394\n    Sequences:  118267\n    Ranks:      9\n    Seq length: 1500\n    Assembled:  16.1 %\n  Bins:        14535\n  NC genomes:  10\n  Precisions: [0.6, 0.7, 0.8, 0.9, 0.95, 0.99]\n  Recalls:    [0.6, 0.7, 0.8, 0.9, 0.95, 0.99]\n  Recoverable genomes: [193, 142, 98, 51, 21, 0]\n  Reconstruction (assemblies):\n    P\\R   0.6  0.7  0.8  0.9 0.95 0.99\n    0.6   160  135   87   39   14    4\n    0.7   143  122   80   36   13    3\n    0.8   118  101   71   33   13    3\n    0.9    79   71   55   29   12    3\n    0.95   70   63   50   26   12    3\n    0.99   57   51   39   23   10    2","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"A wealth of information is readily available:","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"binning.ref gives the underlying Reference\nbinning.recalls and binning.precision gives the recall/precision thresholds used in benchmarking\nYou can get the number of genomes that are assembled at the various recall levels with recoverable_genomes: This sets an upper limit of how many genomes can be reconstructed given a hypothetical perfect binning","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"julia> println(binning.recoverable_genomes)\n[193, 142, 98, 51, 21, 0]","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"The function print_matrix will display the number of recovered genomes/assemblies. It takes two optional keyword: level, the taxonomic rank (defaults to 0, meaning strain level), and assembly which defaults to true. If set to false, benchmark number of recovered genomes, not number of recovered assemblies.","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"julia> print_matrix(binning; level=1, assembly=false)\nP\\R   0.6  0.7  0.8  0.9 0.95 0.99\n0.6    62   48   31   10    1    0\n0.7    59   46   30   10    1    0\n0.8    59   46   30   10    1    0\n0.9    59   46   30   10    1    0\n0.95   56   44   30   10    1    0\n0.99   46   37   25    9    1    0","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"You can also get the number of genomes or assemblies reconstructed at a given precision/recall level directly with n_recovered:","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"julia> n_recovered(binning, 0.75, 0.9; assembly=true)\n55\n\njulia> n_recovered(binning, 0.66, 0.91; level=1)\n44","category":"page"},{"location":"walkthrough/#Bins","page":"Walkthrough","title":"Bins","text":"","category":"section"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"The Binning object obviously contains our bins. Let's pick a particularly good bin:","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"julia> bin = binning.bins[261]\nBin \"S13Cvae_479\"\n  Sequences: 47\n  Breadth:   2482216\n  Intersecting 1 genome","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"The \"breadth\" here is the sum of the length of its sequences. bin.sequences gets an iterable of its sequences:","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"julia> collect(bin.sequences)\n47-element Vector{Sequence}:\n Sequence(\"S13C2748\", 73726)\n Sequence(\"S13C13383\", 18077)\n Sequence(\"S13C59849\", 95631)\n Sequence(\"S13C29443\", 13209)\n Sequence(\"S13C56904\", 44356)\n Sequence(\"S13C14948\", 9756)\n Sequence(\"S13C37848\", 19070)\n Sequence(\"S13C2782\", 22054)\n Sequence(\"S13C40546\", 46579)\n Sequence(\"S13C12978\", 15372)\n ⋮\n Sequence(\"S13C33782\", 34963)\n Sequence(\"S13C27111\", 21961)\n Sequence(\"S13C10862\", 55887)\n Sequence(\"S13C54667\", 69970)\n Sequence(\"S13C19757\", 199299)\n Sequence(\"S13C64535\", 73792)\n Sequence(\"S13C2686\", 20376)\n Sequence(\"S13C62986\", 19619)\n Sequence(\"S13C14459\", 39940)","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"The \"Intersecting 1 genome\" means that the sequences map to 1 genome. We can get that with the function intersecting, then get the precision/recall with recall_precision:","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"julia> intersecting_genome = only(intersecting(bin))\nGenome \"OTU_97.9674.1\"\n  Parent:        \"Corynebacterium humireducens\"\n  Genome size:   2714910\n  Assembly size: 2575299 (94.9 %)\n  Sources:       1\n  Flags:         1 (organism)\n\njulia> recall_precision(only(intersecting(bin)), bin)\n(recall = 0.9614992278566489, precision = 1.0)","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"Or do the same for a higher clade - let's say a genus:","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"julia> genus = only(Iterators.filter(i -> i.rank == 2, intersecting(Clade, bin)));\n\njulia> recall_precision(genus, bin)\n(recall = 0.9614992278566489, precision = 1.0)","category":"page"},{"location":"#VambBenchmarks","page":"Home","title":"VambBenchmarks","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"VambBenchmarks.jl is a package for efficient benchmarking and interactive exploration of a set of bins against a reference.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Install Julia - preferably using juliaup: https://github.com/JuliaLang/juliaup\nLaunch Julia: julia\nPress ] to enter package mode. You can exit package mode with backspace.\nIn package mode, type add https://github.com/jakobnissen/VambBenchmarks.jl to download and install the benchmarking software","category":"page"},{"location":"#Quickstart","page":"Home","title":"Quickstart","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"using VambBenchmarks\nref =  open(i -> Reference(i), \"files/ref.json\")\nbins = open(i -> Binning(i, ref), \"files/clusters.tsv\")\nprint_matrix(bins)","category":"page"},{"location":"#Concepts","page":"Home","title":"Concepts","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A Sequence is a sequence (e.g. contig) clustered by the binner\nA Genome is a target genome that should be reconstructed by the binner. It can be a virus, organism, plasmid etc. Every Genome have several Sources, and one parent Clade.\nA Flag marks the certaincy about a boolean attribute of a genome, like \"is this a virus?\".\nSources are the sequences that Genomes are composed of. These are typically the reference genome sequences originally obtained by assembly of a purified genome (e.g. clonal colony). Sequences map to zero or more Sources at particular spans, i.e. locations.\nA Clade contain one or more Genomes or Clades. Clades containing genomes are rank 1, and clades containing rank N clades are rank N+1 clades. All genomes descend from a chain of exactly N ranks of clades, where N > 0.\nA Bin is a set of Sequences created by the binner. Every bin is benchmarked against all genomes and clades in the reference.\nA Reference is composed of:\nThe genomes, a set of Genomes, each with a set of Sources and Flags\nThe taxmaps, a full set of Clades that encompasses every Genome at N ranks (where N > 0)\nThe sequences, a list of Sequences, each with zero or more mappings to Sources.\nA Binning is a set of Bins benchmarked against a Reference","category":"page"},{"location":"","page":"Home","title":"Home","text":"See the Reference in the left sidebar.","category":"page"}]
}
