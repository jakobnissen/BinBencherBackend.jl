var documenterSearchIndex = {"docs":
[{"location":"reference/#Reference","page":"Reference","title":"Reference","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"Modules = [BinBencherBackend]\nOrder   = [:type, :function]","category":"page"},{"location":"reference/#BinBencherBackend.Bin","page":"Reference","title":"BinBencherBackend.Bin","text":"Bin(name::AbstractString, ref::Reference, sequences)\n\nBins each represent a bin created by the binner. Conceptually, they are simply a set of Sequence with a name attached. Practically, every Bin is benchmarked against all Genomes and Clades of a given Reference, so each Bin stores data about its intersection with every genome/clade, e.g. its purity and recall.\n\nLike Sources, Bins also have an assembly size for a given Genome. This is the number of base pairs in the genomes covered by any sequence in the Bin, which is always a subset of the genome's assembly size.\n\nBenchmark statistics for a Bin/Genome can be done with either assemblies or genomes as the ground truth.\n\nTrue positives (TP) are defined as the sum of assembly sizes over all sources in the genome\nFalse positives (FP) are the sum of length of sequences in the bin not mapping to the genome\nFalse negatives (FN) is either the genome assembly size or genome size minus TP.\n\nFor Bin/Clade pairs B/C, recall is the maximal recall of B/Ch for all children Ch of C. Precision is the sum of lengths of sequences mapping to any child of the clade divided by the sum of lengths of all sequences in the bin.\n\nSee also: Binning, Genome, Clade\n\nExamples\n\njulia> bin = first(binning.bins)\nBin \"C1\"\n  Sequences: 2\n  Breadth:   65\n  Intersecting 1 genome\n\njulia> first(bin.sequences)\nSequence(\"s1\", 25)\n\njulia> f1(first(ref.genomes), bin)\n0.5714285714285715\n\n\n\n\n\n","category":"type"},{"location":"reference/#BinBencherBackend.Binning","page":"Reference","title":"BinBencherBackend.Binning","text":"Binning(::Union{IO, AbstractString}, ::Reference; kwargs...)\n\nA Binning represents a set of Bins benchmarked against a Reference. Binnings can be created given a set of Bins and a Reference, where the bins may potentially be loaded from a .tsv file. The fields recovered_asms and recovered_genomes are used for benchmarking, these are normally output using the print_matrix function.\n\nA Binning is loaded from a tsv file, which is specified either as an IO, or its path as an AbstractString. If the path ends with .gz, automatically gzip decompress when reading the file.\n\nSee also: print_matrix, Bin, Reference\n\nExamples\n\njulia> bins = Binning(path_to_bins_file, ref);\n\n\njulia> bins isa Binning\ntrue\n\njulia> BinBencherBackend.n_nc(binning)\n0\n\nExtended help\n\nCreate with:\n\nopen(file) do io\n    Binning(\n        io::Union{IO, AbstractString},\n        ref::Reference;\n        min_size::Integer=1,\n        min_seqs::Integer=1,\n        binsplit_separator::Union{AbstractString, Char, Nothing}=nothing,\n        disjoint::Bool=true,\n        recalls=DEFAULT_RECALLS,\n        precisions=DEFAULT_PRECISIONS,\n        filter_genomes=Returns(true)\n)\n\nmin_size: Filter away bins with breadth lower than this\nmin_seqs: Filter away bins with fewer sequences that this\nbinsplit_separator: Split bins based on this separator (nothing means no binsplitting)\ndisjoint: Throw an error if the same sequence is seen in multiple bins\nrecalls and precision: The thresholds to benchmark with\nfilter_genomes: A function f(genome)::Bool. Genomes for which it returns  false are ignored in benchmarking.\n\n\n\n\n\n","category":"type"},{"location":"reference/#BinBencherBackend.Clade","page":"Reference","title":"BinBencherBackend.Clade","text":"Clade{Genome}(name::AbstractString, child::Union{Clade{Genome}, Genome})\n\nA Clade represents any clade above Genome. Every Genome is expected to belong to the same number of clades, e.g. there may be exactly 7 levels of clades above every Genome. Clades always have at least one child (which is either a Genome or a Clade one rank lower), and a parent, unless it's the unique top clade from which all other clades and genomes descend from. The rank of a Genome is 0, clades that contain genomes have rank 1, and clades containing rank-1 clades have rank 2 etc. By default, zero-indexed ranks correspond to OTU, species, genus, family, order, class, phylum and domain.\n\nExamples\n\njulia> top_clade(ref)\nGenus \"F\", 3 genomes\n├─ Species \"D\", 2 genomes\n│  ├─ Genome(gA)\n│  └─ Genome(gB)\n└─ Species \"E\", 1 genome\n   └─ Genome(gC)\n\njulia> top_clade(ref).children\n2-element Vector{Clade{Genome}}:\n Species \"D\", 2 genomes\n Species \"E\", 1 genome\n\n\n\n\n\n","category":"type"},{"location":"reference/#BinBencherBackend.FlagSet","page":"Reference","title":"BinBencherBackend.FlagSet","text":"FlagSet <: AbstractSet{Flag}\n\nFlags are compact sets of Flag associated to a Genome. You can construct them from an iterable of Flag, e.g. a 1-element tuple. FlagSet support most set operations efficiently.\n\nSee also: Flag, Genome\n\nExamples\n\njulia> flags = FlagSet((Flags.organism, Flags.virus));\n\n\njulia> Flags.virus in flags\ntrue\n\njulia> isdisjoint(flags, FlagSet((Flags.organism,)))\nfalse\n\n\n\n\n\n","category":"type"},{"location":"reference/#BinBencherBackend.Flags.Flag","page":"Reference","title":"BinBencherBackend.Flags.Flag","text":"Flag\n\nA flag is a boolean associated to a Genome, stored in a Flags object. A flag may be e.g. Flag.organism, signaling that the genome is known to be an organism.\n\nSee also: FlagSet, Genome\n\nExamples\n\njulia> tryparse(Flag, \"organism\") == Flags.organism\ntrue\n\njulia> tryparse(Flag, \"Canada\") === nothing\ntrue\n\n\n\n\n\n","category":"type"},{"location":"reference/#BinBencherBackend.Genome","page":"Reference","title":"BinBencherBackend.Genome","text":"Genome(name::AbstractString [flags::FlagSet])\n\nGenomes represent individual target genomes (organisms, plasmids, viruses etc), and are conceptually the lowest-level clade that can be reconstructed. Genomes contain one or more Sources, and belong to a single parent Clade. They are identified uniquely among genomes by their name.\n\nA genome have a genome size, which is the sum of the length of all its sources. We consider this to be the true size of the biological genome (assuming its full sequence is contained in its sources), as well as an assembly size, which represent the sum of the assembly sizes of each source.\n\nSee also: Clade, Source, mrca\n\nExamples\n\njulia> gA, gB, gC = collect(ref.genomes);\n\n\njulia> flags(gA)\nFlagSet with 1 element:\n  BinBencherBackend.Flags.organism\n\njulia> mrca(gA, gB)\nSpecies \"D\", 2 genomes\n├─ Genome(gA)\n└─ Genome(gB)\n\n\n\n\n\n","category":"type"},{"location":"reference/#BinBencherBackend.Reference","page":"Reference","title":"BinBencherBackend.Reference","text":"Reference(::Union{IO, AbstractString})\n\nA Reference contains the ground truth to benchmark against. Conceptually, it consists of the following parts:\n\nA list of genomes, each with sources\nThe full taxonomic tree, as lists of clades\nA list of sequences, each with a list of (source, span) to where it maps.\n\nNormally, the types FlagSet Genome, Source, Clade and Sequence do not need to be constructed manually, but are constructed when the Reference is loaded from a JSON file.\n\nA Reference is loaded from a JSON file, which is specified either as an IO, or its path as an AbstractString. If the path ends with .gz, automatically gzip decompress when reading the file.\n\nExamples\n\njulia> ref = Reference(path_to_ref_file);\n\n\njulia> ref isa Reference\ntrue\n\njulia> length(genomes(ref))\n3\n\njulia> nseqs(ref)\n11\n\njulia> first(ref.genomes) isa Genome\ntrue\n\nSee also: subset, Genome, Clade\n\n\n\n\n\n","category":"type"},{"location":"reference/#BinBencherBackend.Sequence","page":"Reference","title":"BinBencherBackend.Sequence","text":"Sequence(name::AbstractString, length::Integer)\n\nType that represents a binnable sequence. Sequences do not contain other information than their name and their length, and are identified by their name.\n\nExamples\n\n```jldoctest julia> Sequence(\"abc\", 5) Sequence(\"abc\", 5)\n\njulia> Sequence(\"abc\", 5) == Sequence(\"abc\", 9) true\n\njulia> Sequence(\"abc\", 0) ERROR: ArgumentError: Cannot instantiate an empty sequence [...]\n\n\n\n\n\n","category":"type"},{"location":"reference/#BinBencherBackend.Source","page":"Reference","title":"BinBencherBackend.Source","text":"Source{Genome}(g::Genome, name::AbstractString, length::Integer)\n\nSources are the \"ground truth\" sequences that the binning attempts to recreate. For example, the assembled contigs of the reference genome (typically full, closed circular contigs) as found in NCBI or elsewhere are each Sources. Many Genomes only contain a single Source namely its full assembled genome. Each Source has a single parent Genome, and a unique name which identifies it.\n\nSources have zero or more mapping Sequences, that each map to the Source at a given span given by a 2-tuple Tuple{Int, Int}.\n\nSources have an assembly size, which is the number of base pairs where any sequence map to.\n\n\n\n\n\n","category":"type"},{"location":"reference/#BinBencherBackend.assembly_size!-Tuple{Function, Vector{Tuple{Int64, Int64}}, Vector, Int64}","page":"Reference","title":"BinBencherBackend.assembly_size!","text":"Compute -> (breadth, totalbp), where breadth is the number of positions in v covered at least once, and totalbp the sum of the lengths of the sequences. v must be a Vector such that all(by(i) isa Tuple{Integer, Integer} for i in v). The scratch input is mutated.\n\n\n\n\n\n","category":"method"},{"location":"reference/#BinBencherBackend.flags-Tuple{Genome}","page":"Reference","title":"BinBencherBackend.flags","text":"flags(g::Genome)::FlagSet\n\nReturns the Flags of the Genome as a FlagSet.\n\nSee also: Flag, FlagSet\n\nExample\n\njulia> flags(genome)\nFlagSet with 1 element:\n  BinBencherBackend.Flags.organism\n\n\n\n\n\n","category":"method"},{"location":"reference/#BinBencherBackend.gold_standard-Tuple{Reference}","page":"Reference","title":"BinBencherBackend.gold_standard","text":"gold_standard(\n    ref::Reference\n    [sequences, a Binning or an iterable of Sequence];\n    disjoint=true,\n    recalls=DEFAULT_RECALLS,\n    precisions=DEFAULT_PRECISIONS\n)::Binning\n\nCreate the optimal Binning object given a Reference, by the optimal binning of the Sequences in sequences. If disjoint, assign each sequence to only a single genome.\n\nIf sequences is not passed, use all sequences in ref. If a Binning is passed, use all sequences in any of its bins. Else, pass an iterable of Sequence.\n\nExtended help\n\nCurrently, the disjoint option uses a simple greedy algorithm to assign sequences to genomes.\n\n\n\n\n\n","category":"method"},{"location":"reference/#BinBencherBackend.intersecting-Tuple{Bin}","page":"Reference","title":"BinBencherBackend.intersecting","text":"intersecting([Genome, Clade]=Genome, x::Bin)\n\nGet an iterator of the Genomes or Clades that bin x intersects with. intersecting(::Bin) defaults to genomes.\n\nExample\n\njulia> collect(intersecting(bin))\n1-element Vector{Genome}:\n Genome(gA)\n\njulia> sort!(collect(intersecting(Clade, bin)); by=i -> i.name)\n2-element Vector{Clade{Genome}}:\n Species \"D\", 2 genomes\n Genus \"F\", 3 genomes\n\n\n\n\n\n","category":"method"},{"location":"reference/#BinBencherBackend.is_organism-Tuple{Genome}","page":"Reference","title":"BinBencherBackend.is_organism","text":"is_organism(g::Genome)::Bool\n\nCheck if g is known to be an organism.\n\nExample\n\njulia> is_organism(genome)\ntrue\n\n\n\n\n\n","category":"method"},{"location":"reference/#BinBencherBackend.is_plasmid-Tuple{Genome}","page":"Reference","title":"BinBencherBackend.is_plasmid","text":"is_plasmid(g::Genome)::Bool\n\nCheck if g is known to be a plasmid.\n\nExample\n\njulia> is_plasmid(genome)\nfalse\n\n\n\n\n\n","category":"method"},{"location":"reference/#BinBencherBackend.is_virus-Tuple{Genome}","page":"Reference","title":"BinBencherBackend.is_virus","text":"is_virus(g::Genome)::Bool\n\nCheck if g is known to be a virus.\n\nExample\n\njulia> is_virus(genome)\nfalse\n\n\n\n\n\n","category":"method"},{"location":"reference/#BinBencherBackend.mrca-Tuple{Union{Clade{Genome}, Genome}, Union{Clade{Genome}, Genome}}","page":"Reference","title":"BinBencherBackend.mrca","text":"mrca(a::Node, b::Node)::Node\n\nCompute the most recent common ancestor (MRCA) of a and b.\n\n\n\n\n\n","category":"method"},{"location":"reference/#BinBencherBackend.n_passing_bins-Tuple{Binning, Real, Real}","page":"Reference","title":"BinBencherBackend.n_passing_bins","text":"n_passing_bins(::Binning, recall, precision; level=0, assembly::Bool=false)::Integer\n\nReturn the number of bins which correspond to any genome or clade at the given recall and precision levels. If assembly is set, a recall of 1.0 means a bin corresponds to a whole assembly, else it corresponds to a whole genome. The argument level sets the taxonomic rank: 0 for Genome (or assemblies).\n\nExamples\n\njulia> n_passing_bins(binning, 0.4, 0.71)\n1\n\njulia> n_passing_bins(binning, 0.65, 0.71)\n0\n\n\n\n\n\n","category":"method"},{"location":"reference/#BinBencherBackend.n_recovered-Tuple{Binning, Real, Real}","page":"Reference","title":"BinBencherBackend.n_recovered","text":"n_recovered(::Binning, recall, precision; level=0, assembly=false)::Integer\n\nReturn the number of genomes or clades reconstructed in the Binning at the given recall and precision levels. If assembly is set, return the number of assemblies reconstructed instead. The argument level sets the taxonomic rank: 0 for Genome (or assemblies).\n\nExamples\n\njulia> n_recovered(binning, 0.4, 0.71)\n1\n\njulia> n_recovered(binning, 0.4, 0.71; assembly=true)\n2\n\njulia> n_recovered(binning, 0.4, 0.71; assembly=true, level=2)\n1\n\n\n\n\n\n","category":"method"},{"location":"reference/#BinBencherBackend.passes_f1-Tuple{Bin, Real}","page":"Reference","title":"BinBencherBackend.passes_f1","text":"passes_f1(bin::Bin, threshold::Real; assembly::Bool=false)::Bool\n\nComputes if bin has an F1 score equal to, or higher than threshold for any genome.\n\nExamples\n\njulia> obs_f1 = f1(only(intersecting(bin)), bin)\n0.5714285714285715\n\njulia> passes_f1(bin, obs_f1)\ntrue\n\njulia> passes_f1(bin, obs_f1 + 0.001)\nfalse\n\n\n\n\n\n","category":"method"},{"location":"reference/#BinBencherBackend.passes_recall_precision-Tuple{Bin, Real, Real}","page":"Reference","title":"BinBencherBackend.passes_recall_precision","text":"passes_recall_precision(bin::Bin, recall::Real, precision::Real; assembly::Bool=false)::Bool\n\nComputes if bin intersects with any Genome with at least the given recall and precision thresholds.\n\nExamples\n\njulia> (r, p) = recall_precision(only(intersecting(bin)), bin)\n(recall = 0.4, precision = 1.0)\n\njulia> passes_recall_precision(bin, 0.40, 1.0)\ntrue\n\njulia> passes_recall_precision(bin, 0.41, 1.0)\nfalse\n\n\n\n\n\n","category":"method"},{"location":"reference/#BinBencherBackend.print_matrix-Tuple{Binning}","page":"Reference","title":"BinBencherBackend.print_matrix","text":"print_matrix(::Binning; level=0, assembly=false)\n\nPrint the number of reconstructed assemblies or genomes at the given taxonomic level (rank). Level 0 corresponds to genomes, level 1 to species, etc. If assembly, print the number of reconstructed assemblies, else print the level of reconstructed genomes.\n\nSee also: Binning\n\n\n\n\n\n","category":"method"},{"location":"reference/#BinBencherBackend.recall_precision-Tuple{Genome, Bin}","page":"Reference","title":"BinBencherBackend.recall_precision","text":"recall_precision(x::Union{Genome, Clade}, bin::Bin; assembly::Bool=true)\n\nGet the recall, precision as a 2-tuple of Float64 for the given genome/bin pair. See the docstring for Bin for how this is computed.\n\nSee also: Bin, Binning\n\nExamples\n\njulia> bingenome = only(intersecting(bin));\n\n\njulia> recall_precision(bingenome, bin)\n(recall = 0.4, precision = 1.0)\n\njulia> recall_precision(bingenome, bin; assembly=false)\n(recall = 0.4, precision = 1.0)\n\njulia> recall_precision(bingenome.parent, bin; assembly=false)\n(recall = 0.4, precision = 1.0)\n\n\n\n\n\n","category":"method"},{"location":"reference/#BinBencherBackend.recursively_delete_child!-Tuple{T} where T<:Union{Clade{Genome}, Genome}","page":"Reference","title":"BinBencherBackend.recursively_delete_child!","text":"Delete a child from the clade tree.\n\n\n\n\n\n","category":"method"},{"location":"reference/#BinBencherBackend.subset!-Tuple{Reference}","page":"Reference","title":"BinBencherBackend.subset!","text":"subset!(\n        ref::Reference;\n        sequences::Function=Returns(true),\n        genomes::Function=Returns(true)\n)::Reference\n\nMutate ref in place, removing genomes and sequences. Keep only sequences S where sequences(S) returns true and genomes G for which genomes(G) returns true.\n\nSee also: subset, Reference\n\nExamples\n\njulia> ref\nReference\n  Genomes:    3\n  Sequences:  11\n  Ranks:      3\n  Seq length: 10\n  Assembled:  61.9 %\n\njulia> subset(ref; genomes=g -> Flags.organism in flags(g))\nReference\n  Genomes:    2\n  Sequences:  11\n  Ranks:      3\n  Seq length: 10\n  Assembled:  91.3 %\n\njulia> BinBencherBackend.subset(ref; sequences=s -> length(s) ≥ 25)\nReference\n  Genomes:    3\n  Sequences:  9\n  Ranks:      3\n  Seq length: 25\n  Assembled:  56.2 %\n\n\n\n\n\n","category":"method"},{"location":"reference/#BinBencherBackend.subset-Tuple{Reference}","page":"Reference","title":"BinBencherBackend.subset","text":"subset(ref::Reference; kwargs...)\n\nNon-mutating copying version of subset!. This is currently much slower than subset!.\n\nSee also: subset!\n\n\n\n\n\n","category":"method"},{"location":"reference/#BinBencherBackend.update_matrix!-Tuple{Matrix{<:Integer}, Vector{<:AbstractFloat}, Vector{Float64}}","page":"Reference","title":"BinBencherBackend.update_matrix!","text":"For each precision column in the matrix, add one to the correct row given by the recall value at the given precision\n\n\n\n\n\n","category":"method"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"DocTestSetup = quote\n    using BinBencherBackend\n\n    (path_to_ref_file, path_to_bins_file) = let\n        dir = joinpath(Base.pkgdir(BinBencherBackend), \"files\")\n        (joinpath(dir, \"ref.json\"), joinpath(dir, \"clusters.tsv\"))\n    end\nend","category":"page"},{"location":"walkthrough/#Example-walkthrough","page":"Walkthrough","title":"Example walkthrough","text":"","category":"section"},{"location":"walkthrough/#Loading-the-reference","page":"Walkthrough","title":"Loading the reference","text":"","category":"section"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"First, let's load the reference:","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"julia> using BinBencherBackend\n\njulia> ref = Reference(path_to_ref_file)\nReference\n  Genomes:    3\n  Sequences:  11\n  Ranks:      3\n  Seq length: 10\n  Assembled:  61.9 %","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"This gives us a few statistics about the reference:","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"Number of genomes\nNumber of sequences\nNumber of taxonomic ranks (strain, species, genus...)\nLength of shortest sequence\nTotal length of genomes that are assembled","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"The genomes here contain both plasmids and organisms. Let's filter the reference using subset! to only retain organisms, and sequences of length 10 or more:","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"julia> subset!(ref;\n           genomes = is_organism,\n           sequences = s -> length(s) >= 10\n       )\nReference\n  Genomes:    2\n  Sequences:  11\n  Ranks:      3\n  Seq length: 10\n  Assembled:  91.3 %","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"note: Note\nThe function subset! will mutate the reference, whereas the function subset will create a new independent reference. At the moment, the latter is much slower.","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"We removed a single genome, namely one labeled as virus.","category":"page"},{"location":"walkthrough/#Genomes","page":"Walkthrough","title":"Genomes","text":"","category":"section"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"We can explore the genomes contained in the reference with the genomes function, which returns an iterable of Genome (in this case, a Set):","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"julia> genomes(ref)\nSet{Genome} with 2 elements:\n  Genome(gA)\n  Genome(gB)","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"Let's look at a Genome in more detail:","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"julia> genome, genome2 = genomes(ref);\n\njulia> genome\nGenome \"gA\"\n  Parent:        \"D\"\n  Genome size:   100\n  Assembly size: 88 (88.0 %)\n  Sources:       1\n  Flags:         1 (organism)","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"The flags can be extracted with the flags(genome) function - each genome contains zero or more flags:","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"julia> flags(genome)\nFlagSet with 1 element:\n  BinBencherBackend.Flags.organism","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"... in this case, this genome is an organism as opposed to a plasmid or virus. You can see all possible flags with instances(Flags.Flag).","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"See also the helper functions is_organism, is_virus and is_plasmid","category":"page"},{"location":"walkthrough/#Sources","page":"Walkthrough","title":"Sources","text":"","category":"section"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"The genome has one source - let's look at that","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"julia> source = only(genome.sources)\nSource \"subjA1\"\ngenome:          Genome(gA)\n  Length:        100\n  Assembly size: 88\n  Sequences:     6","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"A Source is one of the genomic sequences that genomes are composed of. This is distinct from the assembled sequences that we will be binning - a Source represents the reference sequence, typically the full genome, assembled from a sequencing run on a clonal colony. For this genome, we can see it has a length of 100 bp, and that 5 sequences map to this source, covering 88 bp.","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"We can get the sequences mapping to this source:","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"julia> source.sequences\n6-element Vector{Tuple{Sequence, Tuple{Int64, Int64}}}:\n (Sequence(\"s1\", 25), (5, 29))\n (Sequence(\"s1\", 25), (10, 34))\n (Sequence(\"s2\", 40), (1, 40))\n (Sequence(\"s3\", 50), (51, 98))\n (Sequence(\"s7\", 20), (21, 40))\n (Sequence(\"s8\", 25), (2, 26))","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"Where, e.g. the first entrance tells us that the sequence \"s2\" with a length of 40 maps to positions 1:40 (both inclusive).","category":"page"},{"location":"walkthrough/#Clades","page":"Walkthrough","title":"Clades","text":"","category":"section"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"Genomes are organised into a taxonomic hierarchy. We can find the immediate parent of a genome by accessing the field genome.parent. Let's look at another genome:","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"julia> genome2\nGenome \"gB\"\n  Parent:        \"D\"\n  Genome size:   50\n  Assembly size: 49 (98.0 %)\n  Sources:       2\n  Flags:         1 (organism)\n\njulia> clade = genome2.parent\nSpecies \"D\", 2 genomes\n├─ Genome(gA)\n└─ Genome(gB)","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"The parent is an instance of a Clade. Clades are at a specific rank: Rank 1 for species, 2 for genus, 3 for family, etc. Every clade has one or more children: These are the clades one rank lower. Conceptually, rank zero corresponds to Genomes (OTUs, for this reference dataset)","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"julia> clade.children\n2-element Vector{Genome}:\n Genome(gA)\n Genome(gB)","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"We can find the most recent common ancestor (MRCA) of genome and genome2 like this:","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"julia> mrca(genome, genome2)\nSpecies \"D\", 2 genomes\n├─ Genome(gA)\n└─ Genome(gB)","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"They are very distantly related, so the domain \"Bacteria\", one of the highest ranked Clades, are their most recent common ancestor.","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"The top clade can be found with the top_clade(ref) function, which is the universal ancestor of all clades in the reference.","category":"page"},{"location":"walkthrough/#Binnings","page":"Walkthrough","title":"Binnings","text":"","category":"section"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"A Binning is a set of bins benchmarked against a reference. We can load a set of Vamb bins and turn it into a Binning object like this:","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"julia> binning = Binning(path_to_bins_file, ref)\nBinning\n  Reference\n    Genomes:    2\n    Sequences:  11\n    Ranks:      3\n    Seq length: 10\n    Assembled:  91.3 %\n  Bins:        6\n  NC genomes:  0\n  HQ bins:     0\n  Mean bin genome   R/P/F1: 0.51 / 1.0 / 0.672\n  Mean bin assembly R/P/F1: 0.546 / 1.0 / 0.704\n  Precisions: [0.6, 0.7, 0.8, 0.9, 0.95, 0.99]\n  Recalls:    [0.6, 0.7, 0.8, 0.9, 0.95, 0.99]\n  Reconstruction (genomes):\n    P\\R   0.6  0.7  0.8  0.9 0.95 0.99\n    0.6     1    0    0    0    0    0\n    0.7     1    0    0    0    0    0\n    0.8     1    0    0    0    0    0\n    0.9     1    0    0    0    0    0\n    0.95    1    0    0    0    0    0\n    0.99    1    0    0    0    0    0","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"A wealth of information is readily available:","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"binning.ref gives the underlying Reference\nbinning.recalls and binning.precisions gives the recall/precision thresholds used in benchmarking","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"The function print_matrix will display the number of recovered genomes/assemblies. It takes two optional keyword: level, the taxonomic rank (defaults to 0, meaning strain level), and assembly which defaults to true. If set to false, benchmark number of recovered genomes, not number of recovered assemblies.","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"julia> print_matrix(binning; level=1, assembly=false)\nP\\R   0.6  0.7  0.8  0.9 0.95 0.99\n0.6     1    0    0    0    0    0\n0.7     1    0    0    0    0    0\n0.8     1    0    0    0    0    0\n0.9     1    0    0    0    0    0\n0.95    1    0    0    0    0    0\n0.99    1    0    0    0    0    0","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"You can also get the number of genomes or assemblies reconstructed at a given precision/recall level directly with n_recovered:","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"julia> n_recovered(binning, 0.6, 0.7; assembly=true)\n1\n\njulia> n_recovered(binning, 0.66, 0.91; level=1)\n0","category":"page"},{"location":"walkthrough/#Bins","page":"Walkthrough","title":"Bins","text":"","category":"section"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"The Binning object obviously contains our bins. Let's pick a random bin:","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"julia> bin = binning.bins[4]\nBin \"C4\"\n  Sequences: 3\n  Breadth:   55\n  Intersecting 2 genomes","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"The \"breadth\" here is the sum of the length of its sequences. bin.sequences gets an iterable of its sequences:","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"julia> collect(bin.sequences)\n3-element Vector{Sequence}:\n Sequence(\"s5\", 25)\n Sequence(\"s6\", 10)\n Sequence(\"s7\", 20)","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"The \"Intersecting 2 genomes\" means that the sequences map to 2 different genomes - the only two in the reference. We can get that with the function intersecting, then get the precision/recall with recall_precision:","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"julia> Set(intersecting(bin)) == genomes(ref)\ntrue\n\njulia> recall_precision(genome2, bin)\n(recall = 0.6, precision = 1.0)","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"Or do the same for a higher clade - let's say a genus. In this case, we get the same result.","category":"page"},{"location":"walkthrough/","page":"Walkthrough","title":"Walkthrough","text":"julia> genus = only(Iterators.filter(i -> i.rank == 2, intersecting(Clade, bin)));\n\njulia> recall_precision(genus, bin)\n(recall = 0.6, precision = 1.0)","category":"page"},{"location":"#BinBencherBackend","page":"Home","title":"BinBencherBackend","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"BinBencherBackend.jl is a package for efficient benchmarking and interactive exploration of a set of metagenomic assembled genomes (MAGs) against a reference. This is designed to be used for benchmarking metagenomic binners against a simulated metagenome.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Install Julia - preferably using juliaup: https://github.com/JuliaLang/juliaup\nLaunch Julia: julia\nPress ] to enter package mode. You can exit package mode with backspace.\nIn package mode, type add BinBencherBackend to download and install the benchmarking software","category":"page"},{"location":"#Quickstart","page":"Home","title":"Quickstart","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"using BinBencherBackend\nref =  Reference(\"files/ref.json\")\nbins = Binning(\"files/clusters.tsv\", ref)\nprint_matrix(bins)","category":"page"},{"location":"#Concepts","page":"Home","title":"Concepts","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A Sequence is a sequence (e.g. contig) clustered by the binner\nA Genome is a target genome that should be reconstructed by the binner. It can be a virus, organism, plasmid etc. Every Genome have several Sources, and one parent Clade.\nA Flag marks the certaincy about a boolean attribute of a genome, like \"is this a virus?\".\nSources are the sequences that Genomes are composed of. These are typically the reference genome sequences originally obtained by assembly of a purified genome (e.g. clonal colony). Sequences map to zero or more Sources at particular spans, i.e. locations.\nA Clade contain one or more Genomes or Clades. Clades containing genomes are rank 1, and clades containing rank N clades are rank N+1 clades. All genomes descend from a chain of exactly N ranks of clades, where N > 0.\nA Bin is a set of Sequences created by the binner. Every bin is benchmarked against all genomes and clades in the reference.\nA Reference is composed of:\nThe genomes, a set of Genomes, each with a set of Sources and Flags\nThe taxmaps, a full set of Clades that encompasses every Genome at N ranks (where N > 0)\nThe sequences, a list of Sequences, each with zero or more mappings to Sources.\nA Binning is a set of Bins benchmarked against a Reference","category":"page"},{"location":"","page":"Home","title":"Home","text":"See the Reference in the left sidebar.","category":"page"}]
}
