# See the comment on GenomeStats
struct ReferenceStats
    shortest_seq_len::Int
    fraction_assembled::Float64
end

mutable struct Reference
    const genomes::Set{Genome}
    # Sequence name => index into targets
    const target_index_by_name::Dict{String, UInt32}
    const targets::Vector{Tuple{Sequence, Vector{Target}}}
    const clades::Vector{Vector{Clade}}
    stats::Union{ReferenceStats, Nothing}
end

"""
    Reference(::Union{IO, AbstractString})

A `Reference` contains the ground truth to benchmark against.
Conceptually, it consists of the following parts:
* A list of genomes, each with sources
* The full taxonomic tree, as lists of clades
* A list of sequences, each with a list of (source, span) to where it maps.

Normally, the types `FlagSet` `Genome`, `Source`, `Clade` and `Sequence` do not
need to be constructed manually, but are constructed when the `Reference` is loaded
from a JSON file.

A `Reference` is loaded from a JSON file, which is specified either as an `IO`,
or its path as an `AbstractString`. If the path ends with `.gz`, automatically
gzip decompress when reading the file.

# Examples
```jldoctest
julia> ref = Reference(path_to_ref_file);


julia> ref isa Reference
true

julia> length(genomes(ref))
3

julia> n_seqs(ref)
11

julia> first(ref.genomes) isa Genome
true
```

See also: [`subset`](@ref), [`Genome`](@ref), [`Clade`](@ref)
"""
Reference

function new_reference()
    return Reference(
        Set{Genome}(),
        Dict{String, Int}(),
        Vector{Tuple{Sequence, Vector{Target}}}(),
        Vector{Clade}[],
        nothing,
    )
end

Base.show(io::IO, ::Reference) = print(io, "Reference()")
function Base.show(io::IO, ::MIME"text/plain", x::Reference)
    return if get(io, :compact, false)
        show(io, x)
    else
        print(
            io,
            "Reference",
            "\n  Genomes:    ",
            length(genomes(x)),
            "\n  Sequences:  ",
            n_seqs(x),
            "\n  Ranks:      ",
            nranks(x),
            "\n  Seq length: ",
            shortest_seq_len(x),
            "\n  Assembled:  ",
            round(fraction_assembled(x) * 100; digits = 1),
            " %",
        )
    end
end

top_clade(x::Reference) = only(last(x.clades))
genomes(x::Reference) = x.genomes
n_seqs(x::Reference) = length(x.targets)
function nranks(x::Reference)
    isempty(x.genomes) && return 0
    return length(x.clades) + 1
end

function finished_stats(x::Reference)::ReferenceStats
    stats = x.stats
    stats === nothing && unfinished_error(x)
    return stats
end

"""
    shortest_seq_len(ref::Reference)::Int

Get the length of the shortest sequence in the `Reference`.

# Examples
```jldoctest
julia> shortest_seq_len(ref)
10
```
"""
shortest_seq_len(x::Reference) = finished_stats(x).shortest_seq_len

"""
    fraction_assembled(ref::Reference)::Float64

Get the fraction of the total genome size of all genomes in the `Reference`
which is covered by at least one sequence.

See also: [`assembly_size`](@ref), [`genome_size`](@ref)

# Examples
```jldoctest
julia> round(fraction_assembled(ref); digits = 3)
0.619
```
"""
fraction_assembled(x::Reference) = finished_stats(x).fraction_assembled

function finish!(ref::Reference)
    ref.stats === nothing || return ref
    # These are checked by the reference parser and by `subset!`, so if this fails, it's a bug
    isempty(ref.genomes) && error("Reference cannot have zero genomes")
    scratch = Tuple{Int, Int}[]
    foreach(g -> finish!(g, scratch), ref.genomes)
    total_assembly_size = total_genome_size = 0
    for genome in ref.genomes
        total_assembly_size += assembly_size(genome)
        total_genome_size += genome_size(genome)
    end
    shortest = minimum(i -> length(first(i)), ref.targets; init = typemax(Int))
    shortest == typemax(Int) && error("Cannot initialize a Reference with no sequences")
    ref.stats = ReferenceStats(shortest, total_assembly_size / total_genome_size)
    return ref
end

"""
    subset!(
            ref::Reference;
            sequences::Function=Returns(true),
            genomes::Function=Returns(true)
    )::Reference

Mutate `ref` in place, removing genomes and sequences.
Keep only sequences S where `sequences(S)` returns `true` and genomes G for which
`genomes(G)` returns `true`.

See also: [`subset`](@ref), [`Reference`](@ref)

# Examples
```jldoctest
julia> ref
Reference
  Genomes:    3
  Sequences:  11
  Ranks:      3
  Seq length: 10
  Assembled:  61.9 %

julia> subset(ref; genomes=g -> Flags.organism in flags(g))
Reference
  Genomes:    2
  Sequences:  11
  Ranks:      3
  Seq length: 10
  Assembled:  91.3 %

julia> BinBencherBackend.subset(ref; sequences=s -> length(s) ≥ 25)
Reference
  Genomes:    3
  Sequences:  9
  Ranks:      3
  Seq length: 25
  Assembled:  56.2 %
```
"""
function subset!(
        ref::Reference;
        sequences::Function = Returns(true),
        genomes::Function = Returns(true),
    )::Reference
    # Evaluate the predicates once, before mutating anything, such that we can
    # error without leaving `ref` in a partially subsetted state.
    mask = BitVector(sequences(first(i)) for i in ref.targets)
    genomes_to_remove = Genome[g for g in ref.genomes if !(genomes(g)::Bool)]
    if length(genomes_to_remove) == length(ref.genomes)
        throw(ArgumentError("Subsetting would leave the reference with zero genomes"))
    end
    if !any(mask)
        throw(ArgumentError("Subsetting would leave the reference with no sequences"))
    end

    ref = unfinish!(ref)
    sources_to_remove = Set{Source}()

    # Update both ref.targets and ref.target_index_by_name
    new_idx = cumsum(mask)
    keepat!(ref.targets, mask)
    filter!(kv -> mask[last(kv)], ref.target_index_by_name)
    map!(i -> new_idx[i], values(ref.target_index_by_name))

    # Remove the genomes from ref.genomes
    setdiff!(ref.genomes, genomes_to_remove)

    # Compute sources to remove, and remove the sources from ref.targets
    for genome in genomes_to_remove
        union!(sources_to_remove, genome.sources)
    end
    for (_, targets) in ref.targets
        filter!(targets) do (source, _)
            !in(source, sources_to_remove)
        end
    end

    # Remove references to deleted genomes from parents
    # (and recursively remove now-empty Clades)
    for genome in genomes_to_remove
        recursively_delete_child!(genome)
    end

    # Remove references to Clades we deleted from the Clade tree,
    # but which may still be present in ref.clades
    for i in (length(ref.clades) - 1):-1:1
        empty!(ref.clades[i])
        for parent in ref.clades[i + 1]
            union!(ref.clades[i], parent.children::Vector{Clade})
        end
    end

    # Remove sequences from sources. We look up the sequences kept above instead
    # of re-evaluating the predicate, which may not be deterministic.
    for genome in ref.genomes, source in genome.sources
        filter!(i -> haskey(ref.target_index_by_name, first(i).name), source.sequences)
    end

    # Re-initialize the now completely filtered Reference
    return finish!(ref)
end

# This deepcopy is quite slow, and it would be nice to optimise it.
# However, manual deepcopying of references is quite error prone, and having
# an incomplete deepcopy could lead to nasty bugs, so I'll just eat it.
"""
    subset(ref::Reference; kwargs...)

Non-mutating copying version of `subset!`.
This is currently much slower than `subset!`.

See also: [`subset!`](@ref)
"""
subset(ref::Reference; kwargs...) = subset!(deepcopy(ref); kwargs...)

function unfinish!(ref::Reference)
    ref.stats = nothing
    for genome in ref.genomes
        genome.stats = nothing
        for source in genome.sources
            source.stats = nothing
        end
    end
    return ref
end

# Only called when parsing a reference file
function add_genome!(ref::Reference, genome::Genome)
    if in(genome, ref.genomes)
        throw(ReferenceError("Genome \"$(genome.name)\" is present more than once"))
    end
    push!(ref.genomes, genome)
    return ref
end

# Only called when parsing a reference file
function add_sequence!(ref::Reference, seq::Sequence, targets::Vector{Target})
    for (source, span) in targets
        add_sequence!(source, seq, span)
    end
    L = length(ref.target_index_by_name)
    # We disallow a sequence with index typemax(UInt32) so we can use that number
    # to encode a missing sequence
    if L == typemax(UInt32) - 1
        throw(ReferenceError("References can only hold 4294967294 sequences"))
    end
    i = get!(ref.target_index_by_name, seq.name, L + 1)
    if i != L + 1
        throw(ReferenceError("Sequence \"$(seq.name)\" is present more than once"))
    end
    push!(ref.targets, (seq, targets))
    return ref
end

const CLUSTER_FILE_HEADER = "clustername\tcontigname"

# Parse a cluster file, returning the indices (into `ref.targets`) of the sequences
# of each bin, and, if the binning is not disjoint, the line and name of the
# first sequence seen in more than one bin.
function parse_bins(
        io::IO,
        ::Type{Dict},
        ref::Reference,
        binsplit_sep::Union{Nothing, AbstractString, Char} = nothing,
    )::Tuple{Dict{SubString{String}, Vector{UInt32}}, Union{Nothing, Tuple{Int, String}}}
    seen_indices = falses(length(ref.targets))
    idxs_by_binname = Dict{SubString{String}, Vector{UInt32}}()
    first_non_disjoint = nothing
    seen_header = false
    for (line_number, full_line) in enumerate(eachline(io))
        line = rstrip(full_line)
        if !seen_header
            if line != CLUSTER_FILE_HEADER
                throw(
                    BinningError(
                        "Expected header $(repr(CLUSTER_FILE_HEADER)), got $(repr(line))", line_number
                    ),
                )
            end
            seen_header = true
            continue
        end
        isempty(line) && continue
        (binname, seqname) = parse_cluster_line(line, line_number)
        if binsplit_sep !== nothing
            (before, _) = @something split_once(seqname, binsplit_sep) throw(
                BinningError(
                    "Binsplit separator $(repr(binsplit_sep)) not found in sequence name $(repr(seqname))",
                    line_number,
                ),
            )
            joined = string(before, binname)
            binname = SubString(joined, 1:lastindex(joined))
        end
        i = get(ref.target_index_by_name, seqname, nothing)
        if i === nothing
            throw(BinningError("Sequence $(repr(seqname)) is not in the reference", line_number))
        end
        if seen_indices[i] && first_non_disjoint === nothing
            first_non_disjoint = (line_number, String(seqname))
        end
        seen_indices[i] = true
        push!(get!(valtype(idxs_by_binname), idxs_by_binname, binname), i)
    end
    seen_header || throw(BinningError("Cluster file is empty, expected a header", 1))
    # A sequence may be in multiple bins (if not disjoint), but never twice in the same bin
    for (binname, idxs) in idxs_by_binname
        if !allunique(idxs)
            seen = Set{UInt32}()
            dup = first(i for i in idxs if in!(seen, i))
            seqname = first(ref.targets[dup]).name
            throw(
                BinningError(
                    "Sequence \"$(seqname)\" is present multiple times in bin \"$(binname)\"", nothing
                ),
            )
        end
    end
    return (idxs_by_binname, first_non_disjoint)
end

# Split a cluster file line into (bin name, sequence name)
function parse_cluster_line(line::SubString{String}, line_number::Int)
    cu = codeunits(line)
    t1 = findfirst(isequal(UInt8('\t')), cu)
    t1 === nothing && throw(BinningError("Expected two tab-separated fields, found one", line_number))
    t2 = findnext(isequal(UInt8('\t')), cu, t1 + 1)
    if t2 !== nothing
        throw(BinningError("Expected two tab-separated fields, found more", line_number))
    end
    binname = SubString(line, 1:prevind(line, t1))
    seqname = SubString(line, (t1 + 1):lastindex(line))
    # Tabs and newlines can't be present, but a \r could be
    for name in (binname, seqname)
        if !is_valid_bb_identifier(name)
            throw(BinningError(invalid_identifier_message(name), line_number))
        end
    end
    return (binname, seqname)
end

const JSON_VERSION = 2

# [(genome_name, flags, [(sourcename, length)...])...]
# Where flags is the integer representation of Flag,
# i.e. each bit from lowest to highest in the integer
# represents one flag, in the order of definition
# Example: A genome called "MyGenome" has one subject, "C1" of length 1000,
# and is both a plasmid and a virus.
# So: Flags are defined as: Organism = 1, virus = 2, plasmid = 4
# virus + plasmid = 2 | 4 = 6
# So: [("MyGenome", 6, [("C1", 1000)])]
const GENOMES_JSON_T = Vector{Tuple{String, Int, Vector{Tuple{String, Int}}}}

# [(sequence_name,  sequence_length, [(subject_name, from_position, to_position)...])...]
# The positions are the 1-indexed starting and ending positions of where the
# sequence maps to in the subject.
# If the sequence 'wraps' around the subject (i.e. subject is circular), then
# you have from_position > to_position
# Example: Seq S1C1 with length 25 maps to C1 at pos 995 to 19 (wrapping around)
# [("S1C1", 25, [("C1", 995, 19)])]
const SEQUENCES_JSON_T = Vector{Tuple{String, Int, Vector{Tuple{String, Int, Int}}}}

# [[(child_name, parent_name)...], [(parent_name, grandparent_name)...] ...]
# One inner vector per taxonomic rank.
# Lowest level has genome name as the left element and the genome's parent as
# the right element.
# Example: Genomes G1 and G2 are both E. colis and are the only genomes:
# [[("G1", "E. coli"), ("G2", "E. coli")], [("E. coli", "Escherichia")]]
const TAXMAPS_JSON_T = Vector{Vector{Tuple{String, String}}}

struct ReferenceJSON
    version::Int
    genomes::GENOMES_JSON_T
    sequences::SEQUENCES_JSON_T
    taxmaps::TAXMAPS_JSON_T
end

Reference(path::AbstractString) = open_perhaps_gzipped(i -> Reference(i), String(path))

function Reference(io::IO)
    str = read(io, String)
    # The JSON library throws its own exceptions on invalid JSON, so we catch
    # everything thrown by these library calls (and only these calls).
    # Their only input is the file content, so any exception they throw means that
    # the file is malformed.
    # First, determine the version of the JSON file to see if this version of BBB
    # can even deserialize it.
    (lazy_json, version) = try
        lazy_json = JSON.lazy(str)
        v = get(lazy_json, :version, nothing)
        (lazy_json, v === nothing ? nothing : v[])
    catch e
        e isa InterruptException && rethrow()
        throw(ReferenceError("Malformed JSON: " * sprint(showerror, e)))
    end
    # Note: If the top-level value is an array, `version` is not `nothing`, but
    # neither is it an integer
    if !(version isa Int)
        throw(
            ReferenceError(
                "Malformed JSON: Expected a top-level object with an integer field \"version\""
            ),
        )
    end
    if version != JSON_VERSION
        throw(
            ReferenceError(
                "Found reference version $(version), but the supported version of the " *
                    "currently loaded version of BinBencherBackend is $(JSON_VERSION).",
            ),
        )
    end
    json_struct = try
        JSON.parse(lazy_json, ReferenceJSON)
    catch e
        e isa InterruptException && rethrow()
        throw(ReferenceError("Malformed JSON: " * sprint(showerror, e)))
    end
    return Reference(json_struct)
end

function Reference(json_struct::ReferenceJSON)
    if json_struct.version != JSON_VERSION
        @warn (
            "Deserializing reference JSON of version $(json_struct.version), " *
                "but the supported version of the currently loaded version of BinBencherBackend is $(JSON_VERSION)."
        )
    end
    ref = new_reference()
    isempty(json_struct.genomes) && throw(ReferenceError("Reference cannot have zero genomes"))
    if isempty(json_struct.sequences)
        throw(ReferenceError("Reference cannot have zero sequences"))
    end

    # Parse genomes
    for (genomename, flags, sourcesdict) in json_struct.genomes
        if !is_valid_bb_identifier(genomename)
            throw(ReferenceError(invalid_identifier_message(genomename)))
        end
        flagset = @something try_from_bits(FlagSet, flags) throw(
            ReferenceError(
                "Genome \"$(genomename)\" has invalid flags $(flags), must be in 0:$(MAX_FLAG_BITS)"
            ),
        )
        genome = Genome(genomename, flagset)
        add_genome!(ref, genome)
        for (source_name, source_length) in sourcesdict
            add_source!(genome, source_name, source_length)
        end
    end

    # Parse taxonomy in a background task. Exceptions thrown in a task are wrapped
    # in a TaskFailedException when fetched, so we return a ReferenceError from the
    # task instead, and throw it from this task.
    tax_task = Threads.@spawn try
        parse_taxonomy(ref.genomes, json_struct.taxmaps)
    catch e
        e isa ReferenceError ? e : rethrow()
    end

    # Check for unique sources
    source_by_name = Dict{String, Source}()
    for genome in ref.genomes, source in genome.sources
        if haskey(source_by_name, source.name)
            throw(
                ReferenceError(
                    "Duplicate source: \"$(source.name)\" belongs to both genome " *
                        "\"$(source_by_name[source.name].genome.name)\" and \"$(genome.name)\".",
                ),
            )
        end
        source_by_name[source.name] = source
    end

    # Parse sequences
    for (seq_name, seq_length, targs) in json_struct.sequences
        if !is_valid_bb_identifier(seq_name)
            throw(ReferenceError(invalid_identifier_message(seq_name)))
        end
        if seq_length < 1
            throw(
                ReferenceError("Sequence \"$(seq_name)\" has length $(seq_length), but must be at least 1")
            )
        end
        targets = Target[]
        for (source_name, from, to) in targs
            source = get(source_by_name, source_name, nothing)
            if source === nothing
                throw(
                    ReferenceError(
                        "Sequence \"$(seq_name)\" maps to source \"$(source_name)\", but no such source in reference",
                    ),
                )
            end
            push!(targets, (source, (Int(from), Int(to))))
        end
        seq = Sequence(seq_name, seq_length)
        add_sequence!(ref, seq, targets)
    end

    # Add taxonomy from the previously spawned task
    taxonomy = fetch(tax_task)::Union{ReferenceError, Vector{Vector{Clade}}}
    taxonomy isa ReferenceError && throw(taxonomy)
    copy!(ref.clades, taxonomy)

    # Finalize the reference
    return finish!(ref)
end

function save(io::IO, ref::Reference)
    json_dict = Dict{String, Any}()
    json_dict["version"] = JSON_VERSION
    # Genomes
    genomes::GENOMES_JSON_T = [
        (genome.name, Int(genome.flags.x), [(s.name, s.length) for s in genome.sources]) for genome in ref.genomes
    ]
    json_dict["genomes"] = genomes

    # Sequences
    sequences::SEQUENCES_JSON_T = [
        (
            seq.name,
            seq.length,
            [(source.name, first(span), last(span)) for (source, span) in targets],
        ) for (seq, targets) in ref.targets
    ]

    json_dict["sequences"] = sequences

    # Taxmaps
    taxmaps::TAXMAPS_JSON_T =
        [[(genome.name, parent_clade(genome).name) for genome in ref.genomes]]
    json_dict["taxmaps"] = taxmaps

    for clades in ref.clades
        length(clades) == 1 && isnothing(only(clades).parent) && break
        v = eltype(taxmaps)()
        push!(taxmaps, v)
        for clade in clades
            parent = clade.parent::Clade
            push!(v, (clade.name, parent.name))
        end
    end

    return JSON.json(io, json_dict)
end

# Invariants for the input list:
# The input is a vector of ranks, consisting of (child, parent) names
# On every rank, the children is a unique set (each child present once) of the parent of the previous rank:
#   - On the first rank, the children is a unique set of genome names
#   - On the last rank, the parent names can be arbitrary.
#     If there is more than a single unique name, then a top node will be added
function parse_taxonomy(
        genomes::Set{Genome},
        taxmaps::Vector{Vector{Tuple{String, String}}},
    )::Vector{Vector{Clade}}
    result = Vector{Vector{Clade}}()
    isempty(genomes) && isempty(taxmaps) && return result
    child_by_name = Dict{String, Node}(g.name => g for g in genomes)
    assigned_children_names = Set{String}()
    parent_by_name = empty(child_by_name)
    for (rank, taxmap) in enumerate(taxmaps)
        cladeset = Clade[]
        empty!(assigned_children_names)
        for (child_name, parent_name) in taxmap
            if child_name ∈ assigned_children_names
                throw(ReferenceError("At rank $(rank), child \"$(child_name)\" is present more than once."))
            end
            # Get child
            child = get(child_by_name, child_name, nothing)
            if child === nothing
                throw(
                    ReferenceError(
                        "At rank $(rank), found child name \"$(child_name)\", but this does not exist " *
                            "on the previous rank.",
                    ),
                )
            end
            # Create parent if it does not already exist
            parent =
                get(parent_by_name, parent_name, nothing)::Union{Clade, Nothing}
            if parent === nothing
                if !is_valid_bb_identifier(parent_name)
                    throw(ReferenceError(invalid_identifier_message(parent_name)))
                end
                parent = Clade(parent_name, child)
                parent_by_name[parent_name] = parent
                push!(cladeset, parent)
            else
                add_child!(parent, child)
            end
            push!(assigned_children_names, child_name)
        end
        if length(assigned_children_names) != length(child_by_name)
            missing_child = first(setdiff(keys(child_by_name), assigned_children_names))
            throw(ReferenceError("At rank $(rank), child \"$(missing_child)\" has no parent"))
        end
        parent_by_name, child_by_name = child_by_name, parent_by_name
        empty!(parent_by_name)
        push!(result, cladeset)
    end
    if isempty(taxmaps)
        gns = collect(genomes)
        top = Clade("top", gns[1])
        for genome in @view(gns[2:end])
            add_child!(top, genome)
        end
        return [[top]]
    elseif length(last(result)) > 1
        top = Clade("top", first(last(result)))
        for child in @view(result[end][2:end])
            add_child!(top, child)
        end
        push!(result, [top])
    end
    foreach(i -> sort!(i; by = j -> j.name), result)
    return result
end
