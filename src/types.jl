# Genome, Clade and Source mutually refer to each other, so they are defined
# together in a typegroup. Their docstrings are in their respective files.

# Values derived from the reference, computed when the reference is finished.
# Each object stores these in a single `stats` field which is `nothing` until
# the object is finished, so the derived values are always either all present or all absent.
struct GenomeStats
    genome_size::Int
    assembly_size::Int
end

struct SourceStats
    assembly_size::Int
    total_bp::Int
end

typegroup
    mutable struct Genome
        const name::String
        const sources::Set{Source}
        const flags::FlagSet
        # Nothing if not yet added to a clade. A finished genome always has a parent.
        parent::Union{Clade, Nothing}
        stats::Union{GenomeStats, Nothing}

        function Genome(name::AbstractString, flags::FlagSet)
            return new(
                check_valid_identifier(String(name)),
                Set{Source}(),
                flags,
                nothing,
                nothing,
            )
        end
    end

    mutable struct Clade
        const name::String
        const rank::Int
        const children::Union{Vector{Clade}, Vector{Genome}}
        ngenomes::Int
        parent::Union{Clade, Nothing}

        function Clade(name::String, child::Union{Clade, Genome})
            check_valid_identifier(name)
            existing_parent = child.parent
            if existing_parent !== nothing
                existing_parent_error(child.name, name, existing_parent.name)
            end
            (rank, ngenomes) = child isa Genome ? (1, 1) : (child.rank + 1, child.ngenomes)
            instance = new(name, rank, [child], ngenomes, nothing)
            child.parent = instance
            return instance
        end
    end

    mutable struct Source
        const name::String
        const genome::Genome
        const length::Int
        const sequences::Vector{Tuple{Sequence, Tuple{Int, Int}}}
        stats::Union{SourceStats, Nothing}

        function Source(genome::Genome, name::AbstractString, length::Integer)
            length ≤ 0 && throw(ArgumentError("Source length must be at least 1"))
            return new(
                check_valid_identifier(String(name)),
                genome,
                Int(length),
                Vector{Tuple{Sequence, Tuple{Int, Int}}}(),
                nothing,
            )
        end
    end
end

@noinline function unfinished_error(x)
    throw(
        ArgumentError(
            "$(x) is not finished. Its derived values are computed when its `Reference` is loaded.",
        ),
    )
end
