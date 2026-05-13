module MatingModels

using ..GeneticArchitectures: OneLocusDiploid

export MatingModel, RandomMating
export mate!

"""
    MatingModel

Abstract supertype for mating model specifications. Concrete subtypes implement
`mate!` to compute post-mating genotype frequencies from pre-mating genotype
frequencies, typically restoring Hardy-Weinberg proportions within each deme.
"""
abstract type MatingModel end

"""
    RandomMating <: MatingModel

Random mating within each deme, restoring Hardy-Weinberg (HW) proportions
each generation.

This is a zero-field struct: random mating is parameter-free. The type exists
to dispatch `mate!`, not to carry data.

Under random mating, the three diploid genotype frequencies at any deme are
fully determined by the allele frequency `p` (frequency of A1):

    A1A1 = p²
    A1A2 = 2p(1-p)
    A2A2 = (1-p)²

Migration mixing across demes introduces Wahlund-effect heterozygote
deficits relative to HW expectations. Applying `mate!` each generation
restores HW within demes, matching the biological convention used in
Mallet's Pascal simulators and the broader population genetics literature.

# Examples
```jldoctest
julia> using HybridZones

julia> mat = RandomMating()
RandomMating()

julia> mat isa MatingModel
true
```
"""
struct RandomMating <: MatingModel end

"""
    mate!(g_next, g_current, model::RandomMating, arch::OneLocusDiploid)

Apply one generation of random mating in-place, restoring Hardy-Weinberg
proportions within each deme.

For each deme `j`, computes the A1 allele frequency:

    p = g_current[1, j] + 0.5 × g_current[2, j]

then writes HW genotype frequencies into `g_next`:

    g_next[1, j] = p²
    g_next[2, j] = 2p(1-p)
    g_next[3, j] = (1-p)²

Random mating is dominance-independent: the HW restoration formula is the
same regardless of whether the locus is `:codominant`, `:dominant_first`,
or `:recessive_first`. Dispatching on `OneLocusDiploid` rather than the
abstract `GeneticArchitecture` mirrors the `select!` pattern and leaves
room for architecture-specific mating logic in future subtypes.

## Contract

- `g_next` is overwritten; `g_current` is read-only.
- Both matrices must have the same size.
- The first dimension of `g_current` must be 3 (one row per genotype).
- After return, each column of `g_next` sums to 1.0 within floating-point
  precision.

# Examples
```jldoctest
julia> using HybridZones

julia> mat = RandomMating();

julia> arch = OneLocusDiploid();

julia> g = reshape([0.0, 1.0, 0.0], 3, 1);  # all heterozygotes, p=0.5

julia> g_next = similar(g);

julia> mate!(g_next, g, mat, arch);

julia> g_next
3×1 Matrix{Float64}:
 0.25
 0.5
 0.25
```
"""
function mate!(
        g_next::AbstractMatrix,
        g_current::AbstractMatrix,
        ::RandomMating,
        ::OneLocusDiploid
)
    size(g_next) == size(g_current) ||
        throw(DimensionMismatch("g_next and g_current must have the same size"))
    size(g_current, 1) == 3 ||
        throw(ArgumentError(
            "expected 3 genotype rows (A1A1, A1A2, A2A2), got $(size(g_current, 1))"
        ))
    n_demes = size(g_current, 2)
    @inbounds for j in 1:n_demes
        p = g_current[1, j] + 0.5 * g_current[2, j]
        q = 1.0 - p
        g_next[1, j] = p * p
        g_next[2, j] = 2.0 * p * q
        g_next[3, j] = q * q
    end
    return g_next
end

end # module MatingModels
