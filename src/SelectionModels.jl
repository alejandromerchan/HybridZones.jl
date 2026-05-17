module SelectionModels

using ..GeneticArchitectures: OneLocusDiploid, dominance

export SelectionModel, FrequencyDependentSelection, SemiDominantFrequencyDependentSelection
export select!, selection_coefficient

"""
    SelectionModel

Abstract supertype for selection model specifications. Concrete subtypes carry
selection parameters and implement `select!` to compute post-selection genotype
frequencies from pre-selection genotype frequencies.
"""
abstract type SelectionModel end

"""
    FrequencyDependentSelection <: SelectionModel

Frequency-dependent selection model following Mallet & Barton (1989, eq. 1–3).

Fitness of each phenotype is:

    w_i = 1 - s × (1 - p_i)

where `p_i` is the frequency of phenotype `i` in the local deme and `s` is the
selection coefficient. Rare phenotypes have lower fitness because predators do
not recognise them as aposematic warning patterns. Common phenotypes are
favoured, creating positive frequency-dependent selection that stabilises
Müllerian mimicry clines in hybrid zones.

# Fields
- `s::Float64`: selection coefficient in [0, 1]

# Examples
```jldoctest
julia> using HybridZones

julia> model = FrequencyDependentSelection(s=0.3)
FrequencyDependentSelection(s=0.3)

julia> model isa SelectionModel
true

julia> selection_coefficient(model)
0.3
```
"""
struct FrequencyDependentSelection <: SelectionModel
    s::Float64
end

function FrequencyDependentSelection(; s)
    0 ≤ s ≤ 1 ||
        throw(ArgumentError("s must be in [0, 1], got $s"))
    return FrequencyDependentSelection(Float64(s))
end

function Base.show(io::IO, model::FrequencyDependentSelection)
    return print(io, "FrequencyDependentSelection(s=$(model.s))")
end

"""
    selection_coefficient(model::FrequencyDependentSelection) -> Float64

Return the selection coefficient `s` of `model`.
"""
selection_coefficient(model::FrequencyDependentSelection) = model.s

"""
    select!(freq_next, freq_current, model::FrequencyDependentSelection, arch::OneLocusDiploid)

Apply one generation of frequency-dependent selection in-place.

Computes post-selection genotype frequencies from `freq_current` and writes
them into `freq_next`. Both matrices must have the same size.

## Matrix layout convention

State matrices use the convention `rows = genotypes, columns = demes`.
For a `OneLocusDiploid` architecture, this means a 3×n_demes matrix:

| Row | Genotype |
|-----|----------|
| 1   | A1A1     |
| 2   | A1A2     |
| 3   | A2A2     |


This layout is chosen for cache efficiency: Julia's column-major storage
means each deme's three genotype frequencies are contiguous in memory.
Users coming from R or Python may need to transpose their input data,
since those languages typically use rows-as-observations conventions.

## Genotype-to-phenotype mapping

The mapping from genotypes to phenotypes is determined by `arch.dominance`:

- `:codominant`: each genotype is its own phenotype (three phenotypes total)
- `:dominant_first`: A1A1 and A1A2 share the A1 phenotype; A2A2 is distinct
- `:recessive_first`: A1A1 is distinct; A1A2 and A2A2 share the A2 phenotype

## Contract

- `freq_next` is overwritten; `freq_current` is read-only.
- Each column of `freq_current` must sum to 1.0 (valid genotype frequencies).
- After return, each column of `freq_next` sums to 1.0 within floating-point
  precision.

# Examples
```jldoctest
julia> using HybridZones

julia> model = FrequencyDependentSelection(s=0.0);

julia> arch = OneLocusDiploid();

julia> freq = reshape([0.25, 0.50, 0.25], 3, 1);

julia> freq_next = similar(freq);

julia> select!(freq_next, freq, model, arch);

julia> freq_next ≈ freq  # s=0 → selection has no effect
true
```
"""
function select!(
        freq_next::AbstractMatrix,
        freq_current::AbstractMatrix,
        model::FrequencyDependentSelection,
        arch::OneLocusDiploid
)
    size(freq_next) == size(freq_current) ||
        throw(DimensionMismatch("freq_next and freq_current must have the same size"))
    size(freq_current, 1) == 3 ||
        throw(ArgumentError(
            "expected 3 genotype rows (A1A1, A1A2, A2A2), got $(size(freq_current, 1))"
        ))
    n_demes = size(freq_current, 2)
    s = model.s
    dom = dominance(arch)
    @inbounds for j in 1:n_demes
        p11 = freq_current[1, j]
        p12 = freq_current[2, j]
        p22 = freq_current[3, j]
        if dom === :codominant
            w11 = 1 - s * (1 - p11)
            w12 = 1 - s * (1 - p12)
            w22 = 1 - s * (1 - p22)
            w_mean = w11 * p11 + w12 * p12 + w22 * p22
            freq_next[1, j] = w11 * p11 / w_mean
            freq_next[2, j] = w12 * p12 / w_mean
            freq_next[3, j] = w22 * p22 / w_mean
        elseif dom === :dominant_first
            pA = p11 + p12
            pa = p22
            wA = 1 - s * (1 - pA)
            wa = 1 - s * (1 - pa)
            w_mean = wA * pA + wa * pa
            freq_next[1, j] = wA * p11 / w_mean
            freq_next[2, j] = wA * p12 / w_mean
            freq_next[3, j] = wa * p22 / w_mean
        else  # :recessive_first
            p1 = p11
            p2 = p12 + p22
            w1 = 1 - s * (1 - p1)
            w2 = 1 - s * (1 - p2)
            w_mean = w1 * p1 + w2 * p2
            freq_next[1, j] = w1 * p11 / w_mean
            freq_next[2, j] = w2 * p12 / w_mean
            freq_next[3, j] = w2 * p22 / w_mean
        end
    end
    return freq_next
end

"""
    SemiDominantFrequencyDependentSelection(s)

Frequency-dependent selection with semi-dominant phenotype mixing.

This model implements the selection regime used in Mallet's WC1SEDO
Pascal program. Unlike `FrequencyDependentSelection`, which treats
each genotype as its own distinct phenotype under `:codominant`
dominance, this model treats heterozygotes as having partial
similarity to both homozygous parental phenotypes.

For each genotype, a "perceived phenotype frequency" is computed
that includes contributions from related phenotypes via semi-dominance
mixing:

- Each homozygote receives a contribution of half the heterozygote
  frequency (the heterozygote looks half-like each parent)
- The heterozygote receives a contribution of half of each homozygote
  frequency (symmetric, since heterozygote looks half-like both
  parents)

The fitness of each genotype is then computed from its perceived
phenotype frequency:

    w_G = 1 - s * (1 - sim_G)

Selection acts proportionally on each genotype, with normalization
by mean fitness.

# Arguments
- `s`: selection coefficient, in [0, 1]

# Biological interpretation
Semi-dominant frequency-dependent selection is biologically appropriate
for warning-color systems where heterozygotes have visually intermediate
phenotypes that share partial similarity with both parental forms,
and predators have partial recognition across similar patterns. This
is the model used by Mallet & Barton (1989) and the broader Pascal
simulator tradition for *Heliconius* hybrid zones where the
heterozygote wing pattern is intermediate between the parental races.

For systems where heterozygotes have a distinctly different phenotype
recognized as its own pattern class, use `FrequencyDependentSelection`
with `:codominant` dominance instead.

# Dispatch
This model dispatches only on `OneLocusDiploid` with `:codominant`
dominance. Combining it with `:dominant_first` or `:recessive_first`
will raise a `MethodError`. The semi-dominance is intrinsic to the
selection model itself.

# Examples
```jldoctest
julia> using HybridZones

julia> arch = OneLocusDiploid(dominance = :codominant);

julia> sel = SemiDominantFrequencyDependentSelection(0.3);
```
"""
struct SemiDominantFrequencyDependentSelection <: SelectionModel
    s::Float64

    function SemiDominantFrequencyDependentSelection(s::Real)
        0 ≤ s ≤ 1 ||
            throw(ArgumentError("s must be in [0, 1], got $s"))
        return new(Float64(s))
    end
end

function Base.show(io::IO, model::SemiDominantFrequencyDependentSelection)
    return print(io, "SemiDominantFrequencyDependentSelection(s=$(model.s))")
end

"""
    selection_coefficient(model::SemiDominantFrequencyDependentSelection) -> Float64

Return the selection coefficient `s` of `model`.
"""
selection_coefficient(model::SemiDominantFrequencyDependentSelection) = model.s

"""
    select!(g_next, g_current, model::SemiDominantFrequencyDependentSelection, arch::OneLocusDiploid)

Apply one generation of semi-dominant frequency-dependent selection in-place.

Computes post-selection genotype frequencies from `g_current` using Pascal's
semi-dominant phenotype mixing (WC1SEDO.PAS lines 720–735) and writes them
into `g_next`. Both matrices must have the same size. Only `:codominant`
dominance is supported; other modes raise an `ArgumentError`.

## Semi-dominant mixing

The perceived phenotype frequency for each genotype includes contributions
from related phenotypes:

    sim_A1A1 = p_A1A1 + 0.5 * p_A1A2
    sim_A1A2 = p_A1A2 + 0.5 * (p_A1A1 + p_A2A2)
    sim_A2A2 = p_A2A2 + 0.5 * p_A1A2

Fitness for each genotype is `w_G = 1 - s * (1 - sim_G)`.

See `FrequencyDependentSelection` for the strict codominant alternative
where each genotype is its own distinct phenotype.

# Examples
```jldoctest
julia> using HybridZones

julia> model = SemiDominantFrequencyDependentSelection(0.0);

julia> arch = OneLocusDiploid(dominance = :codominant);

julia> freq = reshape([0.25, 0.50, 0.25], 3, 1);

julia> freq_next = similar(freq);

julia> select!(freq_next, freq, model, arch);

julia> freq_next ≈ freq  # s=0 → selection has no effect
true
```
"""
function select!(
        g_next::AbstractMatrix,
        g_current::AbstractMatrix,
        model::SemiDominantFrequencyDependentSelection,
        arch::OneLocusDiploid
)
    dominance(arch) === :codominant ||
        throw(ArgumentError(
            "SemiDominantFrequencyDependentSelection requires :codominant " *
            "dominance, got :$(dominance(arch)). Semi-dominance is " *
            "intrinsic to this model; other dominance modes are not meaningful."
        ))
    size(g_next) == size(g_current) ||
        throw(DimensionMismatch("g_next and g_current must have the same size"))
    size(g_current, 1) == 3 ||
        throw(ArgumentError(
            "expected 3 genotype rows (A1A1, A1A2, A2A2), got $(size(g_current, 1))"
        ))

    n_demes = size(g_current, 2)
    s = model.s

    @inbounds for j in 1:n_demes
        p11 = g_current[1, j]
        p12 = g_current[2, j]
        p22 = g_current[3, j]

        # Semi-dominant phenotype mixing (Pascal SimilarPhenotypes)
        sim_11 = p11 + 0.5 * p12
        sim_12 = p12 + 0.5 * (p11 + p22)
        sim_22 = p22 + 0.5 * p12

        # Fitnesses based on perceived phenotype frequencies
        w11 = 1 - s * (1 - sim_11)
        w12 = 1 - s * (1 - sim_12)
        w22 = 1 - s * (1 - sim_22)

        # Mean fitness
        w_mean = w11 * p11 + w12 * p12 + w22 * p22

        # Post-selection frequencies (proportional, normalized)
        g_next[1, j] = w11 * p11 / w_mean
        g_next[2, j] = w12 * p12 / w_mean
        g_next[3, j] = w22 * p22 / w_mean
    end

    return g_next
end

end # module SelectionModels
