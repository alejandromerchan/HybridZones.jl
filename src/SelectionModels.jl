module SelectionModels

using ..GeneticArchitectures: OneLocusDiploid, dominance

export SelectionModel, FrequencyDependentSelection
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

## Matrix layout

The canonical state for a single-locus diploid simulation is a
`3 × n_demes` matrix (rows index genotypes, columns index demes):

| Row | Genotype |
|-----|----------|
| 1   | A1A1     |
| 2   | A1A2     |
| 3   | A2A2     |

This layout is cache-friendly for Julia's column-major memory order: iterating
over demes accesses three contiguous Float64 values per column.

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

end # module SelectionModels
