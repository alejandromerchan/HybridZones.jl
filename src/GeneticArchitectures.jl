module GeneticArchitectures

export GeneticArchitecture, OneLocusDiploid
export dominance, n_loci, n_alleles, allele_names

"""
    GeneticArchitecture

Abstract supertype for genetic architecture specifications. Concrete subtypes
carry locus count, ploidy, allele identities, and dominance relationships used
by selection and migration models to dispatch the correct update rules.
"""
abstract type GeneticArchitecture end

const _VALID_DOMINANCE = (:codominant, :dominant_first, :recessive_first)

"""
    OneLocusDiploid <: GeneticArchitecture

A single autosomal diploid locus with two alleles.

`dominance` specifies how alleles combine in heterozygotes:
- `:codominant`: both alleles contribute equally (additive)
- `:dominant_first`: the first allele is fully expressed; the second is masked
- `:recessive_first`: the first allele is masked; the second is fully expressed

`allele_names` provides human-readable labels for the two alleles, used in
display and output. Defaults to `("A1", "A2")`.

# Examples
```jldoctest
julia> using HybridZones

julia> arch = OneLocusDiploid();

julia> arch isa GeneticArchitecture
true

julia> dominance(arch)
:codominant

julia> allele_names(arch)
("A1", "A2")
```
"""
struct OneLocusDiploid <: GeneticArchitecture
    dominance::Symbol
    allele_names::NTuple{2, String}
end

function OneLocusDiploid(;
        dominance::Symbol = :codominant,
        allele_names::NTuple{2, String} = ("A1", "A2")
)
    dominance ∈ _VALID_DOMINANCE ||
        throw(ArgumentError("dominance must be one of $_VALID_DOMINANCE, got :$dominance"))
    return OneLocusDiploid(dominance, allele_names)
end

"""
    dominance(arch::OneLocusDiploid) -> Symbol

Return the dominance mode of `arch`: one of `:codominant`, `:dominant_first`,
or `:recessive_first`.
"""
dominance(arch::OneLocusDiploid) = arch.dominance

"""
    n_loci(arch::OneLocusDiploid) -> Int

Return the number of loci. Always 1 for `OneLocusDiploid`.
"""
n_loci(::OneLocusDiploid) = 1

"""
    n_alleles(arch::OneLocusDiploid) -> Int

Return the number of alleles per locus. Always 2 for `OneLocusDiploid`.
"""
n_alleles(::OneLocusDiploid) = 2

"""
    allele_names(arch::OneLocusDiploid) -> NTuple{2,String}

Return the two allele labels.
"""
allele_names(arch::OneLocusDiploid) = arch.allele_names

end # module GeneticArchitectures
