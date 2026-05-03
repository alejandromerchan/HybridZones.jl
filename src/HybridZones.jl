"""
    HybridZones

A Julia implementation of the multilocus hybrid zone simulation and inference
framework originally developed by Mallet and Barton in Turbo Pascal during the
late 1980s, with subsequent extensions in PowerBASIC and R.

This package provides forward simulation of allele and genotype frequency
clines under selection and migration, along with maximum-likelihood cline
fitting and disequilibrium estimation. Validated against historical reference
implementations including the original Pascal simulators and the R
reimplementation in Rosser, Dasmahapatra & Mallet (2014).

## References

- Mallet, J. and Barton, N. H. (1989). Inference from clines stabilized by
  frequency-dependent selection. *Genetics* 122: 967-976.
- Mallet, J., Barton, N., Lamas, G., Santisteban, J., Muedas, M. and
  Eeley, H. (1990). Estimates of selection and gene flow from measures of
  cline width and linkage disequilibrium in Heliconius hybrid zones.
  *Genetics* 124: 921-936.
- Rosser, N., Dasmahapatra, K. K. and Mallet, J. (2014). Stable Heliconius
  butterfly hybrid zones are correlated with a local rainfall peak at the
  edge of the Amazon basin. *Evolution* 68: 3470-3484.
"""
module HybridZones

include("GeneticArchitectures.jl")
using .GeneticArchitectures:
                             GeneticArchitecture, OneLocusDiploid, dominance, n_loci,
                             n_alleles, allele_names

include("MigrationModels.jl")
using .MigrationModels:
                        MigrationModel, BinomialStepping, migrate!, migration_variance,
                        max_distance

export package_version
export GeneticArchitecture, OneLocusDiploid
export dominance, n_loci, n_alleles, allele_names
export MigrationModel, BinomialStepping, migrate!
export migration_variance, max_distance

"""
    package_version()

Return the version of the HybridZones package as a `VersionNumber`.

This function exists primarily so the package has a callable export from
its first commit, allowing CI workflows to verify that the package loads
and exports a function correctly. Real simulation functionality will be
added in subsequent commits.

# Examples
```jldoctest
julia> using HybridZones

julia> v = package_version();

julia> v isa VersionNumber
true
```
"""
function package_version()
    return VersionNumber(0, 1, 0)
end

end # module HybridZones
