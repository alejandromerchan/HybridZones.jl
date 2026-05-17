# [Reference](@id reference)

## Contents

```@contents
Pages = ["95-reference.md"]
```

## Index

```@index
Pages = ["95-reference.md"]
```

```@autodocs
Modules = [HybridZones]
```

## Genetic architectures

```@docs
HybridZones.GeneticArchitectures.GeneticArchitecture
HybridZones.GeneticArchitectures.OneLocusDiploid
HybridZones.GeneticArchitectures.dominance
HybridZones.GeneticArchitectures.n_loci
HybridZones.GeneticArchitectures.n_alleles
HybridZones.GeneticArchitectures.allele_names
HybridZones.GeneticArchitectures.n_genotypes
```

## Selection models

```@docs
HybridZones.SelectionModels.SelectionModel
HybridZones.SelectionModels.FrequencyDependentSelection
HybridZones.SelectionModels.select!
HybridZones.SelectionModels.selection_coefficient
```

## Semi-dominant selection

```@docs
HybridZones.SelectionModels.SemiDominantFrequencyDependentSelection
```

## Migration models

```@docs
HybridZones.MigrationModels.MigrationModel
HybridZones.MigrationModels.BinomialStepping
HybridZones.MigrationModels.migrate!
HybridZones.MigrationModels.migration_variance
HybridZones.MigrationModels.max_distance
HybridZones.MigrationModels.kernel
```

## Mating models

```@docs
HybridZones.MatingModels.MatingModel
HybridZones.MatingModels.RandomMating
HybridZones.MatingModels.mate!
```

## Simulation

```@docs
HybridZones.Simulation.simulate
HybridZones.Simulation.secondary_contact
HybridZones.Simulation.allele_frequencies
```
