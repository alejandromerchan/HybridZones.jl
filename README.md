# HybridZones.jl

Forward simulation and inference for multilocus hybrid zones,
implementing the Mallet & Barton (1989) framework in Julia.

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://alejandromerchan.github.io/HybridZones.jl/stable)
[![Development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://alejandromerchan.github.io/HybridZones.jl/dev)
[![Test workflow status](https://github.com/alejandromerchan/HybridZones.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/alejandromerchan/HybridZones.jl/actions/workflows/Test.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/alejandromerchan/HybridZones.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/alejandromerchan/HybridZones.jl)
[![Lint workflow Status](https://github.com/alejandromerchan/HybridZones.jl/actions/workflows/Lint.yml/badge.svg?branch=main)](https://github.com/alejandromerchan/HybridZones.jl/actions/workflows/Lint.yml?query=branch%3Amain)
[![Docs workflow Status](https://github.com/alejandromerchan/HybridZones.jl/actions/workflows/Docs.yml/badge.svg?branch=main)](https://github.com/alejandromerchan/HybridZones.jl/actions/workflows/Docs.yml?query=branch%3Amain)
[![DOI](https://zenodo.org/badge/DOI/FIXME)](https://doi.org/FIXME)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](CODE_OF_CONDUCT.md)
[![BestieTemplate](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JuliaBesties/BestieTemplate.jl/main/docs/src/assets/badge.json)](https://github.com/JuliaBesties/BestieTemplate.jl)

## Project status

HybridZones.jl is under active development. The current state:

- Single-locus simulation under frequency-dependent selection,
  random mating (Hardy-Weinberg restoration), and binomial
  stepping-stone migration is implemented and validated against
  analytical expectations
- Mass conservation, equilibrium convergence, and cline symmetry
  hold to floating-point precision
- Two complementary selection models are implemented:
  `FrequencyDependentSelection` (strict codominance) and
  `SemiDominantFrequencyDependentSelection` (Pascal semi-dominance,
  for exact reproducibility of Mallet's WC1SEDO outputs)
- Validation against Mallet's WC1SEDO Pascal simulator is
  substantially complete: the 1.37× cline-width difference between
  the two implementations is fully diagnosed — it is explained
  entirely by the selection model distinction (strict codominance
  vs. Pascal semi-dominance), not by migration or lifecycle differences
- Multi-locus types, cline fitting, and disequilibrium estimation
  are planned but not yet implemented
- The package is not yet registered with the Julia General registry

See the [architecture document](docs/src/10-architecture.md) for the
full development roadmap and design intent.

## Background

Hybrid zones are regions where genetically distinct populations meet
and exchange migrants, producing characteristic clines in allele
frequencies. The theoretical framework for understanding these clines
was developed by Jim Mallet, Nick Barton, and collaborators starting
in the late 1980s, originally implemented in Turbo Pascal and applied
extensively to *Heliconius* butterfly hybrid zones in the Neotropics.

HybridZones.jl modernizes that framework for current Julia-based
scientific computing. The goal is to preserve the methodological
substance of the original work — including its specific selection
models, migration kernels, and likelihood-based inference — while
making it composable, extensible, and substantially faster than either
the original Pascal or the more recent R reimplementations. Validation
against historical reference implementations is a first-class concern:
the package's correctness is anchored to specific numerical outputs
from Mallet's original Pascal programs and Rosser et al. (2014)'s R
implementation.

The package complements rather than replaces existing tools.
Researchers needing genomic-cline inference on anonymous SNP data
should use bgchm or HZAR; those needing general-purpose forward
simulation should use SLiM; this package occupies the specific niche
of mechanistic forward simulation under classical multilocus cline
genetics with first-class validation against the historical literature.

## Quick start

```julia
using Pkg
Pkg.add(url="https://github.com/alejandromerchan/HybridZones.jl")

using HybridZones

# Build a simulation: single-locus codominant model,
# random mating, frequency-dependent selection, stepping-stone migration
arch = OneLocusDiploid(dominance = :codominant)
mat  = RandomMating()
sel  = FrequencyDependentSelection(s = 0.3)
mig  = BinomialStepping(30.0)

# Standard secondary contact: 101 demes, contact at the center
initial = secondary_contact(arch)

# Simulate 200 generations (enough to reach equilibrium)
trajectory = simulate(arch, mat, sel, mig, initial; n_generations = 200)

# Extract allele frequencies from genotype state
cline = allele_frequencies(trajectory[:, :, end], arch)
```

See the [Getting Started tutorial](https://alejandromerchan.github.io/HybridZones.jl/dev/01-getting-started/)
for a complete walkthrough including visualization and parameter
exploration.

## Documentation

Full documentation is at
[alejandromerchan.github.io/HybridZones.jl](https://alejandromerchan.github.io/HybridZones.jl/dev/).

Key pages:

- [Getting Started](https://alejandromerchan.github.io/HybridZones.jl/dev/01-getting-started/) — your first simulation
- [Architecture](https://alejandromerchan.github.io/HybridZones.jl/dev/10-architecture/) — design intent and development roadmap
- [Performance](https://alejandromerchan.github.io/HybridZones.jl/dev/20-performance/) — benchmarking and optimization notes
- [Validation](https://alejandromerchan.github.io/HybridZones.jl/dev/30-validation/) — comparison against Pascal reference implementations
- [API Reference](https://alejandromerchan.github.io/HybridZones.jl/dev/95-reference/) — full function and type documentation

## Contributing

Contributions are welcome. See [CONTRIBUTING.md](docs/src/90-contributing.md)
for development setup, coding conventions, and the pull request
workflow.

The package follows a small-PR, validation-first development style:
each change is reviewed against the package's regression test suite,
which includes specific numerical references from the historical
Pascal implementations.

## Citation

If you use HybridZones.jl in published research, please cite both
the package itself and the foundational papers it builds on. See
[CITATION.cff](CITATION.cff) for the full bibliography in standard
format.

The most important reference for the simulation framework is:

> Mallet, J. and Barton, N. H. (1989). Inference from clines
> stabilized by frequency-dependent selection. *Genetics* 122:
> 967-976.

## License

MIT. See [LICENSE](LICENSE).
