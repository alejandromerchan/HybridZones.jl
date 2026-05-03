# Architecture

This document describes the architectural intent of HybridZones.jl: what the
package is for, what it is not, and how its internal structure is organized to
support those goals. It is meant to evolve alongside the code, not to specify
the code in advance. When design decisions get made in PRs, they should be
reflected here.

## Purpose

HybridZones.jl provides forward simulation and inference for spatially
structured populations exchanging migrants across hybrid zones, with a
particular focus on the multilocus framework developed by Mallet and Barton in
the late 1980s and extended through the 1990s and 2000s. The package targets
researchers studying speciation, warning-color evolution, ecological clines,
and other systems where the spatial dynamics of allele frequencies, genotype
frequencies, and linkage disequilibrium across a transition zone are the
primary objects of interest.

The package is intended to complement, not replace, existing tools. R packages
like HZAR and bgchm handle data-fitting of cline shapes and genomic-cline
inference well. SLiM provides general-purpose forward simulation in a
domain-specific scripting language. HybridZones.jl occupies the niche of
mechanistic forward simulation under specific genetic architectures, with
inference machinery that exposes underlying biological parameters (selection
coefficients, dispersal variance, dominance, recombination) directly.

## Scope

HybridZones.jl will provide:

- Forward simulation of genotype frequency clines under composable selection
  and migration models
- Configurable genetic architectures with explicit dominance and recombination
  specifications
- Maximum-likelihood and likelihood-profile inference for cline parameters
  from empirical transect data
- Two-locus and multi-locus linkage disequilibrium estimation accounting for
  dominance
- Validation against historical reference implementations (Pascal simulators,
  PowerBASIC disequilibrium estimator, R reimplementations) on canonical
  datasets

The package will not provide:

- Genomic-cline inference on anonymous SNP data (use bgchm)
- Phenomenological cline-shape fitting without genetic mechanism (use HZAR)
- General-purpose forward simulation with arbitrary population structure
  (use SLiM)
- Coalescent or backward-in-time simulation (use msprime)

These are deliberate exclusions. Tools that try to do everything tend to do
nothing especially well, and the existing tools listed handle their respective
niches better than a from-scratch reimplementation would.

## Out of scope for v1

The following modeling capabilities are recognized as valuable
extensions but are deliberately deferred from the initial scope.
Documenting them here helps users with such needs select appropriate
alternative tools and helps future contributors understand which
extensions would be welcome.

**Population density dynamics.** The current framework assumes constant
carrying capacity at every deme — populations are always full, and
migration redistributes frequencies rather than absolute numbers of
individuals. Scenarios involving range expansion into empty habitat,
demographic bottlenecks, founder effects, or density-dependent
dynamics require a different state representation that tracks
population density alongside genotype frequencies. The reaction-
diffusion framework of Sasaki, Kawaguchi & Yoshimori (2002) is the
natural mathematical reference for such extensions.

**Two-dimensional spatial structure.** The framework treats space as
a one-dimensional transect of demes. Real hybrid zones are
two-dimensional phenomena, and 2D systems can produce stable spatial
mosaics rather than smooth clines (Sasaki et al. 2002). Adding 2D
support would require new migration models and visualizations but
the existing genetic architecture and selection model abstractions
would carry over largely unchanged. The architecture document's
"State representation as plain numerical arrays" principle leaves
this extension path open.

**Stochastic dynamics and genetic drift.** The current framework is
fully deterministic. Real populations experience genetic drift,
particularly at small effective population sizes, and stochastic
extinction-recolonization dynamics in metapopulations. A stochastic
extension would track individuals or use diffusion approximations
rather than exact frequency dynamics. Kruuk, Baird, Gale & Barton
(1999) discuss the role of drift in their multilocus framework.

**Time-varying parameters.** Selection coefficients, migration rates,
and dominance are treated as constant across generations. Real
systems may have temporally varying selection (predator learning,
seasonal effects, climate change) or temporally varying migration
(barrier formation/dissolution). Time-varying selection could be
implemented as a wrapper around existing selection models that
modulates parameters per generation; time-varying migration similarly.

**Continuous-space models.** Demes are discrete in this framework,
matching the original Pascal implementations. Some applications
benefit from continuous-space formulations using partial differential
equations (cline width as a continuous function of position rather
than a per-deme array). The discrete-deme approach is simpler to
implement, validate, and reason about, but represents an
approximation when dispersal scales are large relative to deme
spacing.

Users needing any of these capabilities should consider
complementary tools: SLiM for general-purpose stochastic forward
simulation, msprime for backward-in-time coalescent simulation, or
custom implementations of the Sasaki et al. reaction-diffusion
framework for population density dynamics.

## Historical lineage

The framework HybridZones.jl modernizes was originally developed in Turbo
Pascal by Jim Mallet and collaborators starting in the late 1980s, with the
linkage-disequilibrium estimation extended in PowerBASIC through the 1990s
and 2000s. The package aims to preserve the methodological substance of that
work while moving it to a modern language and ecosystem.

Key references:

- Mallet, J. and Barton, N. H. (1989). Inference from clines stabilized by
  frequency-dependent selection. *Genetics* 122: 967-976.
- Mallet, J., Barton, N., Lamas, G., Santisteban, J., Muedas, M. and Eeley, H.
  (1990). Estimates of selection and gene flow from measures of cline width
  and linkage disequilibrium in *Heliconius* hybrid zones. *Genetics* 124:
  921-936.
- Rosser, N., Dasmahapatra, K. K. and Mallet, J. (2014). Stable *Heliconius*
  butterfly hybrid zones are correlated with a local rainfall peak at the
  edge of the Amazon basin. *Evolution* 68: 3470-3484.

The original Pascal sources, PowerBASIC sources, and Excel validation
spreadsheet from the 2003 disequilibrium-estimator debugging session are
preserved in the project's reference materials.

The framework has been extended in several theoretical directions
that are referenced throughout this document but not implemented in
the initial package version:

- Kruuk, L. E. B., Baird, S. J. E., Gale, K. S. and Barton, N. H. (1999).
  A comparison of multilocus clines maintained by environmental
  adaptation or by selection against hybrids. *Genetics* 153: 1959-1971.

- Sasaki, A., Kawaguchi, I. and Yoshimori, A. (2002). Spatial mosaic
  and interfacial dynamics in a Müllerian mimicry system. *Theoretical
  Population Biology* 61: 49-71.

## Architectural principles

The package is organized around a few principles that should be reflected in
all design decisions.

### Composition over hardcoding

The Pascal simulators are monolithic programs where each combination of
selection regime, migration model, and genetic architecture is a separate
executable. HybridZones.jl factors these into independent components that
combine via multiple dispatch. A simulation is the composition of:

- A genetic architecture (number of loci, dominance per locus, recombination
  structure, genotype-to-phenotype map)
- A selection model (frequency-dependent, epistatic, ecological-gradient,
  underdominance, or compositions of these)
- A migration model (binomial stepping-stone, Gaussian, asymmetric, or others)
- A lifecycle specifying the order of operations per generation
- An initial condition (secondary contact, equilibrium, custom)

Users compose simulations from these components rather than choosing among a
fixed set of model variants.

### Lifecycle as data

The order of operations within a generation matters for cline shape and
disequilibrium predictions, particularly at high selection. Different
organisms have different natural orderings: selection-then-migration for
species with adult dispersal acting on locally-grown juveniles, migration-
then-selection for species with passively-dispersed propagules that settle
and experience local selection.

The lifecycle is treated as data (an ordered sequence of operation
specifications) rather than hardcoded into the simulator. This lets users
match the simulation to their organism's biology without modifying source
code.

### State representation as plain numerical arrays

The state of a simulation at any generation is represented as plain numerical
arrays of allele or genotype frequencies, indexed by deme and (where relevant)
genotype. This is for two reasons. First, plain numerical types allow
performance optimizations that complex object-oriented representations would
preclude. Second, it leaves the door open for GPU acceleration via CUDA.jl or
similar without rearchitecting, since GPU array types follow the same
interface.

Genotype-to-phenotype mappings, dominance specifications, and recombination
matrices are stored as small lookup tables computed once at simulation start,
not as functions called per-update.

### Validation as a first-class concern

The package's correctness is anchored to historical reference implementations
(Pascal, PowerBASIC, R) on canonical datasets. Each simulation algorithm and
each inference routine should have at least one regression test comparing its
output to a known reference. Reference data files live in
`test/reference/`. The validation philosophy is that a refactor that breaks
numerical agreement with the reference is treated as a regression unless
explicitly justified.

### Separable inference and simulation

The simulator and the inference machinery are independent modules, even
though they're often used together. Forward simulation produces theoretical
clines; inference fits theoretical clines to observed transect data via
likelihood. A user might do simulation without inference (theoretical
exploration) or inference with externally-generated theoretical clines. The
package supports both pathways.

## Module structure

The package is organized into submodules that may eventually become
independent packages. The boundaries are drawn so that each submodule could
serve users beyond the hybrid-zones use case.

### `GeneticArchitectures`

Defines the genetic structure of organisms being modeled: number of loci,
ploidy, allele identities, dominance relationships, recombination rates
between loci, and genotype-to-phenotype mapping. Concrete types include
`OneLocusDiploid`, with `MultiLocusUnlinked` and `MultiLocusRecombining`
planned for when a working multi-locus simulator requires them.

Dominance is represented as a `Symbol` with three values: `:codominant`
(additive), `:dominant_first` (first allele fully expressed in heterozygotes),
and `:recessive_first` (first allele masked in heterozygotes). The "first"
suffix refers to the first element of the `allele_names` tuple.

### `SelectionModels`

Defines abstract `SelectionModel` and concrete subtypes including
`FrequencyDependentSelection`, `EpistaticSelection`, `EcologicalGradient`,
`Underdominance`, and composition types for combining selection models.
Selection models act on population state and produce post-selection
frequencies.

### `MigrationModels`

Defines abstract `MigrationModel` and concrete subtypes including
`BinomialStepping` (the Mallet-Barton stepping-stone formulation),
`Gaussian`, and `AsymmetricStepping`. Migration models redistribute
frequencies across demes according to a kernel.

### `HybridZones` (top-level)

Composes architecture, selection, and migration models into a `simulate`
function that runs forward simulations under a specified lifecycle. Provides
convenience constructors for common scenarios (warning-color cline, ecological
cline, tension zone).

### `ClineFitting`

Maximum-likelihood inference of cline parameters from empirical transect
data. Implements the multinomial and dominant-phenotype likelihoods used by
the original CLINFIT family, with optimization via Optim.jl. Returns
parameter estimates with likelihood-based confidence intervals.

### `LinkageDisequilibrium`

Two-locus disequilibrium estimation accounting for dominance, with seven
estimators corresponding to the dominance combinations in DIS.BAS
(codominant-codominant through recessive-recessive plus codominant-haploid).
Provides maximum-likelihood point estimates and likelihood-profile support
limits, plus Monte Carlo permutation tests for significance.

## Validation strategy

Each numerical algorithm has at least one regression test against a reference
implementation. Reference sources, in order of preference:

1. Output from the original Pascal or PowerBASIC executables run in DOSBox on
   canonical datasets (Mallet 1990 *erato* data, Mallet et al. 2003 *Anartia*
   data).
2. Output from the R implementation of Rosser et al. 2014, where applicable.
3. Output from the Excel spreadsheet for two-locus disequilibrium estimators
   (the `DiseqEstAlejo.xls` file from the 2003 debugging session).
4. Hand-computed values for small, analytically tractable cases.

Test files in `test/` are organized by subsystem and named `test-*.jl` per
the auto-discover pattern. Reference data files live in `test/reference/`.

## Performance expectations

The package targets per-simulation costs that allow large parameter sweeps to
complete in seconds rather than hours. A three-locus *erato*-style simulation
(100 demes, 100 generations) should run in approximately 10 milliseconds per
parameter combination on modern hardware. A 1000-combination parameter sweep
should complete in 10 seconds or less when threaded across available cores.

This is roughly two orders of magnitude faster than the equivalent R
implementation in Rosser et al. 2014 supplementary code, which is the
expected ratio for the kind of nested-loop scalar-arithmetic workload these
simulators have. The package's value over the R implementation is partly this
performance and partly the architectural composability that enables analyses
the R code would not naturally support.

GPU acceleration (via CUDA.jl with optional `CuArray` state) is a possible
future extension but is not part of the initial scope. The current state
representation is GPU-friendly so that this extension remains feasible.

## Conventions

Internal conventions that should be followed in all PRs:

- Test files match `test-*.jl` and live in `test/`. Reference data files
  live in `test/reference/`.
- Conventional Commits format for commit messages: `feat:`, `fix:`,
  `docs:`, `refactor:`, `test:`, `ci:`, `chore:`.
- Doctests in docstrings where examples are useful and stable. The docs
  workflow runs them as part of CI.
- ExplicitImports compliance: no implicit imports, no stale explicit imports.
- JuliaFormatter compliance: code is auto-formatted on commit.

Architectural decisions of any consequence are documented in this file as
they're made. Trivial implementation choices stay in the source code.
