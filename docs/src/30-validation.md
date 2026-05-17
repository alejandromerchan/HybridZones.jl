# Validation

This page documents validation work comparing HybridZones.jl against
historical reference implementations. The package's correctness is
anchored to specific numerical outputs from those reference
implementations, with documented tolerances and explicit notes about
known algorithmic differences.

## Reference implementations

Three reference sources are used for validation:

1. **Mallet's Pascal simulators** (1988–2000), particularly WC1SEDO.PAS
   for single-locus simulation. Run under DOSBox with Turbo Pascal 7.0
   on canonical parameter sets.

2. **Rosser et al. 2014 R supplementary code**, which reimplemented
   the multilocus framework in R for the *Heliconius* hybrid zone
   analysis published in Evolution.

3. **Analytical predictions** from the underlying theory (Mallet &
   Barton 1989, Kruuk et al. 1999), where closed-form expressions
   exist for cline shape, width, and equilibrium properties.

Reference data files are stored in `test/reference/` with provenance
documentation indicating which reference produced each file and
under which parameters.

## Pascal comparison: WC1SEDO single-locus

The single-locus simulation has been compared against Mallet's
WC1SEDO Pascal program for a canonical parameter set:

- 101 demes (Pascal: indices 0–101 = 102 demes; Julia: indices
  1–101 = 101 demes)
- 200 generations
- Selection coefficient s = 0.3
- Migration grain = 60, giving effective σ² ≈ 30
- Codominant (semi-dominant in Pascal) genetic architecture
- Standard secondary contact initial condition

### Findings

HybridZones with the default `:mate_migrate_select` lifecycle and
`FrequencyDependentSelection` with `:codominant` dominance produces
a cline approximately 1.37× narrower than Pascal:

| Quantity | HybridZones (codominant) | Pascal (semi-dominant) | Ratio |
|----------|--------------------------|------------------------|-------|
| Width (10–90%) | 33.65 demes | 46.02 demes | 1.368 |
| Max element-wise \|Δp\| | — | — | 0.071 |

The structured pattern of frequency differences (small at boundaries,
peaking at ~0.071 in the cline centre) indicates a real algorithmic
difference rather than accumulated floating-point drift.

### Diagnosis

Two algorithmic candidates were investigated and one was ruled out:

**Candidate 1: Fused vs separated migration–selection operations.**
Pascal's `MigrationSelection` procedure (lines 741–771 of WC1SEDO.PAS)
combines migration and selection in a single pass per deme, effectively
averaging selection pressure across the migration neighbourhood.
HybridZones separates the operations: `migrate!` runs first, then
`select!` acts on the post-migration frequencies. These two formulations
are not mathematically equivalent; at typical *Heliconius* parameters
they were expected to diverge.

A controlled experiment tested this by implementing both fused
(Pascal-style) and separated (HybridZones-style) variants in Julia
and running both at the same lifecycle position. The two variants
produced bit-identical results: maximum element-wise frequency
difference across the transect was 0.0000. **Pascal's fusion is a
structural code-organisation choice, not a mathematical distinction.**
It produces no numerical difference compared to separated operations
at the same lifecycle position.

**Candidate 2: The selection model itself.**
Pascal's WC1SEDO uses semi-dominant frequency-dependent selection.
The heterozygote A1A2's "perceived phenotype frequency" — the
frequency of the phenotype class that receives the heterozygote's
fitness — includes contributions from both homozygote classes via
a `SimilarPhenotypes` mixing step (lines 720–735 of WC1SEDO.PAS).
HybridZones' `FrequencyDependentSelection` with `:codominant`
dominance uses strict codominance: A1A2 is treated as its own
distinct phenotype with no cross-similarity to either homozygote.

This selection-model difference produces the full observed width
discrepancy.

#### Isolation experiment: factorial 2×2 design

A fully crossed factorial experiment varied selection model (strict
codominant vs Pascal semi-dominant) against operation structure
(fused vs separated), holding all other parameters fixed:

| | HybridZones codominant | Pascal semi-dominant | Row effect |
|---|---|---|---|
| Separated operations | 33.65 demes | 46.02 demes | — |
| Fused operations | 33.65 demes | 46.02 demes | — |
| **Column effect** | — | — | **1.368×** |

The selection-model effect (column direction) is consistent at 1.368×
regardless of operation structure. The fusion effect (row direction)
is exactly zero. The two effects compose multiplicatively, but since
one factor is zero, only the selection model contributes.

#### Historical note: Hardy-Weinberg restoration

An earlier hypothesis attributed part of the cline-width gap to
missing Hardy-Weinberg (HW) restoration. The previous HybridZones
implementation applied selection and migration each generation but
did not include an explicit mating step. Migration mixing across
demes introduces Wahlund-effect heterozygote deficits, and without
mating those deficits accumulate, artificially sharpening the cline.

Adding `RandomMating` with HW restoration (the `MatingModels`
submodule, PR #25) did change cline width — from approximately
31 demes to 33.65 demes, a shift of roughly 8–9%. This is a real
effect but small: it explains only about 6% of the total 1.37×
discrepancy between HybridZones and Pascal. The dominant difference
is the selection model, not missing mating.

The mating step is nonetheless correct to include. Pascal's
WC1SEDO implicitly restores HW proportions each generation as part
of its state update, so `RandomMating` is needed for the simulation
to correctly represent the Pascal life-cycle convention.

### Biological interpretation

Both selection models are biologically valid for warning-color
systems. The distinction matters because heterozygotes in real hybrid
zones can have two different relationships to the parental phenotypes:

1. **Strict codominance.** Heterozygotes display a visually distinct
   intermediate pattern that predators recognise as its own warning-
   color class. The predator community develops a separate learned
   avoidance response specific to the heterozygote pattern. This is
   appropriate when the heterozygote is visually different from both
   parental types — as occurs for some *Heliconius* colour-pattern
   loci where dominance modifiers are weak or absent.

2. **Semi-dominance.** Heterozygotes share visual similarity with
   both parental phenotypes and are partially recognised as either.
   Predators that have learned one parental pattern partly generalise
   to the heterozygote. This is appropriate when the heterozygote is
   visually intermediate and the predator generalises across similar
   patterns — as occurs for many *Heliconius* wing-pattern elements
   where the heterozygote appears as a blended intermediate.

Mallet's Pascal implementation was designed for the second case,
consistent with the empirical biology of *H. erato* and *H. melpomene*
in the Peruvian systems his framework was originally developed for.

### Design decision

HybridZones implements two complementary selection models:

- `FrequencyDependentSelection` with `:codominant` dominance implements *strict*
  codominance: each of the three diploid genotypes (A1A1, A1A2, A2A2) is treated
  as a distinct phenotype with no cross-similarity. This is appropriate when
  heterozygotes have phenotypically distinct intermediate patterns that predators
  treat as a separate class.

- `SemiDominantFrequencyDependentSelection` implements Pascal's semi-dominance
  mixing: the heterozygote's perceived phenotype frequency includes weighted
  contributions from both homozygote classes (0.5 weighting each way). This
  model reproduces exact Pascal cline widths and is appropriate when
  heterozygotes share visual similarity with both parental phenotypes.

For users requiring exact Pascal reproducibility — for example, to reproduce
specific published analyses or to compare directly against the original
Mallet 1990 *erato* fit values — use `SemiDominantFrequencyDependentSelection`.
See the [Architecture](10-architecture.md) document for the full model
descriptions and biological interpretation.

## Lifecycle ordering effects

Independent of the selection model, lifecycle ordering produces
measurable differences in cline width. With `FrequencyDependentSelection`
and `:codominant` dominance (s = 0.3, σ² = 30, 200 generations):

| Ordering | Width (10–90%) |
|----------|---------------|
| `:mate_migrate_select` (default) | 33.65 demes |
| `:mate_select_migrate` | 38.20 demes |
| `:select_migrate_mate` | 38.20 demes |

`:mate_migrate_select` produces the narrowest cline because selection
acts on the post-migration state, which has Wahlund-effect heterozygote
deficits from mixing across demes. Those deficits amplify the rare-
morph penalty under frequency-dependent selection, driving stronger
divergence. The other two orderings apply selection before or just
after mating to a state closer to Hardy-Weinberg, producing weaker
effective selection and broader clines.

`:mate_select_migrate` and `:select_migrate_mate` produce identical
allele frequency trajectories (maximum element-wise difference 0.0
at floating-point precision) despite producing different genotype
frequency trajectories. This is a mathematical consequence of
`RandomMating` being an idempotent allele-frequency-conserving
projection: mating maps any genotype state to its Hardy-Weinberg
projection without changing allele frequencies, so transposing
mating and a selection–migration sequence leaves allele frequencies
unchanged even though genotype frequencies differ. Users who care
about within-deme genotype dynamics (e.g., heterozygote deficits)
should use the full genotype trajectory to distinguish these two
orderings; users who only track allele frequencies can treat them
as equivalent.

## Dominance and bistability

Under `:dominant_first` or `:recessive_first` dominance,
`FrequencyDependentSelection` has bistable dynamics in single-locus
systems. The unstable interior equilibrium sits at allele frequency
p ≈ 0.707 (where the recessive phenotype frequency q² ≈ 0.5,
making the rare-morph penalty on the recessive allele switch from
advantageous to disadvantageous).

With the standard `secondary_contact` initial condition, the contact
deme starts at allele frequency p = 0.5 for the A1 allele. For
`:recessive_first` dominance, this is *below* the unstable threshold
for the recessive (A1) allele. Selection then drives the recessive
allele toward local extinction, and migration propagates the sweep
across the transect. No stable cline forms.

This is biologically correct: single-locus dominance under
frequency-dependent selection cannot maintain a stable polymorphism
without additional forces. Real warning-color hybrid zones with
dominant patterns are stabilized by multi-locus epistasis across
colour-pattern loci, where heterozygotes at one locus interact with
genotypes at other loci to produce phenotypically distinct hybrid
combinations that experience distinct selection pressures.

Multi-locus extensions (`MultiLocusUnlinked` and
`MultiLocusRecombining`) are planned to support this case.

## Tolerance for regression tests

Regression tests against the Pascal reference data use qualitative
cline properties rather than element-wise numerical equality, given
the documented selection-model difference:

- The cline is monotonic (allele frequency increases consistently
  across the transect)
- The cline is symmetric around its centre
- The cline centre is at the geometric midpoint of the transect
  to within approximately 1 deme
- The cline width is reproducible to floating-point precision
  across runs with identical parameters

Element-wise comparison against Pascal output is documented in
`test/reference/` but is not used as a regression test, given the
known and understood selection-model difference.

## Future validation work

Planned validation activities:

- Comparison against Rosser et al. (2014) R supplementary code,
  which uses a different implementation but should produce results
  equivalent to HybridZones' separated form at the same selection
  model (subject to verification)
- Analytical comparison against Kruuk et al. (1999) closed-form
  expressions for cline shape under simple parameter regimes
- Three-locus validation once `MultiLocusUnlinked` is implemented,
  comparing against Pascal's WC3ARBCD multi-locus program

Each validation will be added to this page as it is completed.
