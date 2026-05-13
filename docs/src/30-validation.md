# Validation

This page documents validation work comparing HybridZones.jl against
historical reference implementations. The package's correctness is
anchored to specific numerical outputs from those reference
implementations, with documented tolerances and explicit notes about
known algorithmic differences.

## Reference implementations

Three reference sources are used for validation:

1. **Mallet's Pascal simulators** (1988-2000), particularly WC1SEDO.PAS
   for single-locus simulation. Run under DOSBox with Turbo Pascal 7.0
   on canonical parameter sets.

2. **Rosser et al. 2014 R supplementary code**, which reimplemented
   the multilocus framework in R for the Heliconius hybrid zone
   analysis published in Evolution.

3. **Analytical predictions** from the underlying theory (Mallet &
   Barton 1989, Kruuk et al. 1999), where closed-form expressions
   exist for cline shape, width, and equilibrium properties.

Reference data files are stored in `test/reference/` with provenance
documentation indicating which reference produced each file and
under which parameters.

## Pascal comparison: WC1SEDO single-locus

The single-locus simulation has been compared against Mallet's
WC1SEDO Pascal program for one canonical parameter set:

- 101 demes (Pascal: indices 0-101 = 102 demes; Julia: indices
  1-101 = 101 demes)
- 200 generations
- Selection coefficient s = 0.3
- Migration grain = 60, giving effective σ² ≈ 30
- Codominant (semi-dominant) genetic architecture
- Standard secondary contact initial condition

### Findings

The comparison revealed a systematic algorithmic difference. Julia
produces clines approximately 1.5× narrower than Pascal:

| Quantity | Julia | Pascal | Ratio |
|----------|-------|--------|-------|
| Width (10-90% of frequency range) | 31.0 demes | 46.0 demes | 1.48 |
| Maximum spatial gradient at center | 0.0345 / deme | 0.0230 / deme | 0.67 |
| Cline center deme | 51.0 (exactly) | ~50.0 (interpolated) | — |

Element-wise allele frequency differences across the transect are
small at the boundaries (~1%) and largest in the cline transition
zone (~9% maximum). This structured pattern indicates a real
algorithmic difference rather than accumulated floating-point drift.

### Diagnosis

Source code analysis identified the difference. Pascal's
`MigrationSelection` procedure (lines 741-771 of WC1SEDO.PAS) combines
migration and selection in a single pass per deme:

1. Compute incoming genotype frequencies via convolution with the
   migration kernel (each deme receives weighted contributions from
   neighbors)
2. Apply phenotype-aware fitness to each genotype
3. Normalize by mean fitness

This fusion effectively averages selection pressure across the
migration neighborhood. The selection coefficient at any given deme
acts on a smoothed mixture of local and neighboring frequencies
rather than on local frequencies alone.

Julia's implementation separates these operations. The `select!`
function applies fitness at each deme based on local frequencies
only; `migrate!` then redistributes post-selection frequencies via
the kernel. Selection is therefore "pointwise" — its effect at each
deme is based purely on what's at that deme.

The two approaches are not mathematically equivalent. In the limit
of weak selection or short generation timesteps they would converge,
but for typical Heliconius parameters (s ≈ 0.1-0.6) the difference
is substantial and produces the observed 1.5× width discrepancy.

### Biological interpretation

The separated approach is also more biologically accurate for the
warning-color systems the framework was originally developed to
model. In real hybrid zones, the two phenomena are temporally and
spatially distinct:

- Adult butterflies disperse from their natal sites (migration)
- They subsequently experience predation in their new location,
  where local predators have learned to recognize the locally
  common warning pattern (selection)

Migration is therefore a property of butterflies; selection is a
property of the local predator community. Pascal's fusion conflates
these into a single weighted-average operation, which doesn't
correspond to a discrete biological process.

The historical reason for the fusion was almost certainly
computational efficiency. On 1988-era hardware (PC AT class machines
with optional Intel 8087 coprocessors), separating the operations
into two passes per generation would have approximately doubled the
floating-point operation count. Mallet's choice traded mathematical
fidelity for runtime performance, accepting the approximation
inherent in the fused form.

### Design decision

HybridZones.jl uses the separated form. This is a deliberate
modernization that prioritizes biological accuracy over historical
fidelity to the Pascal reference. Modern hardware makes the
performance argument irrelevant, and the separated form aligns with
the underlying biology.

Users requiring exact reproduction of Pascal numerical output for
purposes such as reproducing specific published analyses should be
aware that this package will produce different numbers than the
Pascal originals at the same parameter values, with the difference
being approximately a 1.5× factor in cline width.

If exact Pascal reproducibility becomes a recurring requirement,
a future extension could provide a "fused" composite model that
implements Pascal's original `MigrationSelection` algorithm as an
explicit alternative to the default separated form. This is not
currently planned but the architecture supports such an extension.

### Update: Hardy-Weinberg restoration

A second algorithmic difference was subsequently identified. The previous
HybridZones.jl implementation applied selection and migration each generation
but did not include an explicit mating step. Migration mixing across demes
introduces Wahlund-effect heterozygote deficits relative to Hardy-Weinberg
expectations, and without an explicit mating step these deficits accumulate
over generations, artificially sharpening the cline.

Pascal's `WC1SEDO.PAS` includes an implicit HW restoration at each generation
as part of its state update. The `MatingModels` submodule (added in a
subsequent PR) adds `RandomMating` as the initial concrete mating type, which
restores HW proportions within each deme each generation. The default simulate
lifecycle is now `:mate_migrate_select`, matching Pascal's per-generation
ordering. This partially closes the cline-width gap relative to Pascal; the
remaining discrepancy is attributable to the fused vs. separated
selection-migration difference described above.

## Tolerance for regression tests

Regression tests against the Pascal reference data use the cline
*shape* as a qualitative validation rather than element-wise
numerical equality. Specific assertions:

- The cline is monotonic (allele frequency increases or decreases
  consistently across the transect)
- The cline is approximately symmetric around its center
- The cline center is at the geometric midpoint of the transect
  to within a few demes
- The cline width is consistent across runs with identical parameters
  to floating-point precision

Element-wise comparison against Pascal output is documented but
not used as a regression test, given the known algorithmic
difference.

## Future validation work

Planned validation activities:

- Comparison against Rosser et al. (2014) R supplementary code,
  which uses a different implementation but should produce results
  equivalent to Julia's separated form (subject to verification)
- Analytical comparison against Kruuk et al. (1999) closed-form
  expressions for cline shape under simple parameter regimes
- Three-locus validation once `MultiLocusUnlinked` and recombination
  are implemented

Each validation will be added to this page as it's completed.
