# Getting Started

Hybrid zones are regions where genetically distinct populations meet and
exchange migrants. Where selection favours locally common phenotypes — as in
warning-color mimicry systems — allele frequencies shift sharply from one
population to the other, producing a characteristic cline shape across the
contact zone. Measuring cline width and linkage disequilibrium from empirical
transect data makes it possible to estimate the underlying selection
coefficients and dispersal rates.

HybridZones.jl provides forward simulation of allele frequency clines under
user-specified selection and migration models, together with the inference
machinery to fit cline parameters to empirical transect data. A simulation is
built by composing a genetic architecture, a selection model, and a migration
model; the same components can be rearranged without touching source code to
explore different biological scenarios.

The package implements the multilocus cline framework developed by James
Mallet, Nicholas Barton, and collaborators starting in the late 1980s,
modernized in Julia. The framework has been applied extensively to *Heliconius*
butterfly hybrid zones in the Neotropics, where warning-color clines can be
measured precisely in the field.

## Installation

The package is not yet registered in the Julia General registry. Install
directly from the repository:

```julia
using Pkg
Pkg.add(url="https://github.com/alejandromerchan/HybridZones.jl")
```

Once registration is complete, the simpler form will work:

```julia
Pkg.add("HybridZones")
```

## Your first simulation

The following example builds a one-locus codominant cline under
frequency-dependent selection and stepping-stone migration, runs it to
equilibrium, and extracts the allele frequency profile. All of this fits
in a handful of lines.

```julia
using HybridZones

# Build the model components
arch = OneLocusDiploid(dominance = :codominant)
sel  = FrequencyDependentSelection(s = 0.3)
mig  = BinomialStepping(30.0)

# Set up the initial condition
initial = secondary_contact(arch)

# Run 200 generations to equilibrium
trajectory = simulate(arch, sel, mig, initial; n_generations = 200)

# Extract allele frequencies at the final generation
p = allele_frequencies(trajectory[:, :, end], arch)
```

`OneLocusDiploid` represents a single autosomal locus with two alleles and
three diploid genotypes (A1A1, A1A2, A2A2). The `:codominant` dominance mode
means each genotype maps to a distinct phenotype — there is no masking — which
is the appropriate model for *Heliconius* warning colors where heterozygotes
display an intermediate pattern.

`FrequencyDependentSelection(s = 0.3)` implements the Mallet–Barton fitness
function in which the rare phenotype suffers a fitness cost. The parameter `s`
is the selection coefficient: `s = 0.3` means a locally rare phenotype has
fitness `1 - s = 0.7` relative to the common phenotype. Typical values for
*Heliconius* warning-color hybrid zones lie in the range 0.1–0.6.

`BinomialStepping(30.0)` is stepping-stone migration with dispersal variance
σ² = 30 demes per generation. The dispersal kernel is binomial, matching the
Mallet & Barton 1989 Pascal implementation. Typical σ² values for *Heliconius*
lie in the range 5–100 demes.

`secondary_contact(arch)` returns a 101-deme transect in which the left half
is fixed for A1A1, the right half is fixed for A2A2, and the central deme
(deme 51) starts as all heterozygotes. This is the standard initial condition
for modelling secondary contact between two populations that have been
geographically isolated and are now meeting for the first time. An odd number
of demes is conventional so the contact point is a single central deme and the
two flanking regions are equal in length.

`simulate` runs the per-generation loop: selection then migration in the
default `:selection_first` order. The result is a three-dimensional array
indexed as `(genotype, deme, generation)`, where `trajectory[:, :, 1]` is
the initial condition and `trajectory[:, :, g+1]` is the state after
generation `g`. With `n_generations = 200`, the array has 201 time slices.

`allele_frequencies` converts the genotype frequency matrix into a more
interpretable quantity: the A2 allele frequency at each deme. For a
codominant locus this is `p_A2[j] = freq_A2A2[j] + 0.5 × freq_A1A2[j]`.
The result `p` is a 101-element vector. At equilibrium, `p[51] ≈ 0.5` by
symmetry, rising toward 1 on the right and falling toward 0 on the left.

## Visualizing the cline

With the equilibrium allele frequency vector `p` in hand, plotting the cline
takes a few lines with CairoMakie:

```julia
using CairoMakie

fig = Figure(size = (900, 500))
ax  = Axis(fig[1, 1],
    xlabel = "Deme",
    ylabel = "Allele frequency (A2)",
    title  = "Equilibrium cline (s=0.3, σ²=30)")
lines!(ax, 1:101, p, linewidth = 2)
fig
```

The resulting curve crosses 0.5 at deme 51 and is symmetric about the contact
point — a property of the codominant model combined with equal flanking
population sizes. The cline is sigmoidal, steep at the centre and approaching
the boundary values asymptotically. Cline width, defined as the inverse of the
maximum slope, is determined by the balance between selection pulling
frequencies apart and migration homogenizing them.

The plotting code above is shown as text and must be run interactively to
produce the figure. Configuring Documenter to render CairoMakie output in the
docs build would add dependency and configuration complexity not worth
introducing here.

## Exploring parameter dependence

Running the same simulation at several values of `s` illustrates how cline
width depends on selection strength. Stronger selection keeps local phenotypes
more common and narrows the cline; weaker selection allows migration to
homogenize frequencies over a broader spatial range.

```julia
selection_values = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
clines = Vector{Vector{Float64}}()
for s in selection_values
    sel_i  = FrequencyDependentSelection(s = s)
    final  = simulate(arch, sel_i, mig, initial;
                      n_generations = 200,
                      save_trajectory = false)
    push!(clines, allele_frequencies(final, arch))
end
```

Passing `save_trajectory = false` tells `simulate` to discard intermediate
generations and return only the final state as a matrix, rather than
accumulating a three-dimensional array. For parameter sweeps where only the
equilibrium cline shape matters, this reduces memory use substantially.
The final state is then passed directly to `allele_frequencies`, which accepts
either a 2D genotype matrix (returning a vector) or a 3D trajectory array
(returning a matrix of allele frequencies over time).

The relationship between selection strength and cline width follows
approximately width ∝ √(σ²/s) for the simplest cases. Doubling `s` therefore
narrows the cline by roughly a factor of √2; halving it broadens the cline by
the same factor.

To compare the six clines on one figure:

```julia
using CairoMakie

cmap   = cgrad(:viridis)
colors = [cmap[t] for t in LinRange(0.1, 0.9, length(selection_values))]

fig = Figure(size = (900, 500))
ax  = Axis(fig[1, 1],
    xlabel = "Deme",
    ylabel = "Allele frequency (A2)",
    title  = "Cline width vs. selection strength (σ²=30)")
for (i, (s_val, p_i)) in enumerate(zip(selection_values, clines))
    lines!(ax, 1:101, p_i,
        label     = "s = $(s_val)",
        color     = colors[i],
        linewidth = 2)
end
Legend(fig[1, 2], ax, "Selection")
fig
```

As above, run this block interactively; it is not evaluated during the docs
build.

## What's next

The [Reference](95-reference.md) page documents all exported types and
functions with their full parameter lists. The [Architecture](10-architecture.md)
document describes the package's design principles and the rationale for the
composition-over-hardcoding approach. The [Performance](20-performance.md) document
covers benchmarking methodology and identified optimization opportunities.

For further exploration, a few directions are worth trying once the basic
workflow is familiar:

- Change the dominance mode to `:dominant_first` or `:recessive_first` and
  compare the resulting cline shape to the codominant case. Dominance breaks
  the symmetry of the cline: the dominant allele tends to introgress further
  into the territory of the recessive one.
- Swap the lifecycle order to `:migration_first` by passing
  `lifecycle_order = :migration_first` to `simulate`, and compare the cline
  shape at high `s`. The ordering matters most when selection is strong
  relative to migration.
- Run longer simulations (`n_generations = 500` or more) to verify that the
  cline has genuinely reached equilibrium rather than still evolving slowly.
  At low `s` and high σ², equilibration can take many generations.

## Citation

Users who publish results generated with HybridZones.jl should cite both
the package itself and the foundational papers that the simulation framework
implements. The package CITATION.cff file lists the full details. The most
directly relevant methodological references are Mallet & Barton (1989), which
develops the frequency-dependent selection cline framework and the dispersal
variance estimator, and Mallet et al. (1990), which applies the framework to
empirical *Heliconius* transect data and derives the selection and gene flow
estimators.
