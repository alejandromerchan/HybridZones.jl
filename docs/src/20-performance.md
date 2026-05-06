# Performance

This page documents the performance characteristics of HybridZones.jl,
including identified optimization opportunities, profiling results,
and completed performance work. The goal is to track what's known
about the package's performance over time, both as a record of design
decisions and as a guide for future optimization work.

## Current performance

A representative single-locus simulation (101 demes, 200 generations,
σ²=30, s=0.3) on modern hardware:

- ~19 ms per simulation
- 6 allocations regardless of generation count (per-generation loop
  is allocation-free)
- Memory: ~5 KiB working buffers; ~480 KiB additional for trajectory
  storage when `save_trajectory=true`

This performance is adequate for typical research workflows — a
1000-combination parameter sweep completes in roughly 20 seconds.
The architecture document's aspirational target was ~10 ms per
simulation, so there's headroom for optimization, but no urgency.

## Identified optimization opportunities

The following optimizations have been identified through profiling
but deferred. They're documented here for future work.

### Migration inner loop: branch-free boundary handling

**Status:** Deferred. To be addressed after Pascal validation is in
place to ensure optimizations don't introduce subtle correctness
regressions.

**Profile evidence:** Profiling a 10,000-generation simulation showed
migration code consuming approximately 99% of total runtime. The
hottest single source location was the migration inner loop body
(MigrationModels.jl line 165), with concentrated samples in the
boundary-check branch and surrounding arithmetic.

**Current implementation:** The migration inner loop uses an
if-elseif-else cascade to handle deme indices outside `[1, n_demes]`,
substituting boundary values for "ghost" demes:

    if source < 1
        s += w * freq_current[1]
    elseif source > n_demes
        s += w * freq_current[n_demes]
    else
        s += w * freq_current[source]
    end

**Proposed optimization:** Replace the branch with a clamp operation
that's branch-free at the assembly level:

    clamped = clamp(source, 1, n_demes)
    s += w * freq_current[clamped]

This change should allow SIMD vectorization of the surrounding loop,
since there's no control flow inside. Combined with `@inbounds @simd`
annotations, expected speedup is 2–4x on migration code, which
translates to a similar speedup on the overall simulator since
migration is currently the dominant cost.

**Risks:** Numerical results should be identical (clamp produces the
same values as the conditional cascade), but care should be taken to
verify with the Pascal regression test suite once that's in place.
The `@inbounds` annotation eliminates bounds checking, which means
any indexing bug would silently produce wrong results rather than
throwing an error. Tests must cover boundary cases thoroughly before
this annotation is added.

### Per-row migration of matrix state

**Status:** Deferred. May be revisited if multi-locus simulations
show migration call overhead becoming significant.

**Background:** The current implementation calls `migrate!` once per
genotype row when applying migration to a matrix-valued state
(implemented as `for g in 1:n_genotypes; migrate!(@view(...), ...);
end`). For OneLocusDiploid (3 genotypes) this is negligible. For
the eventual three-locus types (27 genotypes for WC3ARBCD-style
simulations) the per-call overhead becomes a meaningful fraction
of total time.

**Proposed optimization:** Add a matrix-aware `migrate!` method
that operates on all rows simultaneously, potentially with better
cache reuse since the migration kernel is the same for all
genotypes.

## Completed performance work

(Empty for now. Optimizations move to this section after they're
implemented and validated.)

## Performance methodology

The package's performance is measured using BenchmarkTools.jl.
Standard benchmark setup:

    using BenchmarkTools
    using HybridZones

    arch = OneLocusDiploid(dominance = :codominant)
    sel = FrequencyDependentSelection(s = 0.3)
    mig = BinomialStepping(30.0)
    initial = secondary_contact(arch)

    @btime simulate($arch, $sel, $mig, $initial;
                    n_generations = 200,
                    save_trajectory = false)

Profiling uses Julia's built-in Profile module, typically with a
1ms sampling interval and a single longer simulation rather than
many short ones:

    using Profile
    Profile.init(delay = 0.001)
    Profile.clear()
    @profile simulate(arch, sel, mig, initial;
                      n_generations = 10000,
                      save_trajectory = false)
    Profile.print(format = :flat, sortedby = :count)

Performance regressions and improvements should be measured against
this standard configuration so results are comparable over time.
