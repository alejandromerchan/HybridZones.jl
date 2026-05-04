using HybridZones
using Test

# Shared fixtures used across multiple testsets
const _arch = OneLocusDiploid(dominance = :codominant)
const _sel = FrequencyDependentSelection(s = 0.3)
const _mig = BinomialStepping(30.0)
const _initial = secondary_contact(_arch)

# ──────────────────────────────────────────────────────────────
# secondary_contact
# ──────────────────────────────────────────────────────────────

@testset "secondary_contact — shape and column sums" begin
    state = secondary_contact(OneLocusDiploid())
    @test size(state) == (3, 101)
    for j in 1:101
        @test isapprox(sum(state[:, j]), 1.0; atol = 1e-12)
    end
end

@testset "secondary_contact — zone structure (default 101 demes)" begin
    state = secondary_contact(OneLocusDiploid())
    # left half: A1A1 only
    for j in 1:50
        @test state[1, j] == 1.0
        @test state[2, j] == 0.0
        @test state[3, j] == 0.0
    end
    # contact deme 51: all heterozygotes
    @test state[1, 51] == 0.0
    @test state[2, 51] == 1.0
    @test state[3, 51] == 0.0
    # right half: A2A2 only
    for j in 52:101
        @test state[1, j] == 0.0
        @test state[2, j] == 0.0
        @test state[3, j] == 1.0
    end
end

@testset "secondary_contact — custom n_demes" begin
    state = secondary_contact(OneLocusDiploid(); n_demes = 11)
    @test size(state) == (3, 11)
    # contact = (11+1)÷2 = 6; left = 1:5; right = 7:11
    for j in 1:5
        @test state[1, j] == 1.0
    end
    @test state[2, 6] == 1.0
    for j in 7:11
        @test state[3, j] == 1.0
    end
end

@testset "secondary_contact — invalid n_demes throws" begin
    @test_throws ArgumentError secondary_contact(OneLocusDiploid(); n_demes = 0)
    @test_throws ArgumentError secondary_contact(OneLocusDiploid(); n_demes = -5)
end

# ──────────────────────────────────────────────────────────────
# n_genotypes
# ──────────────────────────────────────────────────────────────

@testset "n_genotypes — OneLocusDiploid returns 3" begin
    @test n_genotypes(OneLocusDiploid()) == 3
    @test n_genotypes(OneLocusDiploid(dominance = :dominant_first)) == 3
end

# ──────────────────────────────────────────────────────────────
# allele_frequencies
# ──────────────────────────────────────────────────────────────

@testset "allele_frequencies — Hardy-Weinberg arithmetic" begin
    arch = OneLocusDiploid()
    # A2 freq = freq_A2A2 + 0.5 * freq_A1A2
    state = reshape([0.25, 0.50, 0.25], 3, 1)
    @test allele_frequencies(state, arch) ≈ [0.5]

    state2 = reshape([1.0, 0.0, 0.0], 3, 1)
    @test allele_frequencies(state2, arch) ≈ [0.0]

    state3 = reshape([0.0, 0.0, 1.0], 3, 1)
    @test allele_frequencies(state3, arch) ≈ [1.0]

    state4 = reshape([0.0, 1.0, 0.0], 3, 1)
    @test allele_frequencies(state4, arch) ≈ [0.5]
end

@testset "allele_frequencies — multi-deme state" begin
    arch = OneLocusDiploid()
    # 3 demes: A1A1-fixed, het, A2A2-fixed
    state = [1.0 0.0 0.0
             0.0 1.0 0.0
             0.0 0.0 1.0]
    freq = allele_frequencies(state, arch)
    @test freq ≈ [0.0, 0.5, 1.0]
    @test length(freq) == 3
end

@testset "allele_frequencies — 3D trajectory" begin
    arch = OneLocusDiploid()
    traj = simulate(_arch, _sel, _mig, _initial; n_generations = 5)
    freq_traj = allele_frequencies(traj, arch)
    @test size(freq_traj) == (101, 6)
    # first time slice matches allele_frequencies of initial state
    @test freq_traj[:, 1] ≈ allele_frequencies(_initial, arch)
end

# ──────────────────────────────────────────────────────────────
# simulate — argument validation
# ──────────────────────────────────────────────────────────────

@testset "simulate — invalid n_generations throws" begin
    @test_throws ArgumentError simulate(_arch, _sel, _mig, _initial; n_generations = 0)
    @test_throws ArgumentError simulate(_arch, _sel, _mig, _initial; n_generations = -1)
end

@testset "simulate — invalid lifecycle_order throws" begin
    @test_throws ArgumentError simulate(
        _arch, _sel, _mig, _initial; n_generations = 1, lifecycle_order = :random_order
    )
    @test_throws ArgumentError simulate(
        _arch, _sel, _mig, _initial; n_generations = 1, lifecycle_order = :both
    )
end

@testset "simulate — mismatched initial_state rows throws" begin
    wrong_rows = zeros(2, 101)  # 2 rows instead of 3
    @test_throws DimensionMismatch simulate(
        _arch, _sel, _mig, wrong_rows; n_generations = 1
    )
end

@testset "simulate — initial_state columns not summing to 1 throws" begin
    bad_state = copy(_initial)
    bad_state[1, 5] = 0.5  # column 5 now sums to 1.5
    @test_throws ArgumentError simulate(_arch, _sel, _mig, bad_state; n_generations = 1)
end

# ──────────────────────────────────────────────────────────────
# simulate — trajectory storage
# ──────────────────────────────────────────────────────────────

@testset "simulate — trajectory shape and initial slice" begin
    n_gens = 10
    traj = simulate(_arch, _sel, _mig, _initial; n_generations = n_gens)
    @test size(traj) == (3, 101, n_gens + 1)
    @test traj[:, :, 1] == _initial
end

@testset "simulate — save_trajectory=false returns 2D final state" begin
    n_gens = 10
    traj = simulate(_arch, _sel, _mig, _initial; n_generations = n_gens)
    final_traj = traj[:, :, end]
    final_only = simulate(
        _arch, _sel, _mig, _initial; n_generations = n_gens, save_trajectory = false
    )
    @test final_only isa Matrix
    @test size(final_only) == (3, 101)
    @test final_only ≈ final_traj
end

# ──────────────────────────────────────────────────────────────
# simulate — mass conservation
# ──────────────────────────────────────────────────────────────

@testset "simulate — mass conservation across generations" begin
    n_gens = 20
    traj = simulate(_arch, _sel, _mig, _initial; n_generations = n_gens)
    for g in 1:(n_gens + 1)
        for j in 1:101
            @test isapprox(sum(traj[:, j, g]), 1.0; atol = 1e-10)
        end
    end
end

@testset "simulate — mass conservation, migration_first" begin
    n_gens = 20
    traj = simulate(
        _arch, _sel, _mig, _initial;
        n_generations = n_gens, lifecycle_order = :migration_first
    )
    for g in 1:(n_gens + 1)
        for j in 1:101
            @test isapprox(sum(traj[:, j, g]), 1.0; atol = 1e-10)
        end
    end
end

# ──────────────────────────────────────────────────────────────
# simulate — lifecycle ordering produces different results
# ──────────────────────────────────────────────────────────────

@testset "simulate — lifecycle ordering produces different results" begin
    n_gens = 50
    traj_sf = simulate(
        _arch, _sel, _mig, _initial;
        n_generations = n_gens, lifecycle_order = :selection_first
    )
    traj_mf = simulate(
        _arch, _sel, _mig, _initial;
        n_generations = n_gens, lifecycle_order = :migration_first
    )
    # Final states must differ; check that maximum absolute difference is nonzero
    @test maximum(abs.(traj_sf[:, :, end] .- traj_mf[:, :, end])) > 1e-6
end

# ──────────────────────────────────────────────────────────────
# simulate — comparison with manual loop
# ──────────────────────────────────────────────────────────────

@testset "simulate — matches manual selection then migration loop" begin
    n_gens = 8
    arch = _arch
    sel = _sel
    mig = _mig
    initial = secondary_contact(arch)

    # Manual loop: select! then migrate! row-by-row
    manual = copy(initial)
    buf = similar(manual)
    for _ in 1:n_gens
        select!(buf, manual, sel, arch)
        for i in 1:3
            migrate!(@view(manual[i, :]), @view(buf[i, :]), mig)
        end
    end

    simulated = simulate(arch, sel, mig, initial; n_generations = n_gens,
        save_trajectory = false)
    @test simulated ≈ manual
end

@testset "simulate — matches manual migration then selection loop" begin
    n_gens = 8
    arch = _arch
    sel = _sel
    mig = _mig
    initial = secondary_contact(arch)

    # Manual loop: migrate! row-by-row then select!
    manual = copy(initial)
    buf = similar(manual)
    for _ in 1:n_gens
        for i in 1:3
            migrate!(@view(buf[i, :]), @view(manual[i, :]), mig)
        end
        select!(manual, buf, sel, arch)
    end

    simulated = simulate(arch, sel, mig, initial; n_generations = n_gens,
        lifecycle_order = :migration_first, save_trajectory = false)
    @test simulated ≈ manual
end

# ──────────────────────────────────────────────────────────────
# simulate — equilibrium convergence
# ──────────────────────────────────────────────────────────────

@testset "simulate — equilibrium convergence (s=0.3, σ²=30, 200 gen)" begin
    arch = OneLocusDiploid(dominance = :codominant)
    sel = FrequencyDependentSelection(s = 0.3)
    mig = BinomialStepping(30.0)
    initial = secondary_contact(arch)

    state200 = simulate(arch, sel, mig, initial; n_generations = 200,
        save_trajectory = false)
    freq200 = allele_frequencies(state200, arch)

    # Sigmoid shape: A2 frequency increases monotonically across transect
    for j in 1:100
        @test freq200[j] ≤ freq200[j + 1] + 1e-10
    end

    # Centred at deme 51: A2 frequency ≈ 0.5 at the contact deme
    @test isapprox(freq200[51], 0.5; atol = 0.05)

    # Symmetric around centre: freq[51-k] ≈ 1 - freq[51+k]
    for k in 1:10
        @test isapprox(freq200[51 - k], 1.0 - freq200[51 + k]; atol = 0.02)
    end
end

@testset "simulate — near-equilibrium after 200 generations" begin
    arch = OneLocusDiploid(dominance = :codominant)
    sel = FrequencyDependentSelection(s = 0.3)
    mig = BinomialStepping(30.0)
    initial = secondary_contact(arch)

    state200 = simulate(arch, sel, mig, initial; n_generations = 200,
        save_trajectory = false)
    state400 = simulate(arch, sel, mig, initial; n_generations = 400,
        save_trajectory = false)

    freq200 = allele_frequencies(state200, arch)
    freq400 = allele_frequencies(state400, arch)

    # Running 200 more generations changes frequencies by less than 0.1%
    @test maximum(abs.(freq400 .- freq200)) < 1e-3
end

# ──────────────────────────────────────────────────────────────
# simulate — initial_state not mutated
# ──────────────────────────────────────────────────────────────

@testset "simulate — does not mutate initial_state" begin
    initial = secondary_contact(_arch)
    initial_copy = copy(initial)
    simulate(_arch, _sel, _mig, initial; n_generations = 20, save_trajectory = false)
    @test initial == initial_copy
end
