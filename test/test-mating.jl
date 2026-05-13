using HybridZones
using Test

# ──────────────────────────────────────────────────────────────
# Construction and type hierarchy
# ──────────────────────────────────────────────────────────────

@testset "RandomMating — construction and type hierarchy" begin
    mat = RandomMating()
    @test mat isa MatingModel
    @test mat isa RandomMating
    @test RandomMating <: MatingModel
end

# ──────────────────────────────────────────────────────────────
# HW restoration
# ──────────────────────────────────────────────────────────────

@testset "mate! — Hardy-Weinberg restoration from non-HW state" begin
    arch = OneLocusDiploid()
    mat = RandomMating()
    # Non-HW state: excess heterozygotes (Wahlund-effect deficit)
    # p = 0.6 + 0.5*0.3 = 0.75; HW expects (0.5625, 0.375, 0.0625)
    g = reshape([0.6, 0.3, 0.1], 3, 1)
    g_next = similar(g)
    mate!(g_next, g, mat, arch)
    p = 0.6 + 0.5 * 0.3  # = 0.75
    q = 1.0 - p
    @test g_next[1, 1] ≈ p^2
    @test g_next[2, 1] ≈ 2 * p * q
    @test g_next[3, 1] ≈ q^2
end

# ──────────────────────────────────────────────────────────────
# Mass conservation
# ──────────────────────────────────────────────────────────────

@testset "mate! — mass conservation (column sums = 1)" begin
    arch = OneLocusDiploid()
    mat = RandomMating()
    n_demes = 10
    g = zeros(3, n_demes)
    for j in 1:n_demes
        p = j / (n_demes + 1)
        q = 1.0 - p
        g[1, j] = p^2
        g[2, j] = 2 * p * q
        g[3, j] = q^2
    end
    g_next = similar(g)
    mate!(g_next, g, mat, arch)
    for j in 1:n_demes
        @test isapprox(sum(g_next[:, j]), 1.0; atol = 1e-14)
    end
end

# ──────────────────────────────────────────────────────────────
# Allele frequency preservation
# ──────────────────────────────────────────────────────────────

@testset "mate! — allele frequency preservation per deme" begin
    arch = OneLocusDiploid()
    mat = RandomMating()
    # Non-HW state with known allele frequencies
    # Deme 1: p_A1 = 0.7; deme 2: p_A1 = 0.2
    g = [0.6 0.1
         0.2 0.2
         0.2 0.7]
    g_next = similar(g)
    mate!(g_next, g, mat, arch)
    # A1 frequency = g[1] + 0.5*g[2]
    p1_before = g[1, 1] + 0.5 * g[2, 1]
    p1_after = g_next[1, 1] + 0.5 * g_next[2, 1]
    @test p1_before ≈ p1_after atol = 1e-14

    p2_before = g[1, 2] + 0.5 * g[2, 2]
    p2_after = g_next[1, 2] + 0.5 * g_next[2, 2]
    @test p2_before ≈ p2_after atol = 1e-14
end

# ──────────────────────────────────────────────────────────────
# Idempotence
# ──────────────────────────────────────────────────────────────

@testset "mate! — idempotence: applying twice equals applying once" begin
    arch = OneLocusDiploid()
    mat = RandomMating()
    g = reshape([0.6, 0.3, 0.1], 3, 1)
    g1 = similar(g)
    g2 = similar(g)
    mate!(g1, g, mat, arch)
    mate!(g2, g1, mat, arch)
    @test g2 ≈ g1
end

# ──────────────────────────────────────────────────────────────
# Edge cases
# ──────────────────────────────────────────────────────────────

@testset "mate! — edge case: all A1A1 (p=1)" begin
    arch = OneLocusDiploid()
    mat = RandomMating()
    g = reshape([1.0, 0.0, 0.0], 3, 1)
    g_next = similar(g)
    mate!(g_next, g, mat, arch)
    @test g_next[1, 1] ≈ 1.0
    @test g_next[2, 1] ≈ 0.0
    @test g_next[3, 1] ≈ 0.0
end

@testset "mate! — edge case: all A2A2 (p=0)" begin
    arch = OneLocusDiploid()
    mat = RandomMating()
    g = reshape([0.0, 0.0, 1.0], 3, 1)
    g_next = similar(g)
    mate!(g_next, g, mat, arch)
    @test g_next[1, 1] ≈ 0.0
    @test g_next[2, 1] ≈ 0.0
    @test g_next[3, 1] ≈ 1.0
end

@testset "mate! — edge case: all heterozygotes (p=0.5)" begin
    arch = OneLocusDiploid()
    mat = RandomMating()
    g = reshape([0.0, 1.0, 0.0], 3, 1)
    g_next = similar(g)
    mate!(g_next, g, mat, arch)
    @test g_next[1, 1] ≈ 0.25
    @test g_next[2, 1] ≈ 0.50
    @test g_next[3, 1] ≈ 0.25
end

# ──────────────────────────────────────────────────────────────
# Dominance-independence
# ──────────────────────────────────────────────────────────────

@testset "mate! — dominance-independent: all modes give identical output" begin
    mat = RandomMating()
    g = reshape([0.6, 0.3, 0.1], 3, 1)

    arch_co = OneLocusDiploid(dominance = :codominant)
    arch_d1 = OneLocusDiploid(dominance = :dominant_first)
    arch_r1 = OneLocusDiploid(dominance = :recessive_first)

    g_co = similar(g)
    g_d1 = similar(g)
    g_r1 = similar(g)

    mate!(g_co, g, mat, arch_co)
    mate!(g_d1, g, mat, arch_d1)
    mate!(g_r1, g, mat, arch_r1)

    @test g_co ≈ g_d1
    @test g_co ≈ g_r1
end

# ──────────────────────────────────────────────────────────────
# Input validation
# ──────────────────────────────────────────────────────────────

@testset "mate! — mismatched dimensions throw DimensionMismatch" begin
    arch = OneLocusDiploid()
    mat = RandomMating()
    g = zeros(3, 5)
    g_wrong = zeros(3, 7)
    @test_throws DimensionMismatch mate!(g_wrong, g, mat, arch)
end

@testset "mate! — wrong number of genotype rows throws ArgumentError" begin
    arch = OneLocusDiploid()
    mat = RandomMating()
    g = zeros(2, 5)
    g_next = zeros(2, 5)
    @test_throws ArgumentError mate!(g_next, g, mat, arch)
end
