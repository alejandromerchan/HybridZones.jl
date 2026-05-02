using HybridZones
using Test

# Access the internal kernel accessor via the submodule
kernel = HybridZones.MigrationModels.kernel

@testset "type hierarchy" begin
    m = BinomialStepping(10.0)
    @test m isa MigrationModel
    @test m isa BinomialStepping
end

@testset "construction — kernel sums to 1" begin
    for σ² in (10.0, 30.0, 100.0)
        m = BinomialStepping(σ²)
        @test isapprox(sum(kernel(m)), 1.0; atol = 1e-12)
    end
end

@testset "construction — kernel length is 2n+1" begin
    m = BinomialStepping(10.0)
    @test length(kernel(m)) == 2 * max_distance(m) + 1
end

@testset "construction — n derived from σ²" begin
    # α=0.5: n = ceil(σ² / (2*0.5*0.5)) = ceil(σ² / 0.5) = ceil(2σ²)
    m = BinomialStepping(10.0)
    @test max_distance(m) == ceil(Int, 10.0 / (2 * 0.5 * 0.5))  # = 20
    m2 = BinomialStepping(30.0)
    @test max_distance(m2) == ceil(Int, 30.0 / (2 * 0.5 * 0.5))  # = 60
end

@testset "accessor functions" begin
    m = BinomialStepping(10.0)
    @test migration_variance(m) === 10.0
    @test max_distance(m) == 20
    k = kernel(m)
    @test k isa Vector{Float64}
    @test length(k) == 41  # 2*20+1
end

@testset "validation — σ² ≤ 0 throws" begin
    @test_throws ArgumentError BinomialStepping(0.0)
    @test_throws ArgumentError BinomialStepping(-1.0)
end

@testset "validation — α outside (0,1) throws" begin
    @test_throws ArgumentError BinomialStepping(10.0; α = 0.0)
    @test_throws ArgumentError BinomialStepping(10.0; α = 1.0)
    @test_throws ArgumentError BinomialStepping(10.0; α = -0.1)
    @test_throws ArgumentError BinomialStepping(10.0; α = 1.1)
end

@testset "migrate! — uniform frequency is invariant" begin
    m = BinomialStepping(10.0)
    for p in (0.0, 0.3, 0.5, 1.0)
        freq = fill(p, 100)
        freq_next = similar(freq)
        migrate!(freq_next, freq, m)
        @test all(isapprox.(freq_next, p; atol = 1e-10))
    end
end

@testset "migrate! — step function: mass conservation" begin
    m = BinomialStepping(10.0)
    n_demes = 100
    freq = vcat(zeros(n_demes ÷ 2), ones(n_demes ÷ 2))
    freq_next = similar(freq)
    migrate!(freq_next, freq, m)
    @test isapprox(sum(freq_next), sum(freq); atol = 1e-10)
end

@testset "migrate! — step function: values in [0, 1]" begin
    m = BinomialStepping(10.0)
    n_demes = 100
    freq = vcat(zeros(n_demes ÷ 2), ones(n_demes ÷ 2))
    freq_next = similar(freq)
    migrate!(freq_next, freq, m)
    @test all(x -> -1e-12 ≤ x ≤ 1 + 1e-12, freq_next)
end

@testset "migrate! — step function: midpoint near 0.5" begin
    m = BinomialStepping(10.0)
    n_demes = 100
    freq = vcat(zeros(n_demes ÷ 2), ones(n_demes ÷ 2))
    freq_next = similar(freq)
    migrate!(freq_next, freq, m)
    # Midpoint demes 50 and 51 straddle the step; after one migration step
    # they should be pulled toward 0.5 from their original 0 and 1.
    @test freq_next[50] > 0.0
    @test freq_next[51] < 1.0
end

@testset "migrate! — larger σ² gives more smoothing" begin
    n_demes = 100
    freq = vcat(zeros(n_demes ÷ 2), ones(n_demes ÷ 2))

    m_small = BinomialStepping(5.0)
    m_large = BinomialStepping(20.0)

    freq_next_small = similar(freq)
    freq_next_large = similar(freq)
    migrate!(freq_next_small, freq, m_small)
    migrate!(freq_next_large, freq, m_large)

    # Larger σ² moves the step-50 deme frequency further from 0 toward 0.5.
    @test freq_next_large[50] > freq_next_small[50]
    # And deme 51 further from 1 toward 0.5.
    @test freq_next_large[51] < freq_next_small[51]
end

@testset "migrate! — boundary reflection: left edge near 1.0" begin
    m = BinomialStepping(10.0)
    n_demes = 100
    # Left half = 1.0, right half = 0.0
    freq = vcat(ones(n_demes ÷ 2), zeros(n_demes ÷ 2))
    freq_next = similar(freq)
    migrate!(freq_next, freq, m)
    # Leftmost deme is surrounded by 1.0 (via boundary reflection) so should
    # remain very close to 1.0.
    @test isapprox(freq_next[1], 1.0; atol = 1e-10)
end

@testset "migrate! — boundary reflection: right edge near 0.0" begin
    m = BinomialStepping(10.0)
    n_demes = 100
    freq = vcat(ones(n_demes ÷ 2), zeros(n_demes ÷ 2))
    freq_next = similar(freq)
    migrate!(freq_next, freq, m)
    # Rightmost deme is surrounded by 0.0 via boundary reflection.
    @test isapprox(freq_next[n_demes], 0.0; atol = 1e-10)
end

@testset "migrate! — DimensionMismatch on length mismatch" begin
    m = BinomialStepping(10.0)
    @test_throws DimensionMismatch migrate!(zeros(10), zeros(20), m)
end

@testset "show — informative representation" begin
    m = BinomialStepping(10.0)
    s = sprint(show, m)
    @test contains(s, "BinomialStepping")
    @test contains(s, "10.0")
end
