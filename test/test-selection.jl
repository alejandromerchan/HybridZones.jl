using HybridZones
using Test

@testset "type hierarchy" begin
    model = FrequencyDependentSelection(s = 0.1)
    @test model isa SelectionModel
    @test model isa FrequencyDependentSelection
end

@testset "construction — valid s" begin
    @test FrequencyDependentSelection(s = 0.1) isa FrequencyDependentSelection
    @test FrequencyDependentSelection(s = 0.6) isa FrequencyDependentSelection
    @test FrequencyDependentSelection(s = 0.0) isa FrequencyDependentSelection
    @test FrequencyDependentSelection(s = 1.0) isa FrequencyDependentSelection
end

@testset "construction — invalid s throws" begin
    @test_throws ArgumentError FrequencyDependentSelection(s = -0.1)
    @test_throws ArgumentError FrequencyDependentSelection(s = 1.5)
end

@testset "accessor — selection_coefficient" begin
    model = FrequencyDependentSelection(s = 0.3)
    @test selection_coefficient(model) === 0.3
end

@testset "show — informative representation" begin
    model = FrequencyDependentSelection(s = 0.5)
    s = sprint(show, model)
    @test contains(s, "FrequencyDependentSelection")
    @test contains(s, "0.5")
end

@testset "select! — mass conservation, all dominance modes" begin
    model = FrequencyDependentSelection(s = 0.5)
    # 3×5 matrix; columns represent demes with varied genotype distributions
    freq = [0.6 0.2 0.1 0.0 1.0
            0.3 0.5 0.7 0.5 0.0
            0.1 0.3 0.2 0.5 0.0]
    freq_next = similar(freq)
    for dom in (:codominant, :dominant_first, :recessive_first)
        arch = OneLocusDiploid(dominance = dom)
        select!(freq_next, freq, model, arch)
        for j in 1:size(freq, 2)
            @test isapprox(sum(freq_next[:, j]), 1.0; atol = 1e-12)
        end
    end
end

@testset "select! — no selection (s=0)" begin
    model = FrequencyDependentSelection(s = 0.0)
    freq = reshape([0.2, 0.5, 0.3], 3, 1)
    freq_next = similar(freq)
    for dom in (:codominant, :dominant_first, :recessive_first)
        arch = OneLocusDiploid(dominance = dom)
        select!(freq_next, freq, model, arch)
        @test freq_next ≈ freq
    end
end

@testset "select! — stable polymorphism at uniform frequencies (codominant)" begin
    # At (1/3, 1/3, 1/3), all phenotypes are equally common, so all fitnesses
    # are equal and selection has no directional effect. This is a fixed point.
    model = FrequencyDependentSelection(s = 0.5)
    arch = OneLocusDiploid(dominance = :codominant)
    freq = reshape([1 / 3, 1 / 3, 1 / 3], 3, 1)
    freq_next = similar(freq)
    for _ in 1:1000
        select!(freq_next, freq, model, arch)
        freq, freq_next = freq_next, freq
    end
    @test all(isapprox.(freq, 1 / 3; atol = 1e-10))
end

@testset "select! — rare morph penalty (codominant)" begin
    # A2A2 is rare: its phenotype is rare, so it has low fitness and decreases.
    # A1A1 is dominant: it has high fitness and increases.
    model = FrequencyDependentSelection(s = 0.5)
    arch = OneLocusDiploid(dominance = :codominant)
    freq = reshape([0.99, 0.0, 0.01], 3, 1)
    freq_next = similar(freq)
    select!(freq_next, freq, model, arch)
    @test freq_next[3, 1] < freq[3, 1]
    @test freq_next[1, 1] > freq[1, 1]
end

@testset "select! — s=0 makes all dominance modes produce equal output" begin
    model = FrequencyDependentSelection(s = 0.0)
    freq = reshape([0.2, 0.5, 0.3], 3, 1)
    freq_next_cod = similar(freq)
    freq_next_dom = similar(freq)
    freq_next_rec = similar(freq)
    select!(freq_next_cod, freq, model, OneLocusDiploid(dominance = :codominant))
    select!(freq_next_dom, freq, model, OneLocusDiploid(dominance = :dominant_first))
    select!(freq_next_rec, freq, model, OneLocusDiploid(dominance = :recessive_first))
    @test freq_next_cod ≈ freq
    @test freq_next_dom ≈ freq
    @test freq_next_rec ≈ freq
end

@testset "select! — dominant_first: A1A1 and A1A2 have equal growth rates" begin
    # A1A1 and A1A2 share the A1 phenotype, so they receive the same fitness
    # and therefore have the same growth rate p'_ij / p_ij = w_A / w_mean.
    model = FrequencyDependentSelection(s = 0.5)
    arch = OneLocusDiploid(dominance = :dominant_first)
    freq = reshape([0.4, 0.3, 0.3], 3, 1)
    freq_next = similar(freq)
    select!(freq_next, freq, model, arch)
    growth_11 = freq_next[1, 1] / freq[1, 1]
    growth_12 = freq_next[2, 1] / freq[2, 1]
    @test isapprox(growth_11, growth_12; atol = 1e-12)
end

@testset "select! — recessive_first: A1A2 and A2A2 have equal growth rates" begin
    # A1A2 and A2A2 share the A2 phenotype under recessive_first dominance.
    model = FrequencyDependentSelection(s = 0.5)
    arch = OneLocusDiploid(dominance = :recessive_first)
    freq = reshape([0.3, 0.4, 0.3], 3, 1)
    freq_next = similar(freq)
    select!(freq_next, freq, model, arch)
    growth_12 = freq_next[2, 1] / freq[2, 1]
    growth_22 = freq_next[3, 1] / freq[3, 1]
    @test isapprox(growth_12, growth_22; atol = 1e-12)
end

@testset "select! — DimensionMismatch on size mismatch" begin
    model = FrequencyDependentSelection(s = 0.3)
    arch = OneLocusDiploid()
    @test_throws DimensionMismatch select!(zeros(3, 5), zeros(3, 10), model, arch)
end

@testset "select! — ArgumentError on wrong number of rows" begin
    model = FrequencyDependentSelection(s = 0.3)
    arch = OneLocusDiploid()
    @test_throws ArgumentError select!(zeros(2, 5), zeros(2, 5), model, arch)
end

# ── SemiDominantFrequencyDependentSelection ─────────────────────────────────

@testset "SemiDominantFDS — type hierarchy" begin
    model = SemiDominantFrequencyDependentSelection(0.3)
    @test model isa SelectionModel
    @test model isa SemiDominantFrequencyDependentSelection
end

@testset "SemiDominantFDS — construction valid s" begin
    @test SemiDominantFrequencyDependentSelection(0.0) isa
          SemiDominantFrequencyDependentSelection
    @test SemiDominantFrequencyDependentSelection(0.3) isa
          SemiDominantFrequencyDependentSelection
    @test SemiDominantFrequencyDependentSelection(1.0) isa
          SemiDominantFrequencyDependentSelection
end

@testset "SemiDominantFDS — construction invalid s throws" begin
    @test_throws ArgumentError SemiDominantFrequencyDependentSelection(-0.1)
    @test_throws ArgumentError SemiDominantFrequencyDependentSelection(1.5)
end

@testset "SemiDominantFDS — selection_coefficient accessor" begin
    model = SemiDominantFrequencyDependentSelection(0.3)
    @test selection_coefficient(model) === 0.3
end

@testset "SemiDominantFDS — basic correctness (hand-computed)" begin
    # p11=0.6, p12=0.3, p22=0.1, s=0.3
    # sim_11 = 0.6 + 0.5*0.3 = 0.75,  w11 = 1 - 0.3*(1-0.75) = 0.925
    # sim_12 = 0.3 + 0.5*0.7 = 0.65,  w12 = 1 - 0.3*(1-0.65) = 0.895
    # sim_22 = 0.1 + 0.5*0.3 = 0.25,  w22 = 1 - 0.3*(1-0.25) = 0.775
    # w_mean = 0.925*0.6 + 0.895*0.3 + 0.775*0.1 = 0.901
    model = SemiDominantFrequencyDependentSelection(0.3)
    arch = OneLocusDiploid(dominance = :codominant)
    freq = reshape([0.6, 0.3, 0.1], 3, 1)
    freq_next = similar(freq)
    select!(freq_next, freq, model, arch)
    @test isapprox(freq_next[1, 1], 0.925 * 0.6 / 0.901; atol = 1e-12)
    @test isapprox(freq_next[2, 1], 0.895 * 0.3 / 0.901; atol = 1e-12)
    @test isapprox(freq_next[3, 1], 0.775 * 0.1 / 0.901; atol = 1e-12)
end

@testset "SemiDominantFDS — mass conservation" begin
    model = SemiDominantFrequencyDependentSelection(0.5)
    arch = OneLocusDiploid(dominance = :codominant)
    freq = [0.6 0.2 0.1 0.0 1.0
            0.3 0.5 0.7 0.5 0.0
            0.1 0.3 0.2 0.5 0.0]
    freq_next = similar(freq)
    select!(freq_next, freq, model, arch)
    for j in 1:size(freq, 2)
        @test isapprox(sum(freq_next[:, j]), 1.0; atol = 1e-12)
    end
end

@testset "SemiDominantFDS — symmetry preserved" begin
    # Symmetric input (p11 = p22) stays symmetric after selection
    model = SemiDominantFrequencyDependentSelection(0.3)
    arch = OneLocusDiploid(dominance = :codominant)
    freq = reshape([0.25, 0.50, 0.25], 3, 1)
    freq_next = similar(freq)
    select!(freq_next, freq, model, arch)
    @test isapprox(freq_next[1, 1], freq_next[3, 1]; atol = 1e-14)
end

@testset "SemiDominantFDS — no selection (s=0)" begin
    model = SemiDominantFrequencyDependentSelection(0.0)
    arch = OneLocusDiploid(dominance = :codominant)
    freq = reshape([0.2, 0.5, 0.3], 3, 1)
    freq_next = similar(freq)
    select!(freq_next, freq, model, arch)
    @test freq_next ≈ freq
end

@testset "SemiDominantFDS — dispatch error on non-codominant dominance" begin
    model = SemiDominantFrequencyDependentSelection(0.3)
    freq = reshape([0.4, 0.4, 0.2], 3, 1)
    freq_next = similar(freq)
    @test_throws ArgumentError select!(
        freq_next, freq, model, OneLocusDiploid(dominance = :dominant_first))
    @test_throws ArgumentError select!(
        freq_next, freq, model, OneLocusDiploid(dominance = :recessive_first))
end

@testset "SemiDominantFDS — differs from FrequencyDependentSelection" begin
    # At non-trivial input the two models produce different post-selection states
    fds = FrequencyDependentSelection(s = 0.3)
    sfds = SemiDominantFrequencyDependentSelection(0.3)
    arch = OneLocusDiploid(dominance = :codominant)
    freq = reshape([0.6, 0.3, 0.1], 3, 1)
    out_fds = similar(freq)
    out_sfds = similar(freq)
    select!(out_fds, freq, fds, arch)
    select!(out_sfds, freq, sfds, arch)
    @test !(out_fds ≈ out_sfds)
end

@testset "SemiDominantFDS — reproducibility (identical input → identical output)" begin
    model = SemiDominantFrequencyDependentSelection(0.3)
    arch = OneLocusDiploid(dominance = :codominant)
    freq = reshape([0.6, 0.3, 0.1], 3, 1)
    out1 = similar(freq)
    out2 = similar(freq)
    select!(out1, freq, model, arch)
    select!(out2, freq, model, arch)
    @test out1 == out2
end

@testset "SemiDominantFDS — DimensionMismatch on size mismatch" begin
    model = SemiDominantFrequencyDependentSelection(0.3)
    arch = OneLocusDiploid(dominance = :codominant)
    @test_throws DimensionMismatch select!(zeros(3, 5), zeros(3, 10), model, arch)
end

@testset "SemiDominantFDS — ArgumentError on wrong number of rows" begin
    model = SemiDominantFrequencyDependentSelection(0.3)
    arch = OneLocusDiploid(dominance = :codominant)
    @test_throws ArgumentError select!(zeros(2, 5), zeros(2, 5), model, arch)
end

@testset "Pascal cline width reproducibility" begin
    arch = OneLocusDiploid(dominance = :codominant)
    sel = SemiDominantFrequencyDependentSelection(0.3)
    mig = BinomialStepping(30.0)
    mat = RandomMating()
    initial = secondary_contact(arch)

    trajectory = simulate(arch, mat, sel, mig, initial;
        n_generations = 200,
        lifecycle_order = :mate_migrate_select)

    p_final = allele_frequencies(trajectory[:, :, end], arch)

    function find_crossing(p, target)
        for i in 1:(length(p) - 1)
            if (p[i] - target) * (p[i + 1] - target) < 0
                return i + (p[i] - target) / (p[i] - p[i + 1])
            end
        end
        return nothing
    end

    d10 = find_crossing(p_final, 0.1)
    d90 = find_crossing(p_final, 0.9)
    width = abs(d90 - d10)

    # Pascal width is approximately 46.02 demes; HybridZones with
    # semi-dominant selection should reproduce this within ~0.2 demes
    @test isapprox(width, 46.02; atol = 0.2)
end
