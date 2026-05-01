using HybridZones
using Test

@testset "type hierarchy" begin
    arch = OneLocusDiploid()
    @test arch isa GeneticArchitecture
    @test arch isa OneLocusDiploid
end

@testset "default construction" begin
    arch = OneLocusDiploid()
    @test dominance(arch) == :codominant
    @test allele_names(arch) == ("A1", "A2")
    @test n_loci(arch) == 1
    @test n_alleles(arch) == 2
end

@testset "valid dominance values" begin
    for d in (:codominant, :dominant_first, :recessive_first)
        arch = OneLocusDiploid(dominance = d)
        @test dominance(arch) == d
    end
end

@testset "custom allele names" begin
    arch = OneLocusDiploid(allele_names = ("B1", "B2"))
    @test allele_names(arch) == ("B1", "B2")
end

@testset "invalid dominance throws" begin
    @test_throws ArgumentError OneLocusDiploid(dominance = :dominant)
    @test_throws ArgumentError OneLocusDiploid(dominance = :recessive)
    @test_throws ArgumentError OneLocusDiploid(dominance = :additive)
end
