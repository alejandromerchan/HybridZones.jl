using HybridZones
using Test

@testset "package_version" begin
    v = package_version()
    @test v isa VersionNumber
    @test v == v"0.1.0"
end
