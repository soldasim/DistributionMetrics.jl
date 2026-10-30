using DistributionMetrics
using Test
using Aqua

@testset "DistributionMetrics.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(DistributionMetrics)
    end
    # Write your tests here.
end
