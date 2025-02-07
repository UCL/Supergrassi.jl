using Supergrassi
using Test

@testset "Supergrassi.jl" begin
    # Write your tests here.
end

@testset "Compute demand" begin

    n = 16
    elasticity = (_ = 3.0, a = 2.0)
    log_price = (uk = ones(n), eu = ones(n), world = ones(n))
    goods_consumption = (composite = ones(n), uk = ones(n), eu = ones(n), world = ones(n))

    alpha = Supergrassi.compute_demand(n,
                                       elasticity,
                                       log_price.uk,
                                       log_price.eu,
                                       log_price.world,
                                       goods_consumption)

    @test haskey(alpha, "composite")
    @test haskey(alpha, "uk")
    @test haskey(alpha, "eu")
    @test haskey(alpha, "world")

    for key in keys(alpha)
        @test length(alpha[key]) == n
        @test isfinite(sum(alpha[key]))
    end
    
end

@testset "Compute equilibrium test" begin

    using Distributions
    using Random
    
    n = 16
    log_consumer_price_index = zeros(n)
    rand!(Normal(), log_consumer_price_index)
    
end
