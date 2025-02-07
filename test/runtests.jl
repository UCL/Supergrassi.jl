using Supergrassi
using Test

@testset "Supergrassi.jl" begin
    # Write your tests here.
end

@testset "Compute demand" begin

    n = 16
    elasticity = 2.0
    log_price = (uk = ones(n), eu = ones(n), world = ones(n))
    goods_consumption = (uk = ones(n), eu = ones(n), world = ones(n))

    alpha = Supergrassi.compute_demand.(elasticity,
                                        log_price.uk,
                                        log_price.eu,
                                        log_price.world,
                                        goods_consumption.uk,
                                        goods_consumption.eu,
                                        goods_consumption.world)

    @test length(alpha) == n
    for industry in alpha
        @test isfinite(sum(industry))
    end

    backend = AutoZygote()
    pd = [ones(n), 2*ones(n), 3*ones(n)]

    f(el,peu,pw,cd,ceu,cworld) = pd -> Supergrassi.compute_demand(el,pd,peu,pw,cd,ceu,cworld)

    g = gradient.(f.(elasticity,
                     log_price.eu,
                     log_price.world,
                     goods_consumption.uk,
                     goods_consumption.eu,
                     goods_consumption.world),
                  backend,
                  pd)
    
end

@testset "Compute equilibrium test" begin

    using Distributions
    using Random
    
    n = 16
    log_consumer_price_index = zeros(n)
    rand!(Normal(), log_consumer_price_index)
    
end
