using Supergrassi
using Test
using Enzyme

@testset "Supergrassi.jl" begin
    # Write your tests here.
end

@testset "Compute demand" begin

    n = 16
    elasticity = 2.0
    log_price = (uk = randn(n), eu = randn(n), world = randn(n))
    goods_consumption = (uk = ones(n), eu = ones(n), world = ones(n))

    alpha = Supergrassi.compute_demand.(elasticity,
                                        log_price.uk,
                                        log_price.eu,
                                        log_price.world,
                                        goods_consumption.uk,
                                        goods_consumption.eu,
                                        goods_consumption.world)

    @test length(alpha) == n
    @test all(i -> minimum(i) > 0, alpha)
    
    f(el, pd, peu, pw, fd, feu, fw) = gradient(Forward,
                                               Supergrassi.compute_demand,
                                               Const(el),
                                               pd,
                                               Const(peu),
                                               Const(pw),
                                               Const(fd),
                                               Const(feu),
                                               Const(fw))
    
    grad_alpha = f.(elasticity,
                    log_price.uk,
                    log_price.eu,
                    log_price.world,
                    goods_consumption.uk,
                    goods_consumption.eu,
                    goods_consumption.world)

    # Test grad_alpha is of right shape
    @test length(grad_alpha) == n
    @test all(i -> length(i) == 7, grad_alpha)
    # Test gradients of Consts are `nothing`
    @test all(i -> isnothing(i[1]), grad_alpha)
    @test all(i -> all(j -> isnothing(j), i[3:7]), grad_alpha)
    # Test gradients of `pd` are finite
    @test all(i -> all(j -> isfinite(j), i[2]), grad_alpha)
end

# @testset "Compute equilibrium test" begin

#     using Distributions
#     using Random
    
#     n = 16
#     log_consumer_price_index = zeros(n)
#     rand!(Normal(), log_consumer_price_index)
    
# end
