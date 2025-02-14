using Supergrassi
using Test
using Enzyme
using CSV, DataFrames

@testset "Supergrassi.jl" begin
    # Write your tests here.
end

@testset "Compute demand" begin

    n = 16
    elasticity = 2.0

    alpha = Matrix{Float64}(undef, n, 6)
    grad_alpha = Matrix{Float64}(undef, n, 6)
    df = CSV.read("../data/data_for_compute_demand.csv", DataFrame)

    f(el, pd, peu, pw, fd, feu, fw) = gradient(Forward,
                                               Supergrassi.compute_demand,
                                               Const(el),
                                               pd,
                                               Const(peu),
                                               Const(pw),
                                               Const(fd),
                                               Const(feu),
                                               Const(fw))

    for row in axes(alpha, 1)
        alpha[row,:] .= Supergrassi.compute_demand(elasticity,
                                                   df.logP_uk[row],
                                                   df.logP_eu[row],
                                                   df.logP_w[row],
                                                   df.f_uk[row],
                                                   df.f_eu[row],
                                                   df.f_w[row])

        grad_alpha[row,:] .= f(elasticity,
                               df.logP_uk[row],
                               df.logP_eu[row],
                               df.logP_w[row],
                               df.f_uk[row],
                               df.f_eu[row],
                               df.f_w[row])[2]
    end

    alpha[12,2:6] .= 0.0
    grad_alpha[12,:] .= 0.0
    
    @test isapprox(alpha[:,1:3], [df.alpha_uk df.alpha_eu df.alpha_w], atol = 1e-4)
    @test isapprox(grad_alpha[:,4:6], [df.dlogalpha_uk df.dlogalpha_eu df.dlogalpha_w], atol = 1e-4)
end

# @testset "Compute equilibrium test" begin

#     using Distributions
#     using Random
    
#     n = 16
#     log_consumer_price_index = zeros(n)
#     rand!(Normal(), log_consumer_price_index)
    
# end
