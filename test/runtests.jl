using Supergrassi
using Test
using Enzyme
using CSV, DataFrames

@testset "Supergrassi.jl" begin
    # Write your tests here.
end

@testset "Compute consumption" begin

    n = 16
    elasticity_a = 2.0
    elasticity = 3.0
    tol = 1e-12
    
    alpha = Matrix{Float64}(undef, n, 6)
    grad_alpha = Matrix{Float64}(undef, n, 6)
    df = CSV.read("../data/data_for_compute_demand.csv", DataFrame)

    gradient_function(el, puk, peu, pw, fuk, feu, fw) = gradient(Forward,
                                                                 Supergrassi.consumption_weights,
                                                                 Const(el),
                                                                 puk,
                                                                 Const(peu),
                                                                 Const(pw),
                                                                 Const(fuk),
                                                                 Const(feu),
                                                                 Const(fw))

    for row in axes(alpha, 1)
        alpha[row,:] .= Supergrassi.consumption_weights(elasticity_a,
                                                      df.logP_uk[row],
                                                      df.logP_eu[row],
                                                      df.logP_w[row],
                                                      df.f_uk[row],
                                                      df.f_eu[row],
                                                      df.f_w[row])

        grad_alpha[row,:] .= gradient_function(elasticity_a,
                                               df.logP_uk[row],
                                               df.logP_eu[row],
                                               df.logP_w[row],
                                               df.f_uk[row],
                                               df.f_eu[row],
                                               df.f_w[row])[2]
    end
    
    alpha[12,2:6] .= 0.0
    grad_alpha[12,:] .= 0.0
    
    @test isapprox(alpha[:,1:3], [df.alpha_uk df.alpha_eu df.alpha_w], atol = tol)
    @test isapprox(grad_alpha[:,4:6], [df.dlogalpha_uk df.dlogalpha_eu df.dlogalpha_w], atol = tol)

    Alpha = Vector{Float64}(undef, n)

    logPf = Supergrassi.log_price_index.(elasticity_a,df.logP_uk,df.logP_eu,df.logP_w,df.f_uk,df.f_eu,df.f_w)
    Alpha = jacobian(ForwardWithPrimal, Supergrassi.total_consumption_weight, logPf, Const(df.f), Const(elasticity))

    @test isapprox(Alpha.val, df.alpha, atol = tol)

end

# @testset "Compute equilibrium test" begin

#     using Distributions
#     using Random
    
#     n = 16
#     log_consumer_price_index = zeros(n)
#     rand!(Normal(), log_consumer_price_index)
    
# end
