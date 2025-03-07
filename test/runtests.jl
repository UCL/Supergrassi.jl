using Supergrassi
using Test
using Enzyme
using CSV, DataFrames

@testset "Supergrassi.jl" begin
    # Write your tests here.
end

@testset "Household demand parameters" begin

    n = 16
    elasticity_a = 2.0
    elasticity = 3.0
    tol = 1e-12
    
    alpha = Matrix{Float64}(undef, n, 6)
    grad_alpha = Matrix{Float64}(undef, n, 6)
    df = CSV.read("../data/data_for_household_demand.csv", DataFrame)

    for row in axes(alpha, 1)
        alpha[row,:] .= Supergrassi.consumption_weights(elasticity_a,
                                                      df.logP_uk[row],
                                                      df.logP_eu[row],
                                                      df.logP_w[row],
                                                      df.f_uk[row],
                                                      df.f_eu[row],
                                                      df.f_w[row])
        grad_alpha[row,:] .= gradient(Forward,
                                      Supergrassi.consumption_weights,
                                      Const(elasticity_a),
                                      df.logP_uk[row],
                                      Const(df.logP_eu[row]),
                                      Const(df.logP_w[row]),
                                      Const(df.f_uk[row]),
                                      Const(df.f_eu[row]),
                                      Const(df.f_w[row]))[2]
    end
    
    alpha[12,2:6] .= 0.0
    grad_alpha[12,:] .= 0.0
    
    @test isapprox(alpha[:,1], df.alpha_uk, atol = tol)
    @test isapprox(alpha[:,2], df.alpha_eu, atol = tol)
    @test isapprox(alpha[:,3], df.alpha_w, atol = tol)
    @test isapprox(grad_alpha[:,4], df.dlogalpha_uk, atol = tol)
    @test isapprox(grad_alpha[:,5], df.dlogalpha_eu, atol = tol)
    @test isapprox(grad_alpha[:,6], df.dlogalpha_w, atol = tol)

    Alpha = Vector{Float64}(undef, n)

    logPf = Supergrassi.log_price_index.(elasticity_a,df.logP_uk,df.logP_eu,df.logP_w,df.f_uk,df.f_eu,df.f_w)
    Alpha = jacobian(ForwardWithPrimal, Supergrassi.total_consumption_weight, logPf, Const(df.f), Const(elasticity))

    @test isapprox(Alpha.val, df.alpha, atol = tol)

end

@testset "EU demand parameters" begin

    n = 16
    elasticity_a = 2.0
    elasticity = 3.0
    tol = 1e-12
    df = CSV.read("../data/data_for_eu_demand.csv", DataFrame)

    beta = Matrix{Float64}(undef, n, 6)
    grad_beta = Matrix{Float64}(undef, n, 6)

    for row in axes(beta, 1)

        beta[row,:] .= Supergrassi.consumption_weights(elasticity_a,
                                                       df.logP_uk[row],
                                                       df.logP_eu[row],
                                                       df.logP_w[row],
                                                       df.x1_uk[row],
                                                       df.x1_eu[row],
                                                       df.x1_w[row])
        grad_beta[row,:] .= gradient(Forward,
                                     Supergrassi.consumption_weights,
                                     Const(elasticity_a),
                                     df.logP_uk[row],
                                     Const(df.logP_eu[row]),
                                     Const(df.logP_w[row]),
                                     Const(df.x1_uk[row]),
                                     Const(df.x1_eu[row]),
                                     Const(df.x1_w[row]))[2]


    end

    i = [4,     5,     7,     8,     9,    10,    11,    12,    13,    16]

    beta[i,1] .= 1.0
    beta[i,2:6] .= 0.0
    grad_beta[i,:] .= 0.0
    
    @test isapprox(beta[:,1], df.beta_uk, atol = tol)
    @test isapprox(beta[:,2], df.beta_eu, atol = tol)
    @test isapprox(beta[:,3], df.beta_w, atol = tol)
    @test isapprox(grad_beta[:,4], df.dlogbeta_uk, atol = tol)
    @test isapprox(grad_beta[:,5], df.dlogbeta_eu, atol = tol)
    @test isapprox(grad_beta[:,6], df.dlogbeta_w, atol = tol)
    
end
