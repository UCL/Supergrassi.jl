using Supergrassi
using Test
using Enzyme
using CSV, DataFrames

# @testset "Supergrassi.jl" begin
#     # Write your tests here.
# end

n = 16
elasticity_a = 2.0
elasticity = 3.0
elasticity_tilde = 2.0
tol = 1e-12

@testset "Household demand parameters" begin

    alpha = Matrix{Float64}(undef, n, 3)
    grad_log_alpha = Matrix{Float64}(undef, n, 3)
    df = CSV.read("../data/data_for_household_demand.csv", DataFrame)

    for row in axes(alpha, 1)
        alpha[row,:] .= Supergrassi.parameters_by_region(elasticity_a,
                                                      df.logP_uk[row],
                                                      df.logP_eu[row],
                                                      df.logP_w[row],
                                                      df.f_uk[row],
                                                      df.f_eu[row],
                                                      df.f_w[row])
        grad_log_alpha[row,:] .= gradient(Forward,
                                          Supergrassi.log_parameters_by_region,
                                          Const(elasticity_a),
                                          df.logP_uk[row],
                                          Const(df.logP_eu[row]),
                                          Const(df.logP_w[row]),
                                          Const(df.f_uk[row]),
                                          Const(df.f_eu[row]),
                                          Const(df.f_w[row]))[2]
    end

    mask = 12

    alpha[mask,2:3] .= 0.0
    grad_log_alpha[mask,:] .= 0.0

    @test isapprox(alpha[:,1], df.alpha_uk, atol = tol)
    @test isapprox(alpha[:,2], df.alpha_eu, atol = tol)
    @test isapprox(alpha[:,3], df.alpha_w, atol = tol)
    @test isapprox(grad_log_alpha[:,1], df.dlogalpha_uk, atol = tol)
    @test isapprox(grad_log_alpha[:,2], df.dlogalpha_eu, atol = tol)
    @test isapprox(grad_log_alpha[:,3], df.dlogalpha_w, atol = tol)

    logPf = Supergrassi.log_price_index.(elasticity_a,df.logP_uk,df.logP_eu,df.logP_w,df.f_uk,df.f_eu,df.f_w)
    Alpha = jacobian(ForwardWithPrimal, Supergrassi.total_parameters, logPf, Const(df.f), Const(elasticity))

    @test isapprox(Alpha.val, df.alpha, atol = tol)

end

@testset "EU demand parameters" begin

    df = CSV.read("../data/data_for_eu_demand.csv", DataFrame)

    beta = Matrix{Float64}(undef, n, 3)
    grad_log_beta = Matrix{Float64}(undef, n, 3)
    grad_beta_tilde = Vector{Float64}(undef, n)

    Ex = 1.3953e+00
    ETilde = 2.5449e+01
    PTilde = 1.0
    fx_EUR = 1.1733e+00
    logBeta1Tilde = 5.5323e-01

    for row in axes(beta, 1)

        beta[row,:] .= Supergrassi.parameters_by_region(elasticity_a,
                                                       df.logP_uk[row],
                                                       df.logP_eu[row],
                                                       df.logP_w[row],
                                                       df.x1_uk[row],
                                                       df.x1_eu[row],
                                                       df.x1_w[row])
        grad_log_beta[row,:] .= gradient(Forward,
                                     Supergrassi.log_parameters_by_region,
                                     Const(elasticity_a),
                                     df.logP_uk[row],
                                     Const(df.logP_eu[row]),
                                     Const(df.logP_w[row]),
                                     Const(df.x1_uk[row]),
                                     Const(df.x1_eu[row]),
                                     Const(df.x1_w[row]))[2]

    end

    logPf = Supergrassi.log_price_index.(elasticity_a,df.logP_uk,df.logP_eu,df.logP_w,df.x1_uk,df.x1_eu,df.x1_w)
    Beta = jacobian(ForwardWithPrimal, Supergrassi.total_parameters, logPf, Const(df.x1), Const(elasticity))
    grad_log_beta_tilde = gradient(ForwardWithPrimal,
                                   Supergrassi.log_eu_expenditure_on_uk_exports,
                                   logPf,
                                   Const(df.x1),
                                   Const(Ex),
                                   Const(ETilde),
                                   Const(fx_EUR),
                                   Const(PTilde),
                                   Const(elasticity),
                                   Const(elasticity_tilde))

    mask = [4, 5, 7, 8, 9, 10, 11, 12, 13, 16]

    beta[mask,1] .= 1.0
    beta[mask,2:3] .= 0.0
    grad_log_beta[mask,:] .= 0.0

    @test isapprox(beta[:,1], df.beta_uk, atol = tol)
    @test isapprox(beta[:,2], df.beta_eu, atol = tol)
    @test isapprox(beta[:,3], df.beta_w, atol = tol)
    @test isapprox(grad_log_beta[:,1], df.dlogbeta_uk, atol = tol)
    @test isapprox(grad_log_beta[:,2], df.dlogbeta_eu, atol = tol)
    @test isapprox(grad_log_beta[:,3], df.dlogbeta_w, atol = tol)
    @test isapprox(Beta.val, df.beta, atol = tol)
    @test isapprox(grad_log_beta_tilde.val, logBeta1Tilde, atol = 1e-4)

end

@testset "Rest of World demand parameters" begin

    df = CSV.read("../data/data_for_row_demand.csv", DataFrame)

    beta = Matrix{Float64}(undef, n, 3)
    grad_log_beta = Matrix{Float64}(undef, n, 3)

    Ex = 1.4894e+00
    ETilde = 2.3429e+01
    PTilde = 1.0
    fx_USD = 1.2345e+00
    logBeta2Tilde = 7.4242e-01

    for row in axes(beta, 1)

        beta[row,:] .= Supergrassi.parameters_by_region(elasticity_a,
                                                       df.logP_uk[row],
                                                       df.logP_eu[row],
                                                       df.logP_w[row],
                                                       df.x2_uk[row],
                                                       df.x2_eu[row],
                                                       df.x2_w[row])
        grad_log_beta[row,:] .= gradient(Forward,
                                     Supergrassi.log_parameters_by_region,
                                     Const(elasticity_a),
                                     df.logP_uk[row],
                                     Const(df.logP_eu[row]),
                                     Const(df.logP_w[row]),
                                     Const(df.x2_uk[row]),
                                     Const(df.x2_eu[row]),
                                     Const(df.x2_w[row]))[2]


    end

    logPf = Supergrassi.log_price_index.(elasticity_a,df.logP_uk,df.logP_eu,df.logP_w,df.x2_uk,df.x2_eu,df.x2_w)
    Beta = jacobian(ForwardWithPrimal, Supergrassi.total_parameters, logPf, Const(df.x2), Const(elasticity))
    grad_log_beta2_tilde = gradient(ForwardWithPrimal,
                                   Supergrassi.log_eu_expenditure_on_uk_exports,
                                   logPf,
                                   Const(df.x2),
                                   Const(Ex),
                                   Const(ETilde),
                                   Const(fx_USD),
                                   Const(PTilde),
                                   Const(elasticity),
                                   Const(elasticity_tilde))

    mask = [4, 5, 7, 8, 9, 10, 11, 12, 13, 16]

    beta[mask,1] .= 1.0
    beta[mask,2:3] .= 0.0
    grad_log_beta[mask,:] .= 0.0

    @test isapprox(beta[:,1], df.beta2_uk, atol = tol)
    @test isapprox(beta[:,2], df.beta2_eu, atol = tol)
    @test isapprox(beta[:,3], df.beta2_w, atol = tol)
    @test isapprox(grad_log_beta[:,1], df.dlogbeta2_uk, atol = tol)
    @test isapprox(grad_log_beta[:,2], df.dlogbeta2_eu, atol = tol)
    @test isapprox(grad_log_beta[:,3], df.dlogbeta2_w, atol = tol)
    @test isapprox(Beta.val, df.beta2, atol = tol)
    @test isapprox(grad_log_beta2_tilde.val, logBeta2Tilde, atol = 1e-4)

end

@testset "Capital production function parameters" begin

    df = CSV.read("../data/data_for_capital_production.csv", DataFrame)

    elasticity = 0.4

    rho = Matrix{Float64}(undef, n, 3)
    grad_log_rho = Matrix{Float64}(undef, n, 3)

    for row in axes(rho, 1)

        rho[row,:] .= Supergrassi.parameters_by_region(elasticity_a,
                                                      df.logP_uk[row],
                                                      df.logP_eu[row],
                                                      df.logP_w[row],
                                                      df.I_uk[row],
                                                      df.I_eu[row],
                                                      df.I_w[row])
        grad_log_rho[row,:] .= gradient(Forward,
                                    Supergrassi.log_parameters_by_region,
                                    Const(elasticity_a),
                                    df.logP_uk[row],
                                    Const(df.logP_eu[row]),
                                    Const(df.logP_w[row]),
                                    Const(df.I_uk[row]),
                                    Const(df.I_eu[row]),
                                    Const(df.I_w[row]))[2]

    end

    mask_zero_uk = [2, 6, 7, 11, 13, 14, 16]
    mask_one_uk = [5, 8, 9, 12, 15]
    mask_zero_row = [2, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16]

    mask = Vector{Bool}(undef, n)
    mask .= true
    mask[mask_zero_uk] .= false

    rho[mask_zero_uk,1] .= 0.0
    rho[mask_one_uk,1]  .= 1.0
    rho[mask_zero_row, 2:3] .= 0.0

    grad_log_rho[mask_zero_uk, 1] .= 0.0
    grad_log_rho[mask_one_uk, 1]  .= 0.0
    grad_log_rho[mask_zero_row, 2:3] .= 0.0

    @test isapprox(rho[:,1], df.rho_uk, atol = tol)
    @test isapprox(rho[:,2], df.rho_eu, atol = tol)
    @test isapprox(rho[:,3], df.rho_w, atol = tol)
    @test isapprox(grad_log_rho[:,1], df.dlogrho_uk, atol = tol)
    @test isapprox(grad_log_rho[:,2], df.dlogrho_eu, atol = tol)
    @test isapprox(grad_log_rho[:,3], df.dlogrho_w, atol = tol)

    logPI = Supergrassi.log_price_index.(elasticity_a,df.logP_uk,df.logP_eu,df.logP_w,df.I_uk,df.I_eu,df.I_w)
    Rho = jacobian(ForwardWithPrimal, Supergrassi.total_parameters, logPI[mask], Const(df.I[mask]), Const(elasticity))

    Rho_out = zeros(n)
    Rho_out[mask] .= Rho.val

    @test isapprox(Rho_out, df.rho, atol = tol)

end

@testset "Production function parameters" begin

    xi = 0.4
    xi_a = 2.0

    df = CSV.read("../data/1d_data_for_firm_production.csv", DataFrame)
    df2d = CSV.read("../data/2d_data_for_firm_production.csv", DataFrame)

    gammaM_ref = reshape(df2d.gammaM, (n,n))
    gammaMUK_ref = reshape(df2d.gammaMUK, (n,n))
    gammaMEU_ref = reshape(df2d.gammaMEU, (n,n))
    gammaMW_ref = reshape(df2d.gammaMW, (n,n))

    m = reshape(df2d.m, (n,n))
    mUK = reshape(df2d.mUK, (n,n))
    mEU = reshape(df2d.mEU, (n,n))
    mW = reshape(df2d.mW, (n,n))

    gammaM = zeros(n,n,3)
    grad_gammaM = zeros(n,n,3)

    mask = fill(true, n,n)
    ind = [12 16 96 112 208 209 224 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255 256]
    mask[ind] .= false

    for i in axes(gammaMUK_ref, 1)
        for j in axes(gammaMUK_ref, 2)
            if (mask[i,j])
                gammaM[i,j,:] .= Supergrassi.parameters_by_region(xi_a, df.logP_uk[j], df.logP_eu[j], df.logP_w[j],
                                                                 mUK[i,j], mEU[i,j], mW[i,j])
                grad_gammaM[i,j,:] .= gradient(Forward,
                                               Supergrassi.parameters_by_region,
                                               Const(xi_a),
                                               df.logP_uk[j],
                                               Const(df.logP_eu[j]),
                                               Const(df.logP_w[j]),
                                               Const(mUK[i,j]),
                                               Const(mEU[i,j]),
                                               Const(mW[i,j]))[2]

            end
        end
    end


    @test isapprox(gammaM[:,:,1], gammaMUK_ref, atol = tol)
    @test isapprox(gammaM[:,:,2], gammaMEU_ref, atol = tol)
    @test isapprox(gammaM[:,:,3], gammaMW_ref, atol = tol)

    GammaM = zeros(n,n)
    GammaH = zeros(n)
    GammaK = zeros(n)
    mu = zeros(n)
    grad_GammaM = zeros(n,n,n)
    grad_GammaH = zeros(n,n)
    grad_GammaK = zeros(n,n)

    for i in axes(m,1)

        logPm = Supergrassi.log_price_index.(xi_a,df.logP_uk,df.logP_eu,df.logP_w,mUK[i,:],mEU[i,:],mW[i,:])
        logPm[isinf.(logPm)] .= 0.0

        jacM = jacobian(ForwardWithPrimal,
                        Supergrassi.total_input_parameters,
                        logPm,
                        Const(m[i,:]),
                        Const(df.k[i]),
                        Const(df.k0[i]),
                        Const(df.y[i]),
                        Const(df.h[i]),
                        Const(df.logW[i]),
                        Const(xi),
                        Const(df.tau[i]))

        GammaM[i,:] .= jacM.val
        grad_GammaM[i,:,:] .= jacM.derivs[1]

        jacH = jacobian(ForwardWithPrimal,
                        Supergrassi.total_labor_parameters,
                        logPm,
                        Const(m[i,:]),
                        Const(df.k[i]),
                        Const(df.k0[i]),
                        Const(df.y[i]),
                        Const(df.h[i]),
                        Const(df.logW[i]),
                        Const(xi),
                        Const(df.tau[i]))

        GammaH[i] = jacH.val
        grad_GammaH[i,:] .= jacH.derivs[1]

        jacK = jacobian(ForwardWithPrimal,
                        Supergrassi.total_capital_parameters,
                        logPm,
                        Const(m[i,:]),
                        Const(df.k[i]),
                        Const(df.k0[i]),
                        Const(df.y[i]),
                        Const(df.h[i]),
                        Const(df.logW[i]),
                        Const(xi),
                        Const(df.tau[i]))

        GammaK[i] = jacK.val
        grad_GammaK[i,:] .= jacK.derivs[1]

        mu[i] = Supergrassi.productivity_shock_mean(xi, df.logP_uk[i], logPm, m[i,:], df.k[i], df.k0[i],
                                                    df.y[i], df.h[i], df.logW[i], df.tau[i])
        
    end

    @test isapprox(GammaM, gammaM_ref, atol = tol)
    @test isapprox(GammaH, df.gammaH, atol = tol)
    @test isapprox(GammaK, df.gammaK, atol = tol)

    @test isapprox(mu, df.mu, atol = tol)
end

@testset "Intermediate goods price index" begin

    xi = 0.4
    df = CSV.read("../data/data_for_goods_price_index.csv", DataFrame)

    pdYBar = Supergrassi.intermediate_goods_price_index(df.logP_uk, df.zOC, df.tau, df.mu, df.gammaK, df.K0, xi)

    @test isapprox(pdYBar, df.pdYBar, atol = tol)

end

include("test_filepath_creation.jl")
