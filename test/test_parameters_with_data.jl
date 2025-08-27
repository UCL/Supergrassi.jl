using Supergrassi, DataFrames, CSV, Enzyme, Test

n = 16

df1d = CSV.read(joinpath(data_path, "parms_1d.csv"), DataFrame)
df2d = CSV.read(joinpath(data_path, "parms_2d.csv"), DataFrame)

∂df1d = CSV.read(joinpath(data_path, "dparms_1d.csv"), DataFrame)
∂df2d = CSV.read(joinpath(data_path, "dparms_2d.csv"), DataFrame)

#_, ∂log_params = Supergrassi.compute_all_parameters(clean, prices, true)

@testset "Parameter values" begin

    @test isapprox(params.consumption.uk, df1d.alphaUK, atol = tol)
    @test isapprox(params.consumption.eu, df1d.alphaEU, atol = tol)
    @test isapprox(params.consumption.world, df1d.alphaW, atol = tol)
    @test isapprox(params.consumption.agg, df1d.alpha, atol = tol)

    @test isapprox(params.export_eu.uk, df1d.beta1UK, atol = tol)
    @test isapprox(params.export_eu.eu, df1d.beta1EU, atol = tol)
    @test isapprox(params.export_eu.world, df1d.beta1W, atol = tol)
    @test isapprox(params.export_eu.agg, df1d.beta1, atol = tol)

    @test isapprox(params.export_world.uk, df1d.beta2UK, atol = tol)
    @test isapprox(params.export_world.eu, df1d.beta2EU, atol = tol)
    @test isapprox(params.export_world.world, df1d.beta2W, atol = tol)
    @test isapprox(params.export_world.agg, df1d.beta2, atol = tol)

    @test isapprox(params.investment.uk, df1d.rhoUK, atol = tol)
    @test isapprox(params.investment.eu, df1d.rhoEU, atol = tol)
    @test isapprox(params.investment.world, df1d.rhoW, atol = tol)
    @test isapprox(params.investment.agg, df1d.rho, atol = tol)

    @test isapprox(params.production.uk, reshape(df2d.gammaMUK, (n,n)), atol = tol)
    @test isapprox(params.production.eu, reshape(df2d.gammaMEU, (n,n)), atol = tol)
    @test isapprox(params.production.world, reshape(df2d.gammaMW, (n,n)), atol = tol)
    @test isapprox(params.production.agg, reshape(df2d.gammaM, (n,n)), atol = tol)
    @test isapprox(params.production.human, df1d.gammaH, atol = tol)
    @test isapprox(params.production.capital, df1d.gammaK, atol = tol)

    @test isapprox(params.production.shock_mean, df1d.mu, atol = tol)

    @test isapprox(params.production.low_skill, df1d.gammaHL, atol = tol)
    @test isapprox(params.production.high_skill, df1d.gammaHH, atol = tol)

end

∂α = Matrix{Float64}(undef, n, 3)
∂β1 = Matrix{Float64}(undef, n, 3)
∂β2 = Matrix{Float64}(undef, n, 3)
∂ρ = Matrix{Float64}(undef, n, 3)

∂γH = Matrix{Float64}(undef, n, n)
∂γK = Matrix{Float64}(undef, n, n)
∂μ = Matrix{Float64}(undef, n, n)

logW = Supergrassi.compute_wage_index(clean.household, clean.constants.elasticities.production)
tau = Supergrassi.compute_advalorem_tax(clean.industry)

for row = 1:n
    ∂α[row,:] .= gradient(Forward,
                          Supergrassi.log_parameters_by_region,
                          Const(clean.constants.elasticities.consumption.armington),
                          price_uk[row],
                          Const(price_eu[row]),
                          Const(price_world[row]),
                          Const(clean.industry.regional.consumption.uk[row]),
                          Const(clean.industry.regional.consumption.eu[row]),
                          Const(clean.industry.regional.consumption.world[row]))[2]

    ∂β1[row,:] .= gradient(Forward,
                          Supergrassi.log_parameters_by_region,
                          Const(clean.constants.elasticities.eu_export_demand.armington),
                          price_uk[row],
                          Const(price_eu[row]),
                          Const(price_world[row]),
                          Const(clean.industry.regional.export_eu.uk[row]),
                          Const(clean.industry.regional.export_eu.eu[row]),
                          Const(clean.industry.regional.export_eu.world[row]))[2]

    ∂β2[row,:] .= gradient(Forward,
                          Supergrassi.log_parameters_by_region,
                          Const(clean.constants.elasticities.world_export_demand.armington),
                          price_uk[row],
                          Const(price_eu[row]),
                          Const(price_world[row]),
                          Const(clean.industry.regional.export_world.uk[row]),
                          Const(clean.industry.regional.export_world.eu[row]),
                          Const(clean.industry.regional.export_world.world[row]))[2]

    ∂ρ[row,:] .= gradient(Forward,
                          Supergrassi.log_parameters_by_region,
                          Const(clean.constants.elasticities.investment.armington),
                          price_uk[row],
                          Const(price_eu[row]),
                          Const(price_world[row]),
                          Const(clean.industry.regional.investment.uk[row]),
                          Const(clean.industry.regional.investment.eu[row]),
                          Const(clean.industry.regional.investment.world[row]))[2]

    ∂γH[row,:] = gradient(Forward,
                          Supergrassi.total_labor_parameters,
                          price_uk,
                          Const(price_eu),
                          Const(price_world),
                          Const(clean.industry.regional.input_matrices.uk[row,:]),
                          Const(clean.industry.regional.input_matrices.eu[row,:]),
                          Const(clean.industry.regional.input_matrices.world[row,:]),
                          Const(clean.industry.regional.input_matrices.agg[row,:]),
                          Const(clean.industry.surplus.val[row]),
                          Const(clean.industry.capital.current_year[row]),
                          Const(clean.industry.regional.total_use.agg[row]),
                          Const(clean.household.payments.agg[row]),
                          Const(logW[row]),
                          Const(clean.constants.elasticities.production),
                          Const(tau[row]),
                          Const(true))[1]

    ∂γK[row,:] = gradient(Forward,
                          Supergrassi.total_capital_parameters,
                          price_uk,
                          Const(price_eu),
                          Const(price_world),
                          Const(clean.industry.regional.input_matrices.uk[row,:]),
                          Const(clean.industry.regional.input_matrices.eu[row,:]),
                          Const(clean.industry.regional.input_matrices.world[row,:]),
                          Const(clean.industry.regional.input_matrices.agg[row,:]),
                          Const(clean.industry.surplus.val[row]),
                          Const(clean.industry.capital.current_year[row]),
                          Const(clean.industry.regional.total_use.agg[row]),
                          Const(clean.household.payments.agg[row]),
                          Const(logW[row]),
                          Const(clean.constants.elasticities.production),
                          Const(tau[row]),
                          Const(true))[1]

    ∂μ[row,:] = gradient(Forward,
                         Supergrassi.productivity_shock_mean,
                         Const(clean.constants.elasticities.production),
                         price_uk,
                         Const(price_eu),
                         Const(price_world),
                         Const(clean.industry.regional.input_matrices.uk[row,:]),
                         Const(clean.industry.regional.input_matrices.eu[row,:]),
                         Const(clean.industry.regional.input_matrices.world[row,:]),
                         Const(clean.industry.regional.input_matrices.agg[row,:]),
                         Const(clean.industry.surplus.val[row]),
                         Const(clean.industry.capital.current_year[row]),
                         Const(clean.industry.regional.total_use.agg[row]),
                         Const(clean.household.payments.agg[row]),
                         Const(logW[row]),
                         Const(tau[row]),
                         Const(row),
                         Const(true))[2]



end

@testset "Parameter derivatives" begin

    @test isapprox(∂α[:,1], ∂df1d.logAlphaUK, atol = tol)
    @test isapprox(∂α[:,2], ∂df1d.logAlphaEU, atol = tol)
    @test isapprox(∂α[:,3], ∂df1d.logAlphaW, atol = tol)

    @test isapprox(∂β1[:,1], ∂df1d.logBeta1UK, atol = tol)
    @test isapprox(∂β1[:,2], ∂df1d.logBeta1EU, atol = tol)
    @test isapprox(∂β1[:,3], ∂df1d.logBeta1W, atol = tol)

    @test isapprox(∂β2[:,1], ∂df1d.logBeta2UK, atol = tol)
    @test isapprox(∂β2[:,2], ∂df1d.logBeta2EU, atol = tol)
    @test isapprox(∂β2[:,3], ∂df1d.logBeta2W, atol = tol)

    @test isapprox(∂ρ[:,1], ∂df1d.logRhoUK, atol = tol)
    @test isapprox(∂ρ[:,2], ∂df1d.logRhoEU, atol = tol)
    @test isapprox(∂ρ[:,3], ∂df1d.logRhoW, atol = tol)

    @test isapprox(∂γH, reshape(∂df2d.logGammaH, (n,n)), atol = tol)
    @test isapprox(∂γK, reshape(∂df2d.logGammaK, (n,n)), atol = tol)
    @test isapprox(∂μ, reshape(∂df2d.logMu, (n,n)), atol = tol)

end
