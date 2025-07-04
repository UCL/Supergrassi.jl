using Supergrassi, DataFrames, CSV, Enzyme, Test

tol = 1e-12

# path = joinpath(@__DIR__, "..", "config","settings.yml")
# settings_path = create_filepath(path)
# settings = read_settings(settings_path)
# filepaths = check_file_availability(settings)
# data = read_data(filepaths, settings)

# df = outerjoin(CSV.read(joinpath(@__DIR__,"..","data", "data_for_household_demand.csv"), DataFrame),
#                CSV.read(joinpath(@__DIR__, "..", "data", "data_for_eu_demand.csv"), DataFrame),
#                CSV.read(joinpath(@__DIR__, "..", "data", "data_for_row_demand.csv"), DataFrame),
#                CSV.read(joinpath(@__DIR__, "..", "data", "data_for_capital_production.csv"), DataFrame),
#                on = [:logP_uk, :logP_eu, :logP_w], makeunique = true)

# df1d = CSV.read(joinpath(@__DIR__, "..", "data", "1d_data_for_firm_production.csv"), DataFrame)
# df2d = CSV.read(joinpath(@__DIR__, "..", "data", "2d_data_for_firm_production.csv"), DataFrame)

df1d = CSV.read(joinpath(@__DIR__, "..", "data", "parms_1d.csv"), DataFrame)
df2d = CSV.read(joinpath(@__DIR__, "..", "data", "parms_2d.csv"), DataFrame)

∂df1d = CSV.read(joinpath(@__DIR__, "..", "data", "dparms_1d.csv"), DataFrame)
∂df2d = CSV.read(joinpath(@__DIR__, "..", "data", "dparms_2d.csv"), DataFrame)

n = 16
gammaM_ref = reshape(df2d.gammaM, (n,n))
gammaMUK_ref = reshape(df2d.gammaMUK, (n,n))
gammaMEU_ref = reshape(df2d.gammaMEU, (n,n))
gammaMW_ref = reshape(df2d.gammaMW, (n,n))

prices = DataFrame([df.logP_uk, df.logP_eu, df.logP_w], ["uk", "eu", "world"])

clean = Supergrassi.clean_data(data,settings)
Supergrassi.postprocess_clean_data!(clean)

params, ∂params = Supergrassi.compute_all_parameters(clean, prices, false)
log_params, ∂log_params = Supergrassi.compute_all_parameters(clean, prices, true)

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

    @test isapprox(params.production.input_uk, reshape(df2d.gammaMUK, (n,n)), atol = tol)
    @test isapprox(params.production.input_eu, reshape(df2d.gammaMEU, (n,n)), atol = tol)
    @test isapprox(params.production.input_world, reshape(df2d.gammaMW, (n,n)), atol = tol)
    @test isapprox(params.production.input_agg, reshape(df2d.gammaM, (n,n)), atol = tol)
    @test isapprox(params.production.input_human, df1d.gammaH, atol = tol)
    @test isapprox(params.production.input_capital, df1d.gammaK, atol = tol)

    @test isapprox(params.production.shock_mean, df1d.mu, atol = tol)

end

@testset "Parameter derivatives" begin

    @test isapprox(∂log_params.consumption.uk, ∂df1d.logAlphaUK, atol = tol)
    @test isapprox(∂log_params.consumption.eu, ∂df1d.logAlphaEU, atol = tol)
    @test isapprox(∂log_params.consumption.world, ∂df1d.logAlphaW, atol = tol)

    @test isapprox(∂log_params.export_eu.uk, ∂df1d.logBeta1UK, atol = tol)
    @test isapprox(∂log_params.export_eu.eu, ∂df1d.logBeta1EU, atol = tol)
    @test isapprox(∂log_params.export_eu.world, ∂df1d.logBeta1W, atol = tol)

    @test isapprox(∂log_params.export_world.uk, ∂df1d.logBeta2UK, atol = tol)
    @test isapprox(∂log_params.export_world.eu, ∂df1d.logBeta2EU, atol = tol)
    @test isapprox(∂log_params.export_world.world, ∂df1d.logBeta2W, atol = tol)

    @test isapprox(∂log_params.investment.uk, ∂df1d.logRhoUK, atol = tol)
    @test isapprox(∂log_params.investment.eu, ∂df1d.logRhoEU, atol = tol)
    @test isapprox(∂log_params.investment.world, ∂df1d.logRhoW, atol = tol)

    @test isapprox(∂log_params.production.input_human, reshape(∂df2d.logGammaH, (n,n)), atol = tol)
    @test isapprox(∂log_params.production.input_capital, reshape(∂df2d.logGammaK, (n,n)), atol = tol)
    @test isapprox(∂log_params.production.shock_mean, reshape(∂df2d.logMu, (n,n)), atol = tol)
    
end
