using Supergrassi, DataFrames, CSV, Enzyme, Test

tol = 1e-12

path = joinpath(@__DIR__, "..", "config","settings.yml")
settings_path = create_filepath(path)
settings = read_settings(settings_path)
filepaths = check_file_availability(settings)
data = read_data(filepaths, settings)
df = outerjoin(CSV.read(joinpath(@__DIR__,"..","data", "data_for_household_demand.csv"), DataFrame),
               CSV.read(joinpath(@__DIR__, "..", "data", "data_for_eu_demand.csv"), DataFrame),
               CSV.read(joinpath(@__DIR__, "..", "data", "data_for_row_demand.csv"), DataFrame),
               CSV.read(joinpath(@__DIR__, "..", "data", "data_for_capital_production.csv"), DataFrame),
               on = [:logP_uk, :logP_eu, :logP_w], makeunique = true)

prices = DataFrame([df.logP_uk, df.logP_eu, df.logP_w], ["uk", "eu", "world"])

clean = Supergrassi.clean_data(data,settings)
Supergrassi.postprocess_clean_data!(clean)

@testset "household demand parameters" begin

    params, ∂params = Supergrassi.compute_all_parameters(clean, prices)
    log_params, ∂log_params = Supergrassi.compute_all_parameters(clean, prices, Supergrassi.log_parameters_by_region)

    @test isapprox(params.consumption.uk, df.alpha_uk, atol = tol)
    @test isapprox(params.consumption.eu, df.alpha_eu, atol = tol)
    @test isapprox(params.consumption.world, df.alpha_w, atol = tol)
    @test isapprox(params.consumption.agg, df.alpha, atol = tol)

    @test isapprox(params.export_eu.uk, df.beta_uk, atol = tol)
    @test isapprox(params.export_eu.eu, df.beta_eu, atol = tol)
    @test isapprox(params.export_eu.world, df.beta_w, atol = tol)
    @test isapprox(params.export_eu.agg, df.beta, atol = tol)

    @test isapprox(params.export_world.uk, df.beta2_uk, atol = tol)
    @test isapprox(params.export_world.eu, df.beta2_eu, atol = tol)
    @test isapprox(params.export_world.world, df.beta2_w, atol = tol)
    @test isapprox(params.export_world.agg, df.beta2, atol = tol)

    @test isapprox(params.investment.uk, df.rho_uk, atol = tol)
    @test isapprox(params.investment.eu, df.rho_eu, atol = tol)
    @test isapprox(params.investment.world, df.rho_w, atol = tol)
    # @test isapprox(params.investment.agg, df.rho, atol = tol)

    @test isapprox(∂log_params.consumption.uk, df.dlogalpha_uk, atol = tol)
    @test isapprox(∂log_params.consumption.eu, df.dlogalpha_eu, atol = tol)
    @test isapprox(∂log_params.consumption.world, df.dlogalpha_w, atol = tol)

    @test isapprox(∂log_params.export_eu.uk, df.dlogbeta_uk, atol = tol)
    @test isapprox(∂log_params.export_eu.eu, df.dlogbeta_eu, atol = tol)
    @test isapprox(∂log_params.export_eu.world, df.dlogbeta_w, atol = tol)

    @test isapprox(∂log_params.export_world.uk, df.dlogbeta2_uk, atol = tol)
    @test isapprox(∂log_params.export_world.eu, df.dlogbeta2_eu, atol = tol)
    @test isapprox(∂log_params.export_world.world, df.dlogbeta2_w, atol = tol)

    @test isapprox(∂log_params.investment.uk, df.dlogrho_uk, atol = tol)
    @test isapprox(∂log_params.investment.eu, df.dlogrho_eu, atol = tol)
    @test isapprox(∂log_params.investment.world, df.dlogrho_w, atol = tol)

end
