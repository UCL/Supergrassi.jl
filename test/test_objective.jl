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

params, ∂params = Supergrassi.compute_all_parameters(clean, prices, false)
log_params, ∂log_params = Supergrassi.compute_all_parameters(clean, prices, true)

objective_function = Supergrassi.compute_objective_function(df.logP_uk, clean, params)

@testset "Objective function" begin

    target_obj = 13.58236775806153850965
    @test isapprox(objective_function, target_obj, atol = tol)

end
