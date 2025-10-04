using Supergrassi, CSV, DataFrames, Test

tol = 1e-6

config_path = joinpath(@__DIR__, "..", "config","settings.yml")
settings_path = create_filepath(config_path)
settings = read_settings(settings_path)
filepaths = check_file_availability(settings)
data = read_data(filepaths, settings)

data_path = joinpath(@__DIR__, "..", "data")
df = outerjoin(CSV.read(joinpath(data_path, "data_for_household_demand.csv"), DataFrame),
               CSV.read(joinpath(data_path, "data_for_eu_demand.csv"), DataFrame),
               CSV.read(joinpath(data_path, "data_for_row_demand.csv"), DataFrame),
               CSV.read(joinpath(data_path, "data_for_capital_production.csv"), DataFrame),
               on = [:logP_uk, :logP_eu, :logP_w], makeunique = true)

# prices = DataFrame([df.logP_uk, df.logP_eu, df.logP_w], ["uk", "eu", "world"])

log_price_uk = df.logP_uk
log_price_eu = df.logP_eu
log_price_world = df.logP_w
price_uk = exp.(df.logP_uk)
price_eu = exp.(df.logP_eu)
price_world = exp.(df.logP_w)


clean = Supergrassi.clean_data(data,settings)
Supergrassi.postprocess_clean_data!(clean)

params = Supergrassi.compute_all_parameters(clean, log_price_uk, log_price_eu, log_price_world, false)

@testset "Fixtures Typecheck" begin

    @test isa(clean, Supergrassi.CleanData)
    @test isa(params, Supergrassi.Parameters)
    @test isa(df, DataFrame)
    @test isa(settings, Dict{String, Any})
    @test isa(filepaths, Dict{String, Supergrassi.FilePath})
    @test isa(data, Supergrassi.Data)
    @test isa(data_path, String)
    @test isa(config_path, String)
    @test isa(settings_path, Supergrassi.FilePath)
    @test isa(tol, Float64)

end

