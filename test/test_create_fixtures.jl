using Supergrassi, CSV, DataFrames

tol = 1e-12

config_path = joinpath(@__DIR__, "..", "config","settings.yml")
settings_path = create_filepath(config_path)
settings = read_settings(settings_path)
filepaths = check_file_availability(settings)
data = read_data(filepaths, settings)

data_path = joinpath(@__DIR__, "..", "data")

df_clean_data = CSV.read(joinpath(data_path, "test_load_data.csv"), DataFrame)
df2d_clean_data = CSV.read(joinpath(data_path, "test_load_data_2d.csv"), DataFrame)


df = outerjoin(CSV.read(joinpath(data_path, "data_for_household_demand.csv"), DataFrame),
               CSV.read(joinpath(data_path, "data_for_eu_demand.csv"), DataFrame),
               CSV.read(joinpath(data_path, "data_for_row_demand.csv"), DataFrame),
               CSV.read(joinpath(data_path, "data_for_capital_production.csv"), DataFrame),
               on = [:logP_uk, :logP_eu, :logP_w], makeunique = true)

df1d = CSV.read(joinpath(data_path, "parms_1d.csv"), DataFrame)
df2d = CSV.read(joinpath(data_path, "parms_2d.csv"), DataFrame)

∂df1d = CSV.read(joinpath(data_path, "dparms_1d.csv"), DataFrame)
∂df2d = CSV.read(joinpath(data_path, "dparms_2d.csv"), DataFrame)

prices = DataFrame([df.logP_uk, df.logP_eu, df.logP_w], ["uk", "eu", "world"])

clean = Supergrassi.clean_data(data,settings)
Supergrassi.postprocess_clean_data!(clean)

params, _ = Supergrassi.compute_all_parameters(clean, prices, false)
_, ∂log_params = Supergrassi.compute_all_parameters(clean, prices, true)
