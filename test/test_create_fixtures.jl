using Supergrassi, CSV, DataFrames

tol = 1e-12

path = joinpath(@__DIR__, "..", "config","settings.yml")
settings_path = create_filepath(path)
settings = read_settings(settings_path)
filepaths = check_file_availability(settings)
data = read_data(filepaths, settings)

df_clean_data = CSV.read(joinpath(@__DIR__, "..", "data", "test_load_data.csv"), DataFrame)
df2d_clean_data = CSV.read(joinpath(@__DIR__, "..", "data", "test_load_data_2d.csv"), DataFrame)


df = outerjoin(CSV.read(joinpath(@__DIR__,"..","data", "data_for_household_demand.csv"), DataFrame),
               CSV.read(joinpath(@__DIR__, "..", "data", "data_for_eu_demand.csv"), DataFrame),
               CSV.read(joinpath(@__DIR__, "..", "data", "data_for_row_demand.csv"), DataFrame),
               CSV.read(joinpath(@__DIR__, "..", "data", "data_for_capital_production.csv"), DataFrame),
               on = [:logP_uk, :logP_eu, :logP_w], makeunique = true)

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
