using Supergrassi, DataFrames, CSV, Enzyme, Test

tol = 1e-12

path = joinpath(@__DIR__, "..", "config","settings.yml")
settings_path = create_filepath(path)
settings = read_settings(settings_path)
filepaths = check_file_availability(settings)
data = read_data(filepaths, settings)
df = CSV.read(joinpath(@__DIR__,"..","data", "data_for_household_demand.csv"), DataFrame)
prices = DataFrame([df.logP_uk, df.logP_eu, df.logP_w], ["uk", "eu", "world"])

clean = Supergrassi.clean_data(data,settings)
Supergrassi.postprocess_clean_data!(clean)

@testset "household demand parameters" begin

    α, _ = Supergrassi.compute_parameter(clean.industry.regional.consumption,
                                         clean.constants.elasticities.consumption,
                                         prices)
    _, ∂logα = Supergrassi.compute_parameter(clean.industry.regional.consumption,
                                             clean.constants.elasticities.consumption,
                                             prices, Supergrassi.log_parameters_by_region)
    
    @test isapprox(α.uk, df.alpha_uk, atol = tol)
    @test isapprox(α.eu, df.alpha_eu, atol = tol)
    @test isapprox(α.world, df.alpha_w, atol = tol)
    @test isapprox(α.agg, df.alpha, atol = tol)

    @test isapprox(∂logα.uk, df.dlogalpha_uk, atol = tol)
    @test isapprox(∂logα.eu, df.dlogalpha_eu, atol = tol)
    @test isapprox(∂logα.world, df.dlogalpha_w, atol = tol)

end
