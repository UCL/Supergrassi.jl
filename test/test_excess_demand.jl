using Test
using Supergrassi
using CSV, DataFrames, Enzyme

tol = 1e-12

if (!@isdefined data)
    settings_path = create_filepath("config/settings.yml")
    settings = read_settings(settings_path)
    filepaths = check_file_availability(settings)
    data = read_data(filepaths, settings)
end

if (!@isdefined clean)
    clean = Supergrassi.clean_data(data,settings)
    Supergrassi.postprocess_clean_data!(clean)
end

df = outerjoin(CSV.read(joinpath(@__DIR__,"..","data", "data_for_household_demand.csv"), DataFrame),
               CSV.read(joinpath(@__DIR__, "..", "data", "data_for_eu_demand.csv"), DataFrame),
               CSV.read(joinpath(@__DIR__, "..", "data", "data_for_row_demand.csv"), DataFrame),
               CSV.read(joinpath(@__DIR__, "..", "data", "data_for_capital_production.csv"), DataFrame),
               on = [:logP_uk, :logP_eu, :logP_w], makeunique = true)
log_prices = DataFrame([df.logP_uk, df.logP_eu, df.logP_w], ["uk", "eu", "world"])
prices = DataFrame([exp.(df.logP_uk), exp.(df.logP_eu), exp.(df.logP_w)], ["uk", "eu", "world"])

params, ∂params = Supergrassi.compute_all_parameters(clean, log_prices)
#log_params, ∂log_params = Supergrassi.compute_all_parameters(clean, prices, Supergrassi.log_parameters_by_region)

operating_cost = randn(nrow(prices))
household_expenditure = 7.2742 # Copied from matlab code to avoid numerical issues

@testset "Market Clearing" begin

    logPf = Supergrassi.log_price_index(params.consumption,
                                        prices,
                                        clean.constants.elasticities.consumption.armington)

    logPBar = Supergrassi.log_agg_price_index(params.consumption.agg,
                                              logPf,
                                              clean.constants.elasticities.consumption.substitution)

    logEf = Supergrassi.log_expenditure(params.consumption.agg,
                                        household_expenditure,
                                        clean.constants.elasticities.consumption.substitution,
                                        logPf, logPBar)

    F = Supergrassi.market_clearing_price(prices.uk, operating_cost, household_expenditure,
                                          prices.eu, prices.world, params, clean.industry, clean.constants)

    # Jac = jacobian(ForwardWithPrimal,
    #                Supergrassi.market_clearing_price,
    #                prices.uk,
    #                operating_cost,
    #                household_expenditure,
    #                Const(prices.eu),
    #                Const(prices.world),
    #                Const(params),
    #                Const(clean.industry),
    #                Const(clean.constants))

    # F = Jac.val
    # ∂F_∂Pd = Jac.derivs[1]
    # ∂F_∂zOC = Jac.derivs[2]
    # ∂F_∂E = Jac.derivs[3]

    @test length(F) == nrow(prices)
    @show F

end
