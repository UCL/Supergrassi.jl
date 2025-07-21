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

df = CSV.read(joinpath(@__DIR__,"..","data", "excess_demand_terms.csv"), DataFrame)
log_prices = DataFrame([df.logP_uk, df.logP_eu, df.logP_w], ["uk", "eu", "world"])
price_uk = exp.(df.logP_uk)
price_eu = exp.(df.logP_eu)
price_world = exp.(df.logP_w)

params, ∂params = Supergrassi.compute_all_parameters(clean, log_prices)
#log_params, ∂log_params = Supergrassi.compute_all_parameters(clean, prices, Supergrassi.log_parameters_by_region)

operating_cost = df.zOC
# Copied from matlab code to validate results
household_expenditure = 7.27421398152353670952

logPf = Supergrassi.log_price_index(params.consumption.uk, params.consumption.eu, params.consumption.world,
                                    price_uk, price_eu, price_world,
                                    clean.constants.elasticities.consumption.armington)
logPBar = Supergrassi.log_agg_price_index(params.consumption.agg,
                                          logPf,
                                          clean.constants.elasticities.consumption.substitution)
logEf = Supergrassi.log_expenditure(params.consumption.agg,
                                    household_expenditure,
                                    clean.constants.elasticities.consumption.substitution,
                                    logPf, logPBar)
EF = Supergrassi.expenditure_by_region(params.consumption.uk, params.consumption.eu, params.consumption.world,
                                       price_uk, price_eu, price_world, logEf, clean.constants.elasticities.consumption)    

F = Supergrassi.market_clearing_price(price_uk, operating_cost, household_expenditure,
                                      price_eu, price_world, params, clean.industry, clean.constants)

@testset "Market Clearing" begin
    
    # Jac = jacobian(set_runtime_activity(ForwardWithPrimal),
    #                Supergrassi.market_clearing_price,
    #                price_uk,
    #                operating_cost,
    #                household_expenditure,
    #                Const(price_eu),
    #                Const(price_world),
    #                Const(params),
    #                Const(clean.industry),
    #                Const(clean.constants))

    @test isapprox(F[1][1:15], df.pdYBar[1:15], atol = tol)
    @test isapprox(F[2][1:15], df.EFd[1:15], atol = tol)
    @test isapprox(F[3][1:15], df.EX1d[1:15], atol = tol)
    @test isapprox(F[4][1:15], df.EX2d[1:15], atol = tol)
    @test isapprox(F[5][1:15], df.EId[1:15], atol = tol)
    @test isapprox(F[6][1:15], df.EMd[1:15], atol = tol)
    
    # F = Jac.val
    # ∂F_∂Pd = Jac.derivs[1]
    # ∂F_∂zOC = Jac.derivs[2]
    # ∂F_∂E = Jac.derivs[3]

    # @test length(F) == nrow(df)
    # @test ndims(∂F_∂Pd) == 2
    # @test ndims(∂F_∂zOC) == 2
    # @test ndims(∂F_∂E) == 1
    @show F

end
