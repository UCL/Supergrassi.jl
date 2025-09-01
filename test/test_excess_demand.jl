using Test
using Supergrassi
using CSV, DataFrames, Enzyme

tol = 1e-8

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

log_price_uk = df.logP_uk
log_price_eu = df.logP_eu
log_price_world = df.logP_w
price_uk = exp.(df.logP_uk)
price_eu = exp.(df.logP_eu)
price_world = exp.(df.logP_w)

params = Supergrassi.compute_all_parameters(clean, log_price_uk, log_price_eu, log_price_world)
#log_params = Supergrassi.compute_all_parameters(clean, prices, Supergrassi.log_parameters_by_region)

operating_cost = df.zOC
# Copied from matlab code to validate results
household_expenditure = 7.27421398152353670952
muI = Supergrassi.compute_muI(clean.industry, params.constants.elasticities.investment);

    
F = Supergrassi.market_clearing_price_constraint(price_uk, operating_cost, household_expenditure,
                                                 price_eu, price_world, muI, params, clean.industry, clean.constants)

@testset "Market Clearing" begin

    @test isapprox(F[1][1:15], df.pdYBar[1:15], atol = tol)
    @test isapprox(F[2][1:15], df.EFd[1:15], atol = tol)
    @test isapprox(F[3][1:15], df.EX1d[1:15], atol = tol)
    @test isapprox(F[4][1:15], df.EX2d[1:15], atol = tol)
    @test isapprox(F[5][1:15], df.EId[1:15], atol = tol)
    @test isapprox(F[6][1:15], df.EMd[1:15], atol = tol)

end


Jac = jacobian(set_runtime_activity(ReverseWithPrimal),
               Supergrassi.market_clearing_price_constraint,
               price_uk, operating_cost, household_expenditure,
               Const(price_eu),
               Const(price_world),
               Const(muI),
               Const(params),
               Const(clean.industry),
               Const(clean.constants))

@testset "Jacobian" begin

    @test length(Jac.val) == 6
    @test isapprox(sum(Jac.val), df.pdYBar + df.EFd + df.EX1d + df.EX2d + df.EId + df.EMd, atol = tol)
    @test length(Jac.derivs[1]) == 16
    @test length(Jac.derivs[2]) == 16
    @test length(Jac.derivs[3]) == 6
    
#    @test length(Jac.val) == nrow(df)
#    @test ndims(Jac.derivs[1]) == 2
#    @test ndims(Jac.derivs[2]) == 2
#    @test ndims(Jac.derivs[3]) == 1

end
