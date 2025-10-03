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

log_price_uk = df.logP_uk
log_price_eu = df.logP_eu
log_price_world = df.logP_w
price_uk = exp.(df.logP_uk)
price_eu = exp.(df.logP_eu)
price_world = exp.(df.logP_w)

params = Supergrassi.compute_all_parameters(clean, log_price_uk, log_price_eu, log_price_world)

operating_cost = clean.industry.surplus.val;
household_expenditure = clean.industry.regional.totals.expenditure;
muI = Supergrassi.compute_muI(clean.industry, params.constants.elasticities.investment);

x = deepcopy([log_price_uk;
             operating_cost;
             household_expenditure;
             muI;
             clean.industry.depreciation.val])

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

OC = Supergrassi.compute_operating_cost_constraint(x, log_price_eu, log_price_world, clean, params)

@testset "Operating Cost" begin

    @test isequal(length(OC), 16)

end

CFC = Supergrassi.compute_fixed_capital_consumption_constraint(x, rand(16), clean.industry, params)

@testset "Fixed Capital" begin

    @test isequal(length(CFC), 16)

end

Jac1 = jacobian(set_runtime_activity(ForwardWithPrimal),
                Supergrassi.market_clearing_price_constraint,
                price_uk, operating_cost, household_expenditure,
                Const(price_eu),
                Const(price_world),
                Const(muI),
                Const(params),
                Const(clean.industry),
                Const(clean.constants))

Jac2 = jacobian(set_runtime_activity(ForwardWithPrimal),
                Supergrassi.compute_operating_cost_constraint,
                x,
                Const(log_price_eu),
                Const(log_price_world),
                Const(clean),
                Const(params))

@testset "Jacobian" begin

    @test length(Jac1.val) == 6
    @test isapprox(sum(Jac1.val), df.pdYBar + df.EFd + df.EX1d + df.EX2d + df.EId + df.EMd, atol = tol)
    @test length(Jac1.derivs[1]) == 16
    @test length(Jac1.derivs[2]) == 16
    @test length(Jac1.derivs[3]) == 6

    @test length(Jac2.val) == 16
    # TODO: Improve these tests
    @test isless(sum(Jac2.val), 4.0)
    @test isless(3.9, sum(Jac2.val))
    @test length(Jac2.derivs[1]) == 16*50
    
end
