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

x = deepcopy([log_price_uk;
              clean.industry.surplus.val;
              clean.industry.regional.totals.savings;
              Supergrassi.compute_muI(clean.industry, params.constants.elasticities.investment);
              clean.industry.depreciation.val])

y = zeros(32)

CEQ = Supergrassi.constraint_wrapper(x, log_price_eu, log_price_world, params, clean.industry, clean.constants, y)
DCEQ = Supergrassi.compute_constraint_function(x, log_price_eu, log_price_world, clean, params) 

@testset "Constraint Function" begin

    @test length(CEQ) == 32
    @test length(DCEQ) == 32*50
    
end

