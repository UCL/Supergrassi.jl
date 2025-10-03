using Test
using Supergrassi
using CSV, DataFrames, Enzyme

x = deepcopy([price_uk;
              clean.industry.surplus.val;
              clean.industry.regional.totals.expenditure;
              Supergrassi.compute_muI(clean.industry, params.constants.elasticities.investment);
              clean.industry.depreciation.val])

y = zeros(48)

CEQ = Supergrassi.constraint_function(x, log_price_eu, log_price_world, clean, params, y)
DCEQ = Supergrassi.constraint_jacobian(x, log_price_eu, log_price_world, clean, params) 

@testset "Constraint Function" begin

    @test length(CEQ) == 48
    @test length(DCEQ) == 48*50
    
end

