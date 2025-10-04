using Supergrassi, Test

df = CSV.read(joinpath(@__DIR__,"..","data", "excess_demand_terms.csv"), DataFrame)

price_uk = exp.(df.logP_uk)
operating_cost = df.zOC

KL, KD, FCF = Supergrassi.compute_capital_market(price_uk, operating_cost, clean.industry, params)

@testset "Capital Market" begin

    @test isapprox(Supergrassi.Delta(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0), 1.8591409142295225)
    @test isapprox(Supergrassi.Delta(0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0), 0.0)
    @test isapprox(Supergrassi.B(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0), 0.0)
    @test isapprox(Supergrassi.B(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5), 0.3824981438992028)
    @test isapprox(Supergrassi.b(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0), 0.0)
    @test isapprox(Supergrassi.b(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5), 0.7550813375962907)

    # The outcome of capital market depends on random numbers so we won't compare exactly to the
    # reference Matlab implementation. Check that the results are within sensible bounds.
    @test isless(maximum(KL), 2.0)
    @test isapprox(minimum(KL), 0.0)
    @test isless(maximum(KD), 4.0)
    @test isless(0.0, minimum(KD))
    
end
