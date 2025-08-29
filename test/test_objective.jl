using Supergrassi, Test

zOC = clean.industry.surplus.val
x = vcat(df.logP_uk, zOC)
objective_function = Supergrassi.compute_objective_function(x, clean, df.logP_eu, df.logP_w)

@testset "Objective function" begin

    target_obj = 13.58236775806153850965
    @test isapprox(objective_function, target_obj, atol = tol)

end
