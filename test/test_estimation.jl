using Supergrassi, Test, Random

Random.seed!(1235)

@testset "End to End" begin
    results = Supergrassi.estimate()

    @test !isnothing(results)

end
