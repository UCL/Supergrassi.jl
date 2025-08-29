using Supergrassi, Test


@testset "End to End" begin
    results = Supergrassi.estimate()

    @test !isnothing(results)

end
