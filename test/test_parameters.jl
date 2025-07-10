using Supergrassi, Test


@testset "Unit Tests for Parameters" begin

    @test Supergrassi.parameter_by_region(0, 0, 0, 0) == 0

end