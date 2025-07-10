using Supergrassi, Test


@testset "Unit Tests for Parameters" begin

    @testset "Parameters by region" begin
        @test Supergrassi.parameter_by_region(0, 0, 0, 0) == 0
        @test Supergrassi.parameter_by_region(0, 1, 5, 5) == 1

        @test Supergrassi.log_parameter_by_region(0, 0, 0, 0) == 0
        @test Supergrassi.log_parameter_by_region(0, 1, 5, 5) == log(Supergrassi.parameter_by_region(0, 1, 5, 5))
    end

end