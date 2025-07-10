using Supergrassi, Test


@testset "Unit Tests for Parameters" begin

    @testset "Parameters by region" begin
        @test Supergrassi.parameter_by_region(0, 0, 0, 0) == 0
        @test Supergrassi.parameter_by_region(0, 1, 5, 5) == 1

        @test Supergrassi.log_parameter_by_region(0, 0, 0, 0) == 0
        @test Supergrassi.log_parameter_by_region(0, 1, 5, 5) == log(Supergrassi.parameter_by_region(0, 1, 5, 5))

        @test length(Supergrassi.parameters_by_region(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0)) == 3
        @test length(Supergrassi.log_parameters_by_region(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0)) == 3
    end

end