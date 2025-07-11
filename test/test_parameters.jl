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

    @testset "Total parameters" begin
        
        log_price_index = [1.0, 2.0, 3.0]
        quantity = [4.0, 5.0, 6.0]
        elas = 12.3

        @test length(Supergrassi.total_parameters(log_price_index, quantity, elas)) == 3

        try
            Supergrassi.total_parameters(log_price_index, quantity[1:2], elas)
        catch e
            @test isa(e, ErrorException)
            @test occursin("log_price_index and quantity must have the same length", e.msg)

        end
    end

    @testset "Log EU Expenditure on UK Exports" begin
        log_price_index = [1.0, 2.0, 3.0]
        quantity = [4.0, 5.0, 6.0]
        Ex = 1.2
        ETilde = 0.8
        ePx = 0.5
        PTilde = 0.3
        elas = 12.3
        elasticity_tilde = 0.1

        @test isa(Supergrassi.log_eu_expenditure_on_uk_exports(log_price_index, quantity, Ex, ETilde, ePx, PTilde, elas, elasticity_tilde), Float64)

        try
            Supergrassi.log_eu_expenditure_on_uk_exports(log_price_index, quantity[1:2], Ex, ETilde, ePx, PTilde, elas, elasticity_tilde)
        catch e
            @test isa(e, ErrorException)
            @test occursin("log_price_index and quantity must have the same length", e.msg)
        end
    end

    @testset "Price Index" begin
        
        @test isa(Supergrassi.price_index(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0), Float64)
        @test isa(Supergrassi.log_price_index(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0), Float64)
        @test Supergrassi.price_index(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0) == exp(Supergrassi.log_price_index(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0))
        @test Supergrassi.log_price_index(1.0, 2.0, 3.0, 4.0, 0.0, 0.0, 0.0) == 0.0
        @test isa(Supergrassi.log_total_price_index(1.0, [2.0, 3.0], [4.0, 5.0, 6.0]), Float64)

    end

end