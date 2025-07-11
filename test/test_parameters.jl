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

    prices_uk = [1.0, 2.0, 3.0]
    prices_eu = [4.0, 5.0, 6.0]
    prices_world = [7.0, 8.0, 9.0]

    input_uk = [10.0, 11.0, 12.0]
    input_eu = [13.0, 14.0, 15.0]
    input_world = [16.0, 17.0, 18.0]
    input_agg = [19.0, 20.0, 21.0]

    surplus = 22.0
    capital = 23.0
    output = 24.0
    labour = 25.0
    low_wages = 26.0
    elasticity = Supergrassi.Elasticity(0.5, 0.5, 0.5, 0.5)
    tau = 0.1

    @testset "Total Input Parameters" begin
        @test isa(Supergrassi.total_input_parameters(prices_uk, prices_eu, prices_world, input_uk, input_eu, input_world, input_agg, surplus, capital, output, labour, low_wages, elasticity, tau, false), Vector{Float64})
        @test isa(Supergrassi.total_input_parameters(prices_uk, prices_eu, prices_world, input_uk, input_eu, input_world, input_agg, surplus, capital, output, labour, low_wages, elasticity, tau, true), Vector{Float64})
        @test !any(isinf, Supergrassi.total_input_parameters(prices_uk, prices_eu, prices_world, input_uk, input_eu, input_world, input_agg, surplus, capital, output, labour, low_wages, elasticity, tau, true))

        try
            Supergrassi.total_input_parameters(prices_uk, prices_eu, prices_world, input_uk[1:2], input_eu[1:2], input_world[1:2], input_agg[1:2], surplus, capital, output, labour, low_wages, elasticity, tau, true)
        catch e
            @test isa(e, ErrorException)
            @test occursin("logPm and input_agg must have the same length", e.msg)
        end
    end

    @testset "Total Labour Parameters" begin
        @test isa(Supergrassi.total_labor_parameters(prices_uk, prices_eu, prices_world, input_uk, input_eu, input_world, input_agg, surplus, capital, output, labour, low_wages, elasticity, tau, false), Float64)
        @test isa(Supergrassi.total_labor_parameters(prices_uk, prices_eu, prices_world, input_uk, input_eu, input_world, input_agg, surplus, capital, output, labour, low_wages, elasticity, tau, true), Float64)
        @test !any(isinf, Supergrassi.total_labor_parameters(prices_uk, prices_eu, prices_world, input_uk, input_eu, input_world, input_agg, surplus, capital, output, labour, low_wages, elasticity, tau, true))

        try
            Supergrassi.total_labor_parameters(prices_uk, prices_eu, prices_world, input_uk[1:2], input_eu[1:2], input_world[1:2], input_agg[1:2], surplus, capital, output, labour, low_wages, elasticity, tau, true)
        catch e
            @test isa(e, ErrorException)
            @test occursin("logPm and input_agg must have the same length", e.msg)
        end
    end

    @testset "Total Capital Parameters" begin
        @test isa(Supergrassi.total_capital_parameters(prices_uk, prices_eu, prices_world, input_uk, input_eu, input_world, input_agg, surplus, capital, output, labour, low_wages, elasticity, tau, false), Float64)
        @test isa(Supergrassi.total_capital_parameters(prices_uk, prices_eu, prices_world, input_uk, input_eu, input_world, input_agg, surplus, capital, output, labour, low_wages, elasticity, tau, true), Float64)
        @test !any(isinf, Supergrassi.total_capital_parameters(prices_uk, prices_eu, prices_world, input_uk, input_eu, input_world, input_agg, surplus, capital, output, labour, low_wages, elasticity, tau, true))

        try
            Supergrassi.total_capital_parameters(prices_uk, prices_eu, prices_world, input_uk[1:2], input_eu[1:2], input_world[1:2], input_agg[1:2], surplus, capital, output, labour, low_wages, elasticity, tau, true)
        catch e
            @test isa(e, ErrorException)
            @test occursin("logPm and input_agg must have the same length", e.msg)
        end
    end

    @testset "Productivity Shock Mean" begin
        @test isa(Supergrassi.productivity_shock_mean(elasticity, prices_uk, prices_eu, prices_world, input_uk, input_eu, input_world, input_agg, surplus, capital, output, labour, low_wages, tau, 1, false), Float64)
        @test isa(Supergrassi.productivity_shock_mean(elasticity, prices_uk, prices_eu, prices_world, input_uk, input_eu, input_world, input_agg, surplus, capital, output, labour, low_wages, tau, 1, true), Float64)
        @test !any(isinf, Supergrassi.productivity_shock_mean(elasticity, prices_uk, prices_eu, prices_world, input_uk, input_eu, input_world, input_agg, surplus, capital, output, labour, low_wages, tau, 1, true))
    end

    @testset "Mathematical Functions" begin
        @test isa(Supergrassi.weight_kernel(1.1, 2.2, 3.3), Float64)
        @test isa(Supergrassi.log_weight_kernel(1.1, 2.2, 3.3), Float64)
        @test isa(Supergrassi.sum_kernel([1.2, 3.4], [5.6, 6.7], 7.8), Float64)
        @test isa(Supergrassi.tauPdMu(4.2, [1.2, 3.4], [5.6, 6.7], 7.8, 1.4, 2.8, 3.2, 4.6, 0.5), Float64)
        @test isa(Supergrassi.capital_fun(1.1, 0.2, 3.3, 4.4, 5.5), Float64)
        @test isa(Supergrassi.labor_fun(1.1, 2.2, 3.3), Float64)
    end

end