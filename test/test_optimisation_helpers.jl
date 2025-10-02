using Test
using Supergrassi

@testset "Optimisation Helpers Tests" begin
    
    @testset "unpack_x function" begin
        n = 3
        x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0]
        
        log_price_uk, zOC, expenditure, log_muI, log_Delta = Supergrassi.unpack_x(n, x)
        
        @test length(log_price_uk) == n
        @test length(zOC) == n  
        @test isa(expenditure, Real)
        @test isa(log_muI, Real)
        @test length(log_Delta) == n
        
        @test log_price_uk == [1.0, 2.0, 3.0]
        @test zOC == [4.0, 5.0, 6.0]
        @test expenditure == 7.0
        @test log_muI == 8.0
        @test log_Delta == [9.0, 10.0, 11.0]
    end
    
    @testset "unpack_x edge cases" begin
        n = 1
        x = [1.0, 2.0, 3.0, 4.0, 5.0]
        
        log_price_uk, zOC, expenditure, log_muI, log_Delta = Supergrassi.unpack_x(n, x)
        
        @test log_price_uk == [1.0]
        @test zOC == [2.0]
        @test expenditure == 3.0
        @test log_muI == 4.0
        @test log_Delta == [5.0]
        
        n = 2
        x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
        
        log_price_uk, zOC, expenditure, log_muI, log_Delta = Supergrassi.unpack_x(n, x)
        
        @test log_price_uk == [1.0, 2.0]
        @test zOC == [3.0, 4.0]
        @test expenditure == 5.0
        @test log_muI == 6.0
        @test log_Delta == [7.0, 8.0]
    end
    
    @testset "unpack_x type preservation" begin
        n = 2
        x = Float32[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
        
        log_price_uk, zOC, expenditure, log_muI, log_Delta = Supergrassi.unpack_x(n, x)
        
        @test eltype(log_price_uk) == Float32
        @test eltype(zOC) == Float32
        @test isa(expenditure, Float32)
        @test isa(log_muI, Float32)
        @test eltype(log_Delta) == Float32
    end
    
    @testset "unpack_x input validation" begin
        n = 2
        wrong_size = 3 * n + 1  # One element short
        x_wrong = ones(wrong_size)
        
        @test_throws ErrorException Supergrassi.unpack_x(n, x_wrong)
        
        # Test correct size works
        correct_size = 3 * n + 2
        x_correct = ones(correct_size)
        
        @test_nowarn Supergrassi.unpack_x(n, x_correct)
    end

end
