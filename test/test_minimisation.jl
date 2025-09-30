using Test
using Supergrassi
using CSV, DataFrames

tol = 1e-12

# Load test fixtures if not already defined
if (!@isdefined data)
    settings_path = create_filepath("config/settings.yml")
    settings = read_settings(settings_path)
    filepaths = check_file_availability(settings)
    data = read_data(filepaths, settings)
end

if (!@isdefined clean)
    clean = Supergrassi.clean_data(data, settings)
    Supergrassi.postprocess_clean_data!(clean)
end

# Load test data
df = CSV.read(joinpath(@__DIR__, "..", "data", "data_for_household_demand.csv"), DataFrame)

log_price_uk = df.logP_uk
log_price_eu = df.logP_eu
log_price_world = df.logP_w

# Create test input vector
zOC = clean.industry.surplus.val
x = vcat(log_price_uk, zOC)

@testset "Minimisation Functions" begin
    
    @testset "compute_gradient function" begin
        
        @testset "Basic functionality" begin
            # Test that gradient computation doesn't error
            gradient_result = zeros(Float64, length(x))
            
            @test_nowarn Supergrassi.compute_gradient(x, clean, log_price_eu, log_price_world, gradient_result)
            
            # Test that gradient is computed (not all zeros)
            @test !all(iszero, gradient_result)
            
            # Test that gradient has correct length
            @test length(gradient_result) == length(x)
            
            # Test that gradient elements are finite
            @test all(isfinite, gradient_result)
        end
    end

    @testset "Integration with objective function" begin
        # Test that gradient is related to objective function
        gradient_result = zeros(Float64, length(x))
        Supergrassi.compute_gradient(x, clean, log_price_eu, log_price_world, gradient_result)
        
        # Compute objective function at current point
        obj_current = Supergrassi.compute_objective_function(x, clean, log_price_eu, log_price_world)
        
        # Test finite difference approximation (basic check)
        epsilon = 1e-6
        x_plus = copy(x)
        x_plus[1] += epsilon
        obj_plus = Supergrassi.compute_objective_function(x_plus, clean, log_price_eu, log_price_world)
        
        finite_diff_grad = (obj_plus - obj_current) / epsilon
        
        # The analytical gradient should be close to finite difference
        # (allowing for numerical errors)
        @test abs(gradient_result[1] - finite_diff_grad) < 1e-3
    end
end
