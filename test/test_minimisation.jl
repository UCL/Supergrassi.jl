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
            gradient_result = zeros(Float64, length(x))
            
            @test_nowarn Supergrassi.compute_gradient(x, clean, log_price_eu, log_price_world, gradient_result)            
            @test !all(iszero, gradient_result)
            @test length(gradient_result) == length(x)
            @test all(isfinite, gradient_result)
        end
    end

    @testset "Integration with objective function" begin
        gradient_result = zeros(Float64, length(x))
        Supergrassi.compute_gradient(x, clean, log_price_eu, log_price_world, gradient_result)
        obj_current = Supergrassi.compute_objective_function(x, clean, log_price_eu, log_price_world)
        epsilon = 1e-6
        x_plus = copy(x)
        x_plus[1] += epsilon
        obj_plus = Supergrassi.compute_objective_function(x_plus, clean, log_price_eu, log_price_world)        
        finite_diff_grad = (obj_plus - obj_current) / epsilon
        @test abs(gradient_result[1] - finite_diff_grad) < 1e-3
    end
end
