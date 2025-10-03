using Supergrassi
using Test
using Aqua
using YAML

Aqua.test_all(Supergrassi)  

@testset "Supergrassi.jl" begin
    # Write your tests here.
end

include("test_filepath_creation.jl")
include("test_assets_cleaning.jl")
include("test_weighted_mean.jl")
include("test_round_shares.jl")

include("test_parameters.jl")

settings = YAML.load_file(joinpath(@__DIR__, "..", "config", "settings.yml"))
if(isdir(joinpath(@__DIR__,"..",settings["files"]["input_dir"])))

    include("test_create_fixtures.jl")

    include("test_clean_data.jl")
    include("test_parameters_with_data.jl")
    include("test_objective.jl")
    include("test_constraint.jl")
    include("test_minimisation.jl")
    include("test_optimisation_helpers.jl")
    include("test_capital_market.jl")
    
    include("test_estimation.jl")
end
