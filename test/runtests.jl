using Supergrassi
using Test
using Aqua

Aqua.test_all(Supergrassi)

@testset "Supergrassi.jl" begin
    # Write your tests here.
end

include("test_filepath_creation.jl")
include("test_assets_cleaning.jl")
include("test_weighted_mean.jl")
if(isdir(joinpath(@__DIR__,"..","input", "uk_data")))
    include("test_clean_data.jl")
end
