using Supergrassi
using Test

using Aqua

Aqua.test_all(Supergrassi)

@testset "Supergrassi.jl" begin
    # Write your tests here.
end

include("test_filepath_creation.jl")
include("test_parameters.jl")
