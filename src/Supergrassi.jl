module Supergrassi

# Write your package code here.

include("types.jl")
include("pathvalidator.jl")
include("reader.jl")
include("settings.jl")
include("data.jl")
include("cleanup.jl")

include("parameters.jl")
include("parameters_interface.jl")
include("excess_demand.jl")
# include("capital_market.jl")
include("objective_function.jl")

export create_filepath, read_data, read_settings, check_file_availability
export clean_data

export compute_objective_function

end
