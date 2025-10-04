module Supergrassi

# Write your package code here.

include("types.jl")
include("pathvalidator.jl")
include("reader.jl")
include("settings.jl")
include("data.jl")

include("utils.jl")
include("remap_industries.jl")
include("cleanup.jl")
include("post_processing.jl")

include("parameters.jl")
include("parameters_interface.jl")
include("excess_demand.jl")
include("capital_market.jl")

include("objective_function.jl")
include("constraint_function.jl")

include("minimisation.jl")
include("optimisation_helpers.jl")

include("estimation.jl")


export create_filepath, read_data, read_settings, check_file_availability
export clean_data
export estimate

export compute_objective_function

end
