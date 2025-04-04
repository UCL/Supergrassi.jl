module Supergrassi

# Write your package code here.
include("utility_function_parameters.jl")
include("excess_demand.jl")

include("FilePathValidator.jl")
include("FileReader.jl")
include("ReadSettings.jl")
include("Data.jl")

export create_filepath, create_filepath_from_template, read_data, read_settings, validate_settings, check_file_availability, organise_data

end
