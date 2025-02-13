module Supergrassi

# Write your package code here.

include("FilePathValidator.jl")
include("FileReader.jl")
include("ReadSettings.jl")

export create_filepath, create_filepath_from_template, read_data, read_settings, validate_settings, check_file_availability

end
