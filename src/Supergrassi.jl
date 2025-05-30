module Supergrassi

# Write your package code here.

include("VariableStructure.jl")
include("FilePathValidator.jl")
include("FileReader.jl")
include("ReadSettings.jl")
include("Data.jl")
include("Cleanup.jl")

export create_filepath, read_data, read_settings, check_file_availability
export CleanData

end
