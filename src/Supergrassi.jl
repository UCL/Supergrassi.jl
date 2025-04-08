module Supergrassi

# Write your package code here.

include("FilePathValidator.jl")
include("FileReader.jl")
include("ReadSettings.jl")
include("Data.jl")
include("Cleanup.jl")

export create_filepath, read_data, read_settings, check_file_availability
export cleanup, create_map_105_to_64

end
