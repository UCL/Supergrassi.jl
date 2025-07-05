using Logging
"""
    FilePath(path::String)

A structure for file paths.

# Fields
- `path::String`: The full path.
- `directory::String`: The directory.
- `file::String`: The file name.
- `extension::String`: The file extension.
"""
struct FilePath

    path::String
    directory::String
    file::String
    extension::String

    function FilePath(path::String)

        directory = dirname(path)
        file = basename(path)
        extension = splitext(file)[2]
        new(path, directory, file, extension)

    end

end

"""
    create_filepath_from_template(basepath::String, substitution_dict::Dict{String, Any})

Creates a file path from a template string and a dictionary of substitutions.

# Arguments
- `basepath::String`: The template string.
- `substitution_dict::Dict{String, Any}`: The dictionary of substitutions.

# Returns
- `FilePath`: The file path object.

# Examples
```julia
substitutions = Dict(
    "name" => "data",
    "format" => "csv",
    "directory" => "dir",
    "suffix" => 1,
)

template = "{directory}/{name}{suffix}.{format}"
filepath = create_filepath_from_template(template, substitutions)
println(filepath.path) # "dir/data1.csv"
println(filepath.directory) # "dir"
println(filepath.file) # "data1.csv"
println(filepath.extension) # ".csv"
```
"""
function create_filepath_from_template(basepath::String, substitution_dict::Dict{String, Any})

    for (key, value) in substitution_dict
        var = "{$key}"
        basepath = replace(basepath, var => value)
    end

    return FilePath(basepath)
end


"""
    create_filepath(basepath::String)

Creates a file path.

# Arguments
- `file_path::FilePath`: The file path object.

# Returns
- `Bool`: Whether the file path is valid.

# Examples
```julia
file_path = FilePath("data.csv")
println(is_valid) # true
```
"""
function create_filepath(basepath::String)

    return FilePath(basepath)

end

"""
    check_file_availability(settings::Dict)

Verifies the availability of data files.

# Arguments
- `settings::Dict{<:Any, <:Any}`: The settings dictionary.

# Returns
- `Dict{String, FilePath}`: A dictionary of file paths.

# Examples
```julia
settings = Dict(
    "files" => Dict(
        "input_dir" => "data",
        "input_derived_dir" => "derived",
        "inputs" => Dict(
            "base" => Dict(
                "data" => "data.csv"
            ),
            "derived" => Dict(
                "scenario_independent" => Dict(
                    "data" => "data.csv"
                ),
                "scenario_dependent" => Dict(
                    "data" => "data.csv"
                )
            )
        )
    ),
    "version" => "v1",
    "files" => Dict(
        "scenario" => "1"
    )
)

filepaths = check_file_availability(settings)
```
"""
function check_file_availability(settings::Dict{<:Any, <:Any})
    base_path = joinpath(@__DIR__, "..", settings["files"]["input_dir"])
    input_derived_dir = settings["files"]["input_derived_dir"]
    version = settings["version"]
    scenario = string(settings["files"]["scenario"])

    filepaths = Dict{String, FilePath}()

    @info "Verifying data from $base_path"

    # Define the categories and paths
    sources = [
        ("base", settings["files"]["inputs"]["base"], base_path),
        ("scenario_independent", settings["files"]["inputs"]["derived"]["scenario_independent"], joinpath(base_path, input_derived_dir, version)),
        ("scenario_dependent", settings["files"]["inputs"]["derived"]["scenario_dependent"], joinpath(base_path, input_derived_dir, version))
    ]

    for (category, files, path_prefix) in sources
        for (key, value) in files
            @debug "Verifying datafiles for $key"

            # Handle scenario-dependent files by modifying the filename
            if category == "scenario_dependent"
                value_parts = splitext(value)
                value = value_parts[1] * scenario * value_parts[2]
            end

            source_path = joinpath(path_prefix, value)

            if !isfile(source_path)
                @warn "File not found: $source_path...$(category == "scenario_dependent" ? " skipping" : "")"
            else
                filepaths[key] = create_filepath(source_path)
            end
        end
    end

    return filepaths
end
