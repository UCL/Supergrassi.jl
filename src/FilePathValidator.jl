"""

    FilePathValidator

    A structure for file paths.

    # Fields
    - path::String: The full path.
    - directory::String: The directory.
    - file::String: The file name.
    - extension::String: The file extension.
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
