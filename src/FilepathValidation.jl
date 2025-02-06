struct FilePath

    path::String
    directory::String
    file::String

    function FilePath(path::String)
        if ispathvalid(path)
            directory, file = splitpath(path)
            return new(path, directory, file)
        else
            throw(ArgumentError("Invalid path: $path"))
        end
    end

end

function create_filepath_from_template(basepath::String, substitution_dict::Dict{String, String})
    filename = replace(basepath, r"\{([^}]+)\}" => (m -> substitution_dict[m[1]]))
    return FilePath(filename)
end
