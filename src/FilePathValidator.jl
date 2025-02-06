struct FilePath

    path::String
    directory::String
    file::String

    function FilePath(path::String)

        directory = dirname(path)
        file = basename(path)
        new(path, directory, file)

    end

end

function create_filepath_from_template(basepath::String, substitution_dict::Dict{String, Any})

    for (key, value) in substitution_dict
        var = "{$key}"
        basepath = replace(basepath, var => value)
    end

    return FilePath(basepath)
end
