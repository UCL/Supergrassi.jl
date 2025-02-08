using YAML

function read_settings(file_path::FilePath)
    settings = YAML.load_file(file_path.path)

    settings_keys = collect(keys(settings))
    for key in settings_keys
        if !isa(key, String)
            error("Invalid key: $key")
        end
    end

    settings = Dict{String, Any}(settings)
    return settings
end

function validate_dictionary_keys(settings::Dict{String, Any}, required_keys::Vector{String})
    dict_keys = collect(keys(settings))
    for key in required_keys
        if !(key in dict_keys)
            error("Missing key: $key")
        end
    end
end



function validate_settings(settings::Dict{String, Any})

    base_keys = ["constants", "flags", "elasticities", "files", "initial_params", "experiment_run", "version"]
    validate_dictionary_keys(settings, base_keys)

    println("Settings are valid.")
end