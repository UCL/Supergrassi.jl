using YAML

function classify_dict(dict::Dict{Any, Any})

    dict_keys = collect(keys(dict))
    key_types = typeof.(values(dict_keys))

    dict_values = collect(values(dict))
    value_types = typeof.(dict_values)

    set_key_types = Set(key_types)
    set_value_types = Set(value_types)

    key_type = length(set_key_types) > 1 ? Any : first(set_key_types)
    value_type = length(set_value_types) > 1 ? Any : first(set_value_types)

    new_dict = Dict{key_type, value_type}(dict)

    for (key, value) in dict
        if isa(value, Dict)
            new_dict[key] = classify_dict(value)
        end
    end

    return new_dict

end


function read_settings(file_path::FilePath)
    settings = YAML.load_file(file_path.path)
    settings = classify_dict(settings)
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