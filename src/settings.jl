using YAML

"""
    classify_dict(dict::Dict)

Classifies a dictionary type based on the types of its keys and values.
Used to create a typed dictionary from an untyped dictionary extracted from a YAML file.

# Arguments
- `dict::Dict{Any, Any}`: The dictionary to classify.

# Returns
- `Dict{key_type, value_type}`: The classified dictionary.

# Examples
```julia
untyped_dict = YAML.load_file("settings.yaml")
typed_dict = classify_dict(dict)
```
"""
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

"""
    read_settings(file_path::FilePath)

Reads settings from a YAML file and classifies the dictionary.

# Arguments
- `file_path::FilePath`: The path to the YAML file containing settings.

# Returns
- `Dict{key_type, value_type}`: The classified settings dictionary.
"""
function read_settings(file_path::FilePath)
    settings = YAML.load_file(file_path.path)
    settings = classify_dict(settings)
    return settings
end
