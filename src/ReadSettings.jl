using YAML

data_type_dict = Dict(
    "int" => Int64,
    "float" => Float64,
    "str" => String,
    "bool" => Bool,
    "elas" => Vector{Vector{Float64}},
    "date" => String,
    "list" => Vector{Int64},
)
"""

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


function read_settings(file_path::FilePath)
    settings = YAML.load_file(file_path.path)
    settings = classify_dict(settings)
    return settings
end

function validate_settings(settings::Dict{<:Any, <:Any}, structure::Dict{<:Any, <:Any})

    for key in keys(structure)
        if !haskey(settings, key)
            error("Missing key in settings YAML: $key")
        end
    end

    for (key, value) in settings
        
        if !haskey(structure, key)
            @warn "Unvalidated Key in settings YAML: $key"
            continue
        end

        if isa(value, Dict)
            validate_settings(value, structure[key])
        else
            if !isa(value, data_type_dict[structure[key]])
                error("Invalid value type for key: $key")
            end
        end

    end
    
end
