module Supergrassi

# Write your package code here.

include("types.jl")
include("pathvalidator.jl")
include("reader.jl")
include("settings.jl")
include("data.jl")
include("cleanup.jl")

include("parameters.jl")
include("parameters_interface.jl")
include("excess_demand.jl")

function estimate()
    @info "Estimation started."

    path = joinpath(@__DIR__, "..", "config","settings.yml")
    settings_path = create_filepath(path)
    settings = read_settings(settings_path)
    filepaths = check_file_availability(settings)
    data = read_data(filepaths, settings)

    @info "Data read successfully."
    
    clean = Supergrassi.clean_data(data,settings)
    Supergrassi.postprocess_clean_data!(clean)

    @info "Data cleaned and post-processed."

    df = CSV.read(joinpath(@__DIR__,"..","data", "data_for_household_demand.csv"), DataFrame)

    prices = DataFrame([df.logP_uk, df.logP_eu, df.logP_w], ["uk", "eu", "world"])

    @info "Prices extracted from data."

    params, ∂params = Supergrassi.compute_all_parameters(clean, prices, false)
    log_params, ∂log_params = Supergrassi.compute_all_parameters(clean, prices, true)

    @info "Parameters computed."

    @info "Estimation completed."

    return settings, data, clean, params, ∂params, log_params, ∂log_params
end

export create_filepath, read_data, read_settings, check_file_availability
export clean_data
export estimate

end
