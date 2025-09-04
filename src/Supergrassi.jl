module Supergrassi

# Write your package code here.

include("types.jl")
include("pathvalidator.jl")
include("reader.jl")
include("settings.jl")
include("data.jl")

include("utils.jl")
include("remap_industries.jl")
include("cleanup.jl")
include("post_processing.jl")

include("parameters.jl")
include("parameters_interface.jl")
include("excess_demand.jl")
include("objective_function.jl")
include("constraint_function.jl")

include("minimisation.jl")

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

    df = CSV.read(joinpath(@__DIR__, "..", "data", "data_for_household_demand.csv"), DataFrame)


    log_prices_uk = df.logP_uk
    log_prices_eu = df.logP_eu
    log_prices_world = df.logP_w

    @info "Prices extracted from data."

    params = Supergrassi.compute_all_parameters(clean, log_prices_uk, log_prices_eu, log_prices_world, false)
    log_params = Supergrassi.compute_all_parameters(clean, log_prices_uk, log_prices_eu, log_prices_world, true)

    @info "Parameters computed."

    # TODO:replace with real household final consumption
    household_final_consumption = [1.2]

    x = deepcopy([log_prices_uk; clean.industry.surplus.val; clean.industry.regional.totals.savings; household_final_consumption; clean.industry.depreciation.val])
    println("Starting minimisation with x of length ", length(x))


    obj_val = compute_objective_function(x, clean, log_prices_eu, log_prices_world)
    gradient = compute_gradient(x, clean, log_prices_eu, log_prices_world)

    constraint_value = constraint_wrapper(x, log_prices_eu, log_prices_world, params, clean.industry, clean.constants)
    jacobian = compute_constraint_function(x, log_prices_eu, log_prices_world, clean, params)


    prob = Ipopt.CreateIpoptProblem(
        16,
        [0.0 for i in 1:16],
        [10.0 for i in 1:16],
        49,
        [0.0 for i in 1:49],
        [10.0 for i in 1:49],
        8,
        10,
        compute_objective_function,
        constraint_wrapper,
        compute_gradient,
        compute_constraint_function,
        [0 for i in 1:16],
    )

    @info "Estimation completed."

    # return settings, data, clean, params, log_params, gradient
    return settings, data, clean, params, log_params, obj_val, gradient, constraint_value, jacobian
end

export create_filepath, read_data, read_settings, check_file_availability
export clean_data
export estimate

export compute_objective_function

end
