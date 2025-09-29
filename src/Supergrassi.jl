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
include("optimisation_helpers.jl")


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
    @info "Starting minimisation with x of length $(length(x))"

    global_jacobian = compute_constraint_function(x, log_prices_eu, log_prices_world, clean, params)
    @info "Global Jacobian computed."
    simple_objective, simple_gradient, simple_constraint, simple_jacobian = 
        create_optimization_functions(log_prices_eu, log_prices_world, clean, params, global_jacobian)
    @info "Optimization functions created."

    trial = sparser(global_jacobian)

    prob = Ipopt.CreateIpoptProblem(
        50,
        [0.0 for i in 1:50],
        [10.0 for i in 1:50],
        32,
        [0.0 for i in 1:32],
        [0.0 for i in 1:32],
        length(trial[1]),
        0,
        simple_objective,
        simple_constraint,
        simple_gradient,
        simple_jacobian,
        compute_hessian,
    )

    prob.x = x

    status = Ipopt.IpoptSolve(prob)

    @info "Estimation completed."

    return status
end

export create_filepath, read_data, read_settings, check_file_availability
export clean_data
export estimate

export compute_objective_function

end
