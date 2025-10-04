using Base.Threads

function estimate(;log_results::Bool = false, log_results_filepath::String="log_results.csv")

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

    x = deepcopy([log_prices_uk; clean.industry.surplus.val; clean.industry.regional.totals.expenditure; household_final_consumption; clean.industry.depreciation.val])
    @info "Starting minimisation with x of length $(length(x))"

    global_jacobian = constraint_jacobian(x, log_prices_eu, log_prices_world, clean, params)
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

    @info "Final status: $status"
    @info "Final variables: $(prob.x)"
    @info "Final objective value: $(prob.obj_val)"

    if log_results

        if !isfile(log_results_filepath)
            open(log_results_filepath, "w") do io
                println(io, "timestamp,status,objective_value,variables")
            end
        end

        open(log_results_filepath, "a") do io
            println(io, "$(Dates.now()),$status,$(prob.obj_val),$(join(prob.x, ' '))")
        end

    end

    return status
end


function batch_estimation(;log_errors::Bool = false, log_errors_filepath::String="log_errors.csv", log_results::Bool = false, log_results_filepath::String="log_results.csv")

    @threads for i in 1:100
        try
            results = estimate(log_results=log_results, log_results_filepath=log_results_filepath)
        catch e

            if log_errors
                if !isfile(log_errors_filepath)
                    open(log_errors_filepath, "w") do io
                        println(io, "timestamp,iteration,error")
                    end
                end

                open(log_errors_filepath, "a") do io
                    println(io, "$(Dates.now()),$i,$(e.msg)")
                end
            end

        end
    end

end
