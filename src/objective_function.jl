"""
    function compute_objective_function(log_price_uk::Vector{<:Number}, zOC::Vector{<:Number}, data::CleanData, params::Parameters)

Computes the objective function value based on the log prices and zOC values.

# Arguments
- `log_price_uk::Vector{<:Number}`: Logarithm of UK prices.
- `zOC::Vector{<:Number}`: Vector of zOC values.
- `data::CleanData`: Cleaned data structure containing industry and regional information.
- `params::Parameters`: Parameters structure containing production and constants.

# Returns
- `objective_value::Float64`: The computed objective function value.
"""
function compute_objective_function(log_price_uk::Vector{<:Number}, zOC::Vector{<:Number}, data::CleanData, params::ParameterSubset)

    tau = compute_advalorem_tax(data.industry)
    # tau = rand(length(log_price_uk))

    mu = params.production.shock_mean
    # mu = rand(length(log_price_uk))

    gammaK = params.production.capital
    # gammaK = rand(length(log_price_uk))

    k0 = data.industry.capital.current_year

    xi = params.constants.elasticities.production.substitution

    excess_demand = intermediate_goods_price_index(log_price_uk, zOC, tau, mu, gammaK, k0, xi)

    w = sqrt.(data.industry.regional.total_use.agg)
    e = w.*(excess_demand .- data.industry.regional.total_use.agg)
    objective_value = 0.5 * sum(e.^2)

    return objective_value

    # return sum(abs.(excess_demand))

end

"""
    compute_objective_function(x::Vector{<:Number}, data::CleanData, params::Parameters)

Computes the objective function value based on a vector of parameters.

# Arguments
- `x::Vector{<:Number}`: Vector containing log prices and zOC values.
- `data::CleanData`: Cleaned data structure containing industry and regional information.
- `params::Parameters`: Parameters structure containing production and constants.

# Returns
- `objective_value::Float64`: The computed objective function value.
"""
function compute_objective_function(x::Vector{<:Number}, data::CleanData, prices_eu::Vector{<:Number}, prices_world::Vector{<:Number})

    # n = length(data.industry.regional.total_use.agg)
    n = 16

    log_price_uk = x[1:n]
    zOC = x[(n+1):end]

    # params = compute_all_parameters(data, log_price_uk, prices_eu, prices_world, false)

    param_subset = compute_parameter_subset(data, log_price_uk, prices_eu, prices_world, false)

    println()
    println("Params ")
    println(param_subset.constants)
    println("===========")

    # exit()

    return compute_objective_function(log_price_uk, zOC, data, param_subset)

end


# function compute_objective_function(x::Vector{<:Number}, data::CleanData, prices_eu::Vector{<:Number}, prices_world::Vector{<:Number})

#     log_price_uk = x[1:n]

#     # params = compute_all_parameters(data, prices_df)

#     return compute_objective_function(x, data, params[1])

# end