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
function compute_objective_function(log_price_uk::Vector{<:Number}, zOC::Vector{<:Number}, data::CleanData, params::Parameters)

    tau = compute_advalorem_tax(data.industry)

    mu = params.production.shock_mean
    gammaK = params.production.input_capital
    k0 = data.industry.capital.current_year
    xi = params.constants.elasticities.production.substitution

    excess_demand = intermediate_goods_price_index(log_price_uk, zOC, tau, mu, gammaK, k0, xi)

    w = sqrt.(data.industry.regional.total_use.agg)
    e = w.*(excess_demand .- data.industry.regional.total_use.agg)
    objective_value = 0.5 * sum(e.^2)

    return objective_value

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
function compute_objective_function(x::Vector{<:Number}, data::CleanData, params::Parameters)

    n = length(data.industry.regional.total_use.agg)

    log_price_uk = x[1:n]
    zOC = x[(n+1):end]

    return compute_objective_function(log_price_uk, zOC, data, params)

end
