using Enzyme
using Ipopt


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

    mu = params.production.shock_mean

    gammaK = params.production.capital

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
function compute_objective_function(x::Vector{<:Number}, data::CleanData, prices_eu::Vector{<:Number}, prices_world::Vector{<:Number})

    n = data.constants.number_of_industries

    log_price_uk = x[1:n]
    zOC = x[(n+1):2*n]

    param_subset = compute_parameter_subset(data, log_price_uk, prices_eu, prices_world, false)

    return compute_objective_function(log_price_uk, zOC, data, param_subset)


end

"""
    function compute_gradient(x::Vector{<:Number}, clean_data::CleanData, prices_eu::Vector{<:Number}, prices_world::Vector{<:Number}, gradient_var::Vector{Float64})
        
Computes the gradient of the objective function with respect to the input vector `x`.

# Arguments
- `x::Vector{<:Number}`: Vector containing log prices and zOC values.
- `clean_data::CleanData`: Cleaned data structure containing industry and regional information.
- `prices_eu::Vector{<:Number}`: Vector containing the EU prices.
- `prices_world::Vector{<:Number}`: Vector containing the rest of the world prices.
- `gradient_var::Vector{Float64}`: Vector to store the computed gradient values.

# Returns
- `gradient_var::Vector{Float64}`: The computed gradient vector.
"""
function compute_gradient(x::Vector{<:Number}, clean_data::CleanData, prices_eu::Vector{<:Number}, prices_world::Vector{<:Number}, gradient_var::Vector{Float64})
    
    gradient_var .= gradient(set_runtime_activity(Reverse), compute_objective_function, x, Const(clean_data), Const(prices_eu), Const(prices_world))[1]
    
    @debug("Gradient computed: ", grad[1])
    @debug("Length of gradient: ", length(grad[1]))
    @debug("Length of gradient_var: ", length(gradient_var))

    return gradient_var



end
