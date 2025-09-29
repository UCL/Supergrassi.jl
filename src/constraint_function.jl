using LinearAlgebra

"""
    function compute_constraint_function(x::Vector{<:Number}, data::CleanData, params::Parameters)

# Arguments
- `x::Vector{<:Number}`: Vector containing equilibrium variables
- `data::CleanData`: Cleaned data structure containing industry and regional information.
- `params::Parameters`: Parameters structure containing production and constants.
"""
function constraint_jacobian(x::Vector{<:Number}, log_price_eu::Vector{T}, log_price_world::Vector{T},
                             data::CleanData, params::Parameters) where {T <: Real}

    n = data.constants.number_of_industries

    log_price_uk, zOC, expenditure, log_TFP, log_Delta = unpack_x(n, x)

    params = compute_all_parameters(data, log_price_uk, log_price_eu, log_price_world)

    Jac = jacobian(set_runtime_activity(Forward),
                   constraint_function,
                   x,
                   Const(log_price_eu),
                   Const(log_price_world),
                   Const(params),
                   Const(data.industry),
                   Const(data.constants),
                   Const(zeros(2*n))
                   )

    DCEQ = first(Jac)

    return DCEQ

end


function constraint_function(x::Vector{T}, price_eu::Vector{T}, price_world::Vector{T},
                             params::Parameters, data::IndustryData, constants::Constants, y::Vector{T}) where {T <: Real}

    F = market_clearing_price_constraint(x, price_eu, price_world, params, data, constants)
    CFC = compute_fixed_capital_consumption_constraint(x, data, params)
    y = [F; CFC]
    return y

end
