"""
    function compute_constraint_function(x::Vector{<:Number}, data::CleanData, params::Parameters)

# Arguments
- `x::Vector{<:Number}`: Vector containing equilibrium variables
- `data::CleanData`: Cleaned data structure containing industry and regional information.
- `params::Parameters`: Parameters structure containing production and constants.
"""
function compute_constraint_function(x::Vector{<:Number}, data::CleanData, params::Parameters)

    return
    
end

"""
    function compute_constraint_function(log_price_uk::Vector{<:Number}, mu_I::Vector{<:Number}, log_Delta::Vector{<:Number}, data::CleanData, params::Parameters)

# Arguments
- `log_price_uk::Vector{<:Number}`: Logarithm of UK prices.
- `log_mu_I::Vector{<:Number}`: logarithm of investment TFP.
- `log_Delta::Vector{<:Number}`: logarithm of depreciation parameter
- `data::CleanData`: Cleaned data structure containing industry and regional information.
- `params::Parameters`: Parameters structure containing production and constants.

"""
function compute_constraint_function(log_price_uk::Vector{<:Number}, mu_I::Vector{<:Number}, log_Delta::Vector{<:Number}, data::CleanData, params::Parameters)

    return
    
end
