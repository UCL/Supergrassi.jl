using LinearAlgebra

"""
    function compute_constraint_function(x::Vector{<:Number}, data::CleanData, params::Parameters)

# Arguments
- `x::Vector{<:Number}`: Vector containing equilibrium variables
- `data::CleanData`: Cleaned data structure containing industry and regional information.
- `params::Parameters`: Parameters structure containing production and constants.
"""
function compute_constraint_jacobian(x::Vector{<:Number})

    y = diagm(ones(length(x[1:32])))
    return y

end


function constraint_wrapper(x::Vector{T}, y::Vector{T}) where {T <: Real}

    y = 0.5 .* x[1:32] .^ 2
    return y

end
