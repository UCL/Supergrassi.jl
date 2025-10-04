using SparseArrays

"""
    function compute_gradient(x::Vector{<:Number}, data::CleanData, prices_eu::Vector{<:Number}, prices_world::Vector{<:Number}, gradient_var::Vector{<:Number})

Dummy function to compute the gradient of the objective function with respect to the input vector `x` as required by Ipopt.

# Arguments
- `x::Vector{<:Number}`: Vector containing log prices and zOC values
- `data::CleanData`: Cleaned data structure containing industry and regional information.
- `prices_eu::Vector{<:Number}`: Vector containing the EU prices.
- `prices_world::Vector{<:Number}`: Vector containing the rest of the world
- `gradient_var::Vector{<:Number}`: Vector to store the computed gradient values.

# Returns
- `gradient_var::Vector{<:Number}`: The computed gradient vector.
"""

function compute_hessian(x::Vector{Float64},
                  rows::Vector{Int32},
                  cols::Vector{Int32},
                  obj_factor::Float64,
                  lambda::Vector{Float64},
                  values::Union{Nothing,Vector{Float64}},
                  )

    return
end

"""
    function sparser(matrix::AbstractMatrix)

Converts a dense matrix to a sparse representation suitable for Ipopt.

# Arguments
- `matrix::AbstractMatrix`: The input dense matrix.

# Returns
- `Tuple{Vector{Int}, Vector{Int}, Vector{Float64}}`: A tuple containing the row indices, column indices, and non-zero values of the sparse matrix.
"""

function sparser(matrix)
    sp = sparse(matrix)
    return findnz(sp)
end

"""
    function create_optimization_functions(log_prices_eu, log_prices_world, clean, params, global_jacobian)

Creates the objective, gradient, constraint, and jacobian functions for optimization using Ipopt.

# Arguments
- `log_prices_eu`: Vector containing the EU prices.
- `log_prices_world`: Vector containing the rest of the world prices.
- `clean`: Cleaned data structure containing industry and regional information.
- `params`: Parameters structure containing production and constants.
- `global_jacobian`: Precomputed global jacobian matrix.

# Returns
- `simple_objective`: Function to compute the objective value.
- `simple_gradient`: Function to compute the gradient of the objective.
- `simple_constraint`: Function to compute the constraints.
- `simple_jacobian`: Function to compute the jacobian of the constraints.
"""
function create_optimization_functions(log_prices_eu, log_prices_world, clean, params, global_jacobian)

    function simple_jacobian(x::Vector{Float64}, rows::Vector{Int32}, cols::Vector{Int32}, vals::Vector{Float64})
        sparse_jacobian = sparser(constraint_jacobian(x, log_prices_eu, log_prices_world, clean, params))
        vals[:] = sparse_jacobian[3]
        return 
    end

    function simple_jacobian(x::Vector{Float64}, rows::Vector{Int32}, cols::Vector{Int32}, vals::Nothing)
        sparse_jacobian = findnz(sparse(global_jacobian))
        rows[:] = sparse_jacobian[1]
        cols[:] = sparse_jacobian[2]
        return 
    end

    simple_objective(x) = compute_objective_function(x, clean, log_prices_eu, log_prices_world)
    simple_gradient(x, gradient_var) = compute_gradient(x, clean, log_prices_eu, log_prices_world, gradient_var)
    simple_constraint(x, constraint) = constraint_function(x, log_prices_eu, log_prices_world, clean, params, constraint)
    
    return simple_objective, simple_gradient, simple_constraint, simple_jacobian
end
