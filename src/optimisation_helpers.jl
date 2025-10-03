using SparseArrays


function compute_hessian(x::Vector{Float64},
                  rows::Vector{Int32},
                  cols::Vector{Int32},
                  obj_factor::Float64,
                  lambda::Vector{Float64},
                  values::Union{Nothing,Vector{Float64}},
                  )

    return
end


function sparser(matrix)
    sp = sparse(matrix)
    return findnz(sp)
end


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
