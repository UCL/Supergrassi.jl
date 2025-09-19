module Supergrassi

using Enzyme
using Ipopt
using LinearAlgebra
using SparseArrays


function sparser(matrix)
    sp = sparse(matrix)
    return findnz(sp)
end

function compute_objective_function(x::Vector{<:Number})
    return sum(x)
end

function compute_gradient(x::Vector{<:Number}, gradient_var::Vector{Float64})
    gradient_var .= ones(length(x))
    return gradient_var
end

function compute_constraint_jacobian(x::Vector{<:Number})
    y = diagm(ones(length(x[1:32])))
    return y
end

function constraint_function(x::Vector{T}, y::Vector{T}) where {T <: Real}
    y = 0.5 .* x[1:32] .^ 2
    return y
end



function compute_hessian(x::Vector{Float64},
                  rows::Vector{Int32},
                  cols::Vector{Int32},
                  obj_factor::Float64,
                  lambda::Vector{Float64},
                  values::Union{Nothing,Vector{Float64}},
                  )

    return
end




function estimate()

    x = rand(50)
    println("Starting minimisation with x of length ", length(x))


    function simple_jacobian(x, rows, cols, vals)

        sparse_jacobian = sparser(compute_constraint_jacobian(x))
        
        if vals === nothing
            rows[:] = sparse_jacobian[1]
            cols[:] = sparse_jacobian[2]
        else
            vals[:] = sparse_jacobian[3]
        end

        return 
    end

    function compute_non_zero_elements_in_jacobian(x)
        trial = sparser(compute_constraint_jacobian(x))
        return length(trial[1])
    end

    n_non_zero = compute_non_zero_elements_in_jacobian(x)

    prob = Ipopt.CreateIpoptProblem(
        50,
        [0.0 for i in 1:50],
        [10.0 for i in 1:50],
        32,
        [0.0 for i in 1:32],
        [10.0 for i in 1:32],
        n_non_zero,
        0,
        compute_objective_function,
        constraint_function,
        compute_gradient,
        simple_jacobian,
        compute_hessian,
    )

    prob.x = x

    status = Ipopt.IpoptSolve(prob)

    return status
end

end
