module Supergrassi

using Enzyme
using Ipopt
using LinearAlgebra
using SparseArrays

function compute_objective_function(x::Vector{<:Number})
    return sum(x)
end

function compute_gradient(x::Vector{<:Number}, y::Vector{Float64})
    y[:] = ones(length(x))
    return y
end

# function compute_constraint_jacobian(x::Vector{<:Number})
#     #y = diagm(ones(length(x[1:4])))
#     y = diagm(ones(4))
#     return y
# end

function constraint_function(x::Vector{T}, y::Vector{T}) where {T <: Real}
    y = 0.5 .* x[1:4] .^ 2
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

function simple_jacobian(x::Vector{Float64}, rows::Vector{Int32}, cols::Vector{Int32}, vals::Nothing)

    if ~isdefined(Supergrassi, :Jg)
        error("Jacobian must be defined in the scope that calls this function")
    end
    
    sparse_jacobian = findnz(sparse(Jg))
    rows[:] = sparse_jacobian[1]
    cols[:] = sparse_jacobian[2]

    return
    
end

function simple_jacobian(x, rows, cols, vals::Vector{Float64})
    
    sparse_jacobian = findnz(sparse(diagm(ones(length(x[1:4])))))
    vals[:] = sparse_jacobian[3]
    
    return 
end

Jg = diagm(1 => ones(3), 3 => [1])

function estimate()

    x = rand(5)
    #Jg = diagm(ones(length(x[1:4])))
    println("Starting minimisation with x of length ", length(x))
    
    n_non_zero = 4

    prob = Ipopt.CreateIpoptProblem(
        5,
        [0.0 for i in 1:5],
        [10.0 for i in 1:5],
        4,
        [0.0 for i in 1:4],
        [10.0 for i in 1:4],
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
