using Enzyme
using Ipopt


function compute_gradient(x::Vector{<:Number}, clean_data::CleanData, prices_eu::Vector{<:Number}, prices_world::Vector{<:Number}, gradient_var::Vector{Float64})
    
    # gradient_var .= zeros(length(x))
    # grad = gradient(set_runtime_activity(Reverse), compute_objective_function, x, Const(clean_data), Const(prices_eu), Const(prices_world))
    
    # println("Gradient computed: ", grad[1])
    # println("Length of gradient: ", length(grad[1]))
    # println("Length of gradient_var: ", length(gradient_var))
    
    # gradient_var .= grad[1]

    # return gradient_var

    gradient_var .= ones(length(x))
    return gradient_var

end
