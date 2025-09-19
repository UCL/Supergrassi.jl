using Enzyme
using Ipopt


function compute_gradient(x::Vector{<:Number}, gradient_var::Vector{Float64})

    gradient_var .= ones(length(x))
    return gradient_var

end
