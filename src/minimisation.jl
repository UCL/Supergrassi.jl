using Enzyme


function compute_gradient(x::Vector{<:Number}, clean_data::CleanData, prices_eu::Vector{<:Number}, prices_world::Vector{<:Number})

    grad = gradient(set_runtime_activity(ForwardWithPrimal), compute_objective_function, x, Const(clean_data), Const(prices_eu), Const(prices_world))

    return grad

end
