using JuMP
using Ipopt
using Enzyme


function compute_gradient(x::Vector{<:Number}, clean_data::CleanData, prices_eu::Vector{<:Number}, prices_world::Vector{<:Number})

    grad = gradient(set_runtime_activity(ForwardWithPrimal), compute_objective_function, x, Const(clean_data), Const(prices_eu), Const(prices_world))

    return grad

end

# clean_data = Supergrassi.CleanData()
# params = Supergrassi.Parameters()

# grad = gradient(ForwardWithPrimal, create_objective_function, [1.0, 2.0, 3.0], Supergrassi.CleanData(), Supergrassi.Parameters())

# function dummy_compute_constraint_function(x::Vector{<:Number}, data::CleanData, params::Parameters)
#     # Dummy constraint function for demonstration purposes
#     # Replace with actual constraint logic
#     return true
# end

# function minimisation(x::Vector{<:Number}, data::CleanData, prices_df::DataFrame)
#     model = Model(Ipopt.Optimizer)

#     n = 16

#     xx = @variable(model, y)

#     println(xx.value)
#     println(typeof(xx))
#     return

#     @objective(model, Min, op_compute_objective_function(xx))
#     # @objective(model, Min, compute_objective_function(x, data, prices_df))

#     @constraint(model, [1 > 0, 2 > 0] .== true)

#     optimize!(model)

#     prices_df.uk = x[1:n]

#     return value(x)

# end



# function minimisation(x::Vector{<:Number}, data::CleanData, prices_df::DataFrame)
#     model = Model(Ipopt.Optimizer)

#     n = 16

#     @variable(model, y)
#     @operator(model, op_compute_objective_function, 200, compute_objective_function)

#     xx = serialise(x, data, prices_df)

#     @objective(model, Min, op_compute_objective_function(xx))
#     # @objective(model, Min, compute_objective_function(x, data, prices_df))

#     @constraint(model, [1 > 0, 2 > 0] .== true)

#     optimize!(model)

#     prices_df.uk = x[1:n]

#     return value(x)

# end