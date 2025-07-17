using JuMP
using Ipopt

function minimisation(x::Vector{<:Number}, data::CleanData, params::Parameters)
    model = Model(Ipopt.Optimizer)

    @variable(model, x >= 0)
    # @variable(model, y >= 0)

    @objective(model, Min, compute_objective_function(x, data, params))

    # @constraint(model, x + y == 5)

    optimize!(model)

    return value(x), value(y)

end

# Example usage
# x_opt, y_opt = minimisation()
# println("Optimal solution: x = $x_opt, y = $y_opt")
