using JuMP
using Ipopt

function minimisation()
    model = Model(Ipopt.Optimizer)

    @variable(model, x >= 0)
    @variable(model, y >= 0)

    @objective(model, Min, (x - 1)^2 + (y - 2)^2)

    @constraint(model, x + y == 5)

    optimize!(model)

    return value(x), value(y)

end

# Example usage
x_opt, y_opt = minimisation()
println("Optimal solution: x = $x_opt, y = $y_opt")
