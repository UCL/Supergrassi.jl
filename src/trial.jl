using Ipopt

function hess(x,
            rows,
            cols,
            obj_factor,
            lambda,
            values,
)
    return
end


function obj(x)
    y = (1 .- x) .^ 2
    return y[1]
end


function constraint(x, y)
    y = x
    return y
end

function grad_obj(x, y)
    y .= -2*(1 .- x)
    return y
end

function jac_constraint(x, rows, cols, vals)

    if vals === nothing
        rows[1] = 1
        cols[1] = 1
    else
        vals[1] = 1.0
    end
    
    return
end

prob = Ipopt.CreateIpoptProblem(
    1,
    [0.0],
    [10.0],
    1,
    [0.0],
    [10.0],
    1,
    0,
    obj,
    constraint,
    grad_obj,
    jac_constraint,
    hess,
)

Ipopt.IpoptSolve(prob)