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

