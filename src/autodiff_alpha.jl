using ADTypes
using BenchmarkTools
using DifferentiationInterface
import ForwardDiff, Zygote
using Random
using SparseConnectivityTracer
using SparseMatrixColorings
using Supergrassi

bakcend = AutoForwardDiff()
pd = 

gradient(Supergrassi.compute_demand,
         backend,
         Constant(n),
         Constant(elasticity),
         pd,
         Constant(log_price.eu),
         Constant(log_price.world),
         Constant(goods_consumption))
