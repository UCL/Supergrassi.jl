using Test, Supergrassi
using DataFrames

tol = 1e-12
threshold = 1e-4


function values_are_valid(arr::AbstractArray{<:T}, threshold::T) where {T <: Number}
    all(x -> x == 0 || x > threshold, arr)
end


@testset "Vectors" begin
    df = DataFrame([zeros(3), ones(3), rand(3) * threshold/10], ["uk", "eu", "world"])
    Supergrassi.round_shares!(df, threshold)

    @test values_are_valid(df.uk, threshold)
    @test values_are_valid(df.eu, threshold)
    @test values_are_valid(df.world, threshold)
    @test isapprox(df.uk + df.eu + df.world, ones(3), atol = tol)

end

@testset "Matrices" begin
    a = zeros(3,4)
    b = ones(3,4)
    c = rand(3,4) * threshold/10
    mat = Supergrassi.InputMatrices(a,b,c,a+b+c)
    Supergrassi.round_shares!(mat, threshold)

    @test values_are_valid(mat.uk, threshold)
    @test values_are_valid(mat.eu, threshold)
    @test values_are_valid(mat.world, threshold)
    @test isapprox(mat.uk + mat.eu + mat.world, ones(3,4), atol = tol)
end
