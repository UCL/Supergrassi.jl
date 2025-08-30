using Supergrassi, Test


@testset "End to End" begin
    settings, data, clean, params, log_params = Supergrassi.estimate()

    @test !isnothing(settings)
    @test !isnothing(data)
    @test !isnothing(clean)
    @test !isnothing(params)
    @test !isnothing(log_params)
end
