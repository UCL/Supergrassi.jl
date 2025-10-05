using Supergrassi, Test, Random

Random.seed!(1235)

@testset "End to End" begin
    results = Supergrassi.estimate()

    @test !isnothing(results)

end


@testset "Batch Estimation" begin

    Supergrassi.batch_estimation(batch_size=10, log_errors=true, log_errors_filepath="test_log_errors.csv", log_results=true, log_results_filepath="test_log_results.csv")

    @test isfile("test_log_results.csv")

    rm("test_log_results.csv")

end