using Test
using Supergrassi
using DataFrames

merge_codes = [1, 1, 2, 35, 35, 72, 4, 4, 4]
assets = ["1", "0.5", "NA", "1.5", "3.7", "NA", "10748", "999", "1.0"]
liabilities = ["NA", "17.5", "1199432", "13.85", "0.0", "9", "1142", "95", "1.0"]

mapping = Dict{String, String}("SIC_64_001" => "A",
                               "SIC_64_002" => "A",
                               "SIC_64_004" => "C",
                               "SIC_64_035" => "B",
                               "SIC_64_072" => "B")

df = DataFrame(sic64 = merge_codes, total_assets_2025 = assets, total_liabilities_2025 = liabilities)

@testset "clean assets liabilities" begin

    clean = Supergrassi.clean_assets_liabilities(df, 2025, mapping)
    @test nrow(clean) == 4
    @test all(clean.Ratio .<= 1)

    sampled = Supergrassi.clean_assets_liabilities(df, 2025, mapping, 1)
    @test nrow(sampled) == 2
    @test all(sampled.Ratio .<= 1)
    
end
