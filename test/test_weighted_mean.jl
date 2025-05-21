using Test
using Supergrassi
using DataFrames

merge_codes = DataFrame(x1 = ["foo","1","2","3","4","5","6"],
                        x7 = ["bar","A","A","B","C","C","C"])

map_64 = Supergrassi.create_map_64_to_16(merge_codes)

df = DataFrame(randn(10,6),
               ["SIC_64_001","SIC_64_002","SIC_64_003","SIC_64_004","SIC_64_005","SIC_64_006"])

weights = DataFrame(ones(1,6),
                    names(df))

@testset "weighted mean" begin

    df2 = Supergrassi.reduce_columns_by_group_weighted_mean(df, map_64, weights = weights)
    @test nrow(df2) == 10
    @test ncol(df2) == 3
    @test isapprox((df.SIC_64_001 + df.SIC_64_002) / 2, df2.A)
    @test isapprox(df.SIC_64_003, df2.B)
    @test isapprox((df.SIC_64_004 + df.SIC_64_005 + df.SIC_64_006) / 3, df2.C)
    
end
