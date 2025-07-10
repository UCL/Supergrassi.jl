using DataFrames
using CSV
using Supergrassi
using Test

m = Supergrassi.InputMatrices(Matrix(reshape(df2d_clean_data.mValueUK, (16, 16))),
                              Matrix(reshape(df2d_clean_data.mValueEU, (16, 16))),
                              Matrix(reshape(df2d_clean_data.mValueW, (16, 16))),
                              Matrix(reshape(df2d_clean_data.mValue, (16, 16))))

@testset "Clean 1d Data" begin

    @test isapprox(clean.industry.regional.delta_v.agg, df_clean_data.DeltaVValue)
    @test isapprox(clean.industry.regional.delta_v.uk, df_clean_data.DeltaVValueUK)
    @test isapprox(clean.industry.regional.delta_v.eu, df_clean_data.DeltaVValueEU)
    @test isapprox(clean.industry.regional.delta_v.world, df_clean_data.DeltaVValueW)

    @test isapprox(clean.industry.regional.investment.agg,  df_clean_data.IValue)
    @test isapprox(clean.industry.regional.investment.uk,  df_clean_data.IValueUK)
    @test isapprox(clean.industry.regional.investment.eu,  df_clean_data.IValueEU)
    @test isapprox(clean.industry.regional.investment.world,  df_clean_data.IValueW)

    @test isapprox(clean.industry.regional.consumption.agg,  df_clean_data.fValue)
    @test isapprox(clean.industry.regional.consumption.uk,  df_clean_data.fValueUK)
    @test isapprox(clean.industry.regional.consumption.eu,  df_clean_data.fValueEU)
    @test isapprox(clean.industry.regional.consumption.world,  df_clean_data.fValueW)

    @test isapprox(clean.industry.regional.export_eu.agg,  df_clean_data.x1Value)
    @test isapprox(clean.industry.regional.export_eu.uk,  df_clean_data.x1ValueUK)
    @test isapprox(clean.industry.regional.export_eu.eu,  df_clean_data.x1ValueEU)
    @test isapprox(clean.industry.regional.export_eu.world,  df_clean_data.x1ValueW)

    @test isapprox(clean.industry.regional.export_world.agg,  df_clean_data.x2Value)
    @test isapprox(clean.industry.regional.export_world.uk,  df_clean_data.x2ValueUK)
    @test isapprox(clean.industry.regional.export_world.eu,  df_clean_data.x2ValueEU)
    @test isapprox(clean.industry.regional.export_world.world,  df_clean_data.x2ValueW)

    @test isapprox(clean.industry.regional.total_use.agg,  df_clean_data.yValue)
    @test isapprox(clean.industry.regional.total_use.uk,  df_clean_data.yValueUK)
    @test isapprox(clean.industry.regional.total_use.eu,  df_clean_data.yValueEU)
    @test isapprox(clean.industry.regional.total_use.world,  df_clean_data.yValueW)

    @test isapprox(clean.household.income.low, df_clean_data.income_lo)
    @test isapprox(clean.household.income.high, df_clean_data.income_hi)

    @test isapprox(clean.household.hours.low, df_clean_data.hLO, nans=true)
    @test isapprox(clean.household.hours.high, df_clean_data.hHI, nans=true)

    @test isapprox(clean.household.wages.low, df_clean_data.wLO, nans=true)
    @test isapprox(clean.household.wages.high, df_clean_data.wHI, nans=true)
    
    @test isapprox(clean.household.payments.low, df_clean_data.hValueLO, nans=true)
    @test isapprox(clean.household.payments.high, df_clean_data.hValueHI, nans=true)
    @test isapprox(clean.household.payments.agg, df_clean_data.hValue, nans=true)

    @test isapprox(clean.industry.capital.current_year, df_clean_data.k0)
    @test isapprox(clean.industry.capital.next_year, df_clean_data.k1)

    @test isapprox(clean.industry.tax.products, df_clean_data.taxValue1)
    @test isapprox(clean.industry.tax.production, df_clean_data.taxValue2)

    @test isapprox(clean.industry.surplus.val, df_clean_data.kValue)
    @test isapprox(clean.industry.depreciation.val, df_clean_data.depreciation)
    @test isapprox(clean.industry.shock_stdev.val, df_clean_data.sigmaBar)

end

@testset "Clean 2d data" begin

    @test isapprox(clean.industry.regional.input_matrices.uk, m.uk, atol = tol)
    @test isapprox(clean.industry.regional.input_matrices.eu, m.eu, atol = tol)
    @test isapprox(clean.industry.regional.input_matrices.world, m.world, atol = tol)
    @test isapprox(clean.industry.regional.input_matrices.agg, m.agg, atol = tol)
    
end
