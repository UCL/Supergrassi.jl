using DataFrames
using CSV
using Supergrassi
using Test

# This test assumes we have a data set available in input/uk_data

path = joinpath(@__DIR__, "..", "config","settings.yml")
settings_path = create_filepath(path)
settings = read_settings(settings_path)
filepaths = check_file_availability(settings)
data = read_data(filepaths, settings)

clean = Supergrassi.clean_data(data,settings)
Supergrassi.postprocess_clean_data!(clean)

df = CSV.read(joinpath(@__DIR__, "..", "data", "test_load_data.csv"), DataFrame)
df2d = CSV.read(joinpath(@__DIR__, "..", "data", "test_load_data_2d.csv"), DataFrame)

nms = names(clean.industry.regional.input_matrices.uk)

m = Supergrassi.InputMatrices(DataFrame(reshape(df2d.mValueUK, (16, 16)), nms),
                              DataFrame(reshape(df2d.mValueEU, (16, 16)), nms),
                              DataFrame(reshape(df2d.mValueW, (16, 16)), nms),
                              DataFrame(reshape(df2d.mValueIMP, (16, 16)), nms),
                              DataFrame(reshape(df2d.mValue, (16, 16)), nms))

@testset "Clean 1d Data" begin

    @test isapprox(clean.industry.regional.delta_v.agg, df.DeltaVValue)
    @test isapprox(clean.industry.regional.delta_v.uk, df.DeltaVValueUK)
    @test isapprox(clean.industry.regional.delta_v.eu, df.DeltaVValueEU)
    @test isapprox(clean.industry.regional.delta_v.world, df.DeltaVValueW)

    @test isapprox(clean.industry.regional.investment.agg,  df.IValue)
    @test isapprox(clean.industry.regional.investment.uk,  df.IValueUK)
    @test isapprox(clean.industry.regional.investment.eu,  df.IValueEU)
    @test isapprox(clean.industry.regional.investment.world,  df.IValueW)

    @test isapprox(clean.industry.regional.consumption.agg,  df.fValue)
    @test isapprox(clean.industry.regional.consumption.uk,  df.fValueUK)
    @test isapprox(clean.industry.regional.consumption.eu,  df.fValueEU)
    @test isapprox(clean.industry.regional.consumption.world,  df.fValueW)

    @test isapprox(clean.industry.regional.export_eu.agg,  df.x1Value)
    @test isapprox(clean.industry.regional.export_eu.uk,  df.x1ValueUK)
    @test isapprox(clean.industry.regional.export_eu.eu,  df.x1ValueEU)
    @test isapprox(clean.industry.regional.export_eu.world,  df.x1ValueW)

    @test isapprox(clean.industry.regional.export_world.agg,  df.x2Value)
    @test isapprox(clean.industry.regional.export_world.uk,  df.x2ValueUK)
    @test isapprox(clean.industry.regional.export_world.eu,  df.x2ValueEU)
    @test isapprox(clean.industry.regional.export_world.world,  df.x2ValueW)

    @test isapprox(clean.industry.regional.total_use.agg,  df.yValue)
    @test isapprox(clean.industry.regional.total_use.uk,  df.yValueUK)
    @test isapprox(clean.industry.regional.total_use.eu,  df.yValueEU)
    @test isapprox(clean.industry.regional.total_use.world,  df.yValueW)

    @test isapprox(clean.household.income.low, df.income_lo)
    @test isapprox(clean.household.income.high, df.income_hi)

    @test isapprox(clean.household.hours.low, df.hLO, nans=true)
    @test isapprox(clean.household.hours.high, df.hHI, nans=true)

    @test isapprox(clean.household.wages.low, df.wLO, nans=true)
    @test isapprox(clean.household.wages.high, df.wHI, nans=true)
    
    @test isapprox(clean.household.payments.low, df.hValueLO, nans=true)
    @test isapprox(clean.household.payments.high, df.hValueHI, nans=true)
    @test isapprox(clean.household.payments.agg, df.hValue, nans=true)

    @test isapprox(clean.industry.capital.current_year, df.k0)
    @test isapprox(clean.industry.capital.next_year, df.k1)

    @test isapprox(clean.industry.tax.products, df.taxValue1)
    @test isapprox(clean.industry.tax.production, df.taxValue2)

    @test isapprox(clean.household.income_share.high, df.hi_income_share, nans=true)
    @test isapprox(clean.industry.surplus.val, df.kValue)
    @test isapprox(clean.industry.depreciation.val, df.depreciation)
    @test isapprox(clean.industry.shock_stdev.val, df.sigmaBar)

end

@testset "Clean 2d data" begin

    @test isapprox(clean.industry.regional.input_matrices.uk, m.uk)
    @test isapprox(clean.industry.regional.input_matrices.eu, m.eu)
    @test isapprox(clean.industry.regional.input_matrices.world, m.world)
    @test isapprox(clean.industry.regional.input_matrices.agg, m.agg)
    
end
