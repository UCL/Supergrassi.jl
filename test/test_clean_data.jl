using CSV
using Supergrassi
using Test

# This test assumes we have a data set available in input/uk_data

@testset "Clean Data" begin

    settings_path = create_filepath("../config/settings.yml")
    settings = read_settings(settings_path)
    filepaths = check_file_availability(settings)
    data = read_data(filepaths, settings)

    clean = Supergrassi.clean_data(data,settings["constants"]["data_year"])

    df = CSV.read(joinpath(@__DIR__, "..", "data", "test_load_data.csv"), DataFrame)

    @test isapprox(clean.delta_v.uk, df.DeltaVValueUK)
    @test isapprox(clean.delta_v.eu, df.DeltaVValueEU)
    @test isapprox(clean.delta_v.world, df.DeltaVValueW)
    @test isapprox(clean.delta_v.imports, df.DeltaVValueIMP)

    @test isapprox(clean.capital_formation.uk,  df.IValueUK)
    @test isapprox(clean.capital_formation.eu,  df.IValueEU)
    @test isapprox(clean.capital_formation.world,  df.IValueW)
    @test isapprox(clean.capital_formation.imports,  df.IValueIMP)

    @test isapprox(clean.consumption.uk,  df.fValueUK)
    @test isapprox(clean.consumption.eu,  df.fValueEU)
    @test isapprox(clean.consumption.world,  df.fValueW)
    @test isapprox(clean.consumption.imports,  df.fValueIMP)

    @test isapprox(clean.export_eu.uk,  df.x1ValueUK)
    @test isapprox(clean.export_eu.eu,  df.x1ValueEU)
    @test isapprox(clean.export_eu.world,  df.x1ValueW)
    @test isapprox(clean.export_eu.imports,  df.x1ValueIMP)

    @test isapprox(clean.export_world.uk,  df.x2ValueUK)
    @test isapprox(clean.export_world.eu,  df.x2ValueEU)
    @test isapprox(clean.export_world.world,  df.x2ValueW)
    @test isapprox(clean.export_world.imports,  df.x2ValueIMP)

    @test isapprox(clean.total_use.uk,  df.yValueUK)
    @test isapprox(clean.total_use.eu,  df.yValueEU)
    @test isapprox(clean.total_use.world,  df.yValueW)
    @test isapprox(clean.total_use.imports,  df.yValueIMP)

    @test isapprox(clean.income.low, df.income_lo)
    @test isapprox(clean.income.high, df.income_hi)

    @test isapprox(clean.payments.low, df.hValueLO, nans=true)
    @test isapprox(clean.payments.high, df.hValueHI, nans=true)
    @test isapprox(clean.payments.agg, df.hValue, nans=true)

    @test isapprox(clean.mean_capital.current, df.k0)
    @test isapprox(clean.mean_capital.next, df.k1)

    @test isapprox(clean.tax.products, df.taxValue1)
    @test isapprox(clean.tax.production, df.taxValue2)

    @test isapprox(clean.income_share.high, df.hi_income_share, nans=true)
    @test isapprox(clean.operating_surplus.val, df.kValue)
    @test isapprox(clean.depreciation.val, df.depreciation)
    @test isapprox(clean.shock_stddev.val, df.sigmaBar)

end
