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

    @test isapprox(clean.industry.regional.delta_v.uk, df.DeltaVValueUK)
    @test isapprox(clean.industry.regional.delta_v.eu, df.DeltaVValueEU)
    @test isapprox(clean.industry.regional.delta_v.world, df.DeltaVValueW)
    @test isapprox(clean.industry.regional.delta_v.imports, df.DeltaVValueIMP)

    @test isapprox(clean.industry.regional.capital_formation.uk,  df.IValueUK)
    @test isapprox(clean.industry.regional.capital_formation.eu,  df.IValueEU)
    @test isapprox(clean.industry.regional.capital_formation.world,  df.IValueW)
    @test isapprox(clean.industry.regional.capital_formation.imports,  df.IValueIMP)

    @test isapprox(clean.industry.regional.consumption.uk,  df.fValueUK)
    @test isapprox(clean.industry.regional.consumption.eu,  df.fValueEU)
    @test isapprox(clean.industry.regional.consumption.world,  df.fValueW)
    @test isapprox(clean.industry.regional.consumption.imports,  df.fValueIMP)

    @test isapprox(clean.industry.regional.export_eu.uk,  df.x1ValueUK)
    @test isapprox(clean.industry.regional.export_eu.eu,  df.x1ValueEU)
    @test isapprox(clean.industry.regional.export_eu.world,  df.x1ValueW)
    @test isapprox(clean.industry.regional.export_eu.imports,  df.x1ValueIMP)

    @test isapprox(clean.industry.regional.export_world.uk,  df.x2ValueUK)
    @test isapprox(clean.industry.regional.export_world.eu,  df.x2ValueEU)
    @test isapprox(clean.industry.regional.export_world.world,  df.x2ValueW)
    @test isapprox(clean.industry.regional.export_world.imports,  df.x2ValueIMP)

    @test isapprox(clean.industry.regional.total_use.uk,  df.yValueUK)
    @test isapprox(clean.industry.regional.total_use.eu,  df.yValueEU)
    @test isapprox(clean.industry.regional.total_use.world,  df.yValueW)
    @test isapprox(clean.industry.regional.total_use.imports,  df.yValueIMP)

    @test isapprox(clean.household.income.low, df.income_lo)
    @test isapprox(clean.household.income.high, df.income_hi)

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
