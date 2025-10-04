"""
    convert_to_ratio!(data::RegionalData)

Convert regional data into ratios of region / sum(regions)

# Arguments
- `data::RegionalData`: The RegionalData struct containing regional data to be converted into ratios.
"""
function convert_to_ratio!(data::RegionalData)

    for field in fieldnames(RegionalData)

        df = getfield(data, field)

        if (field != :delta_v && field != :total_use  && field != :totals)

            convert_to_ratio!(df)

        end

    end

end

"""
    convert_to_ratio!(data::HouseholdData)

Convert household data into ratios of low and high payments to the sum of low and high payments.

# Arguments
- `data::HouseholdData`: The HouseholdData struct containing household data to be converted into ratios.
"""
function convert_to_ratio!(data::HouseholdData)

    for field in [:payments]

        df = getfield(data, field)
        scaling_factor = 1.0 ./ (df.low .+ df.high)

        for col in [:low, :high]

            df[!,col] .*= scaling_factor
            replace!(df[!, col], NaN => 0.0)

        end

    end

end

"""
    convert_to_ratio!(df::DataFrame)

Convert regional vector data into fractions of the aggregate value.
Then renormalise the aggregate value by its sum.
Nans are replaced with 0's.

# Arguments
- `df::DataFrame`: The DataFrame containing regional vector data to be converted into ratios.
"""
function convert_to_ratio!(df::DataFrame)

    for col in [:uk, :eu, :world]
        df[!, col] ./= df.agg
        replace!(df[!, col], NaN => 0.0)
    end

    df.agg ./= sum(df.agg)

end

"""
    convert_to_ratio!(data::InputMatrices)

Convert regional matrix data to fractions of the aggregate value
Nans are replaced with 0's

# Arguments
- `data::InputMatrices`: The InputMatrices struct containing regional matrix data to be converted into ratios.
"""
function convert_to_ratio!(data::InputMatrices)

    for field in [:uk, :eu, :world]

        df = getfield(data, field)
        df ./= data.agg

        for c in eachcol(df)
            replace!(c, NaN => 0.0)
        end

    end

end

"""
    round_shares!(data::RegionalData, threshold::Float64 = 1e-4)

Round values below threshold in regional data to 0, then rescale so that regions sum to 1.

# Arguments
- `data::RegionalData`: The RegionalData struct containing regional data to be rounded.
"""
function round_shares!(data::RegionalData, threshold::Float64 = 1e-4)

    for field in fieldnames(RegionalData)

        df = getfield(data, field)

        if (field != :delta_v && field != :total_use  && field != :totals)

            round_shares!(df, threshold)

            if ( any(df.eu == 1.0) || any(df.world == 1.0) )
                error("The EU and World shares must be less than 1.0")
            end
        end
    end

end

"""
    round_shares!(df::DataFrame, threshold = 1e-4)

Round shares in regional vector data

# Arguments
- `df::DataFrame`: The DataFrame containing regional data to be rounded.
- `threshold`: The threshold below which values are set to 0 (default is `1e-4`).
"""
function round_shares!(df::DataFrame, threshold = 1e-4)

    # Set values below threshold to 0
    for col in [:uk, :eu, :world]
        df[!, col] = map(x -> x < threshold ? 0 : x, df[!, col])
    end

    scaling_factor = df.uk + df.eu + df.world
    # Rescale sum to 1
    for col in [:uk, :eu, :world]
        # Avoid divide by zeros by dividing by 1.0
        replace!(scaling_factor, 0.0 => 1.0)
        df[!, col] ./= scaling_factor
    end

end

"""
    round_shares!(data::InputMatrices, threshold = 1e-4)

Round shares in regional matrix data

# Arguments
- `data::InputMatrices`: The InputMatrices struct containing regional matrix data to be rounded.
"""
function round_shares!(data::InputMatrices, threshold = 1e-4)

    # Set values below threshold to 0
    for field in [:uk, :eu, :world]

        df = getfield(data, field)
        df .= map(x -> x < threshold ? 0.0 : x, df)

    end

    scaling_factor = data.uk .+ data.eu .+ data.world
    replace!(scaling_factor, 0.0 => 1.0)

    for field in [:uk, :eu, :world]

        df = getfield(data, field)
        df ./= scaling_factor

    end

end

"""
    rescale_data!(data::CleanData)

Re-scale the data that does not get convereted into a ratio explicitly, following
<https://github.com/UCL/Supergrassi/blob/main/code/matlab/macro_v2/DataCleaning/RescaleData.m>.

# Arguments
- `data::CleanData`: The CleanData struct containing household and industry data to be rescaled.
"""
function rescale_data!(data::CleanData)

    hours_low_scaling_factor = 1.0 / mean(data.household.hours.low)
    hours_high_scaling_factor = 1.0 / mean(data.household.hours.high)
    capital_scaling_factor = 1.0 / mean(data.industry.capital.current_year)
    use_scaling_factor = 1.0 / mean(data.industry.regional.total_use.uk)
    tax_scaling_factor = 1.0 ./ (data.industry.regional.total_use.uk .- (data.industry.tax.products .+ data.industry.tax.production))

    data.household.hours.low *= hours_low_scaling_factor
    data.household.hours.high *= hours_high_scaling_factor
    data.industry.capital.next_year *= capital_scaling_factor
    data.industry.capital.current_year *= capital_scaling_factor

    data.industry.tax.products *= use_scaling_factor
    data.industry.tax.production *= use_scaling_factor
    data.household.payments.low *= use_scaling_factor
    data.household.payments.high *= use_scaling_factor

    data.household.payments.agg .*= tax_scaling_factor
    data.industry.surplus.val .*= tax_scaling_factor

    data.household.wages.low *= (use_scaling_factor / hours_low_scaling_factor)
    data.household.wages.high *= (use_scaling_factor / hours_high_scaling_factor)

    data.industry.regional.input_matrices.agg .*= tax_scaling_factor

    for col in [:uk, :eu, :world, :agg]
        data.industry.regional.delta_v[!, col] .*= use_scaling_factor
        data.industry.regional.total_use[!, col] .*= use_scaling_factor
    end

    for (i,name) in enumerate(unique(data.industry.assets_liabilities.current_year.SIC16))
        for year in [:current_year, :next_year]
            df = getfield(data.industry.assets_liabilities, year)
            k = data.industry.capital[:,year]
            mask = df.SIC16 .== name
            df.Assets[mask] .*= (k[i] / sum(df.Assets[mask]))
        end
    end

end


"""
    postprocess_clean_data!(data::CleanData)
Post-process the cleaned data to convert regional data into ratios, round shares, and rescale the data.

# Arguments
- `data::CleanData`: The CleanData struct containing household and industry data to be post-processed.
"""
function postprocess_clean_data!(data::CleanData)

    # Note: the order these are called in matters (Should be refactored)
    #       That is because the aggregate values are used to normalise the
    #       regional components, before the aggregate values themselves are normalised
    convert_to_ratio!(data.industry.regional)
    round_shares!(data.industry.regional)
    rescale_data!(data)
    convert_to_ratio!(data.household)


end
