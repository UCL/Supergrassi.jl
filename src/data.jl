"""
Data structures for the model.
"""

"""
    struct RawIncomeData

Stores high and low income data as DataFrames.
"""
struct RawIncomeData
    high::DataFrame
    low::DataFrame
end

"""
    struct RawHouseHoldData

Stores household data (income and hours worked).
"""
struct RawHouseHoldData
    income::RawIncomeData
    hours::RawIncomeData
end

"""
    struct RawIndustryData

Stores raw industry data (capital, turnover, and inventory) as DataFrames.
"""
struct RawIndustryData
    capital::DataFrame
    turnover::DataFrame
    inventory::DataFrame
end


"""
    struct InputOutput

Stores input-output data matrix and other related information.

Contains:

- `raw_data`: The raw input-output data as a DataFrame.
- `input_output_matrix`: The input-output matrix as a DataFrame.
- `industry_names`: An array of industry names.
- `final_consumption`: An array of final consumption values.
- `gross_fixed_capital_formation`: An array of gross fixed capital formation values.
- `delta_v_value_uk`: An array of delta V values for the UK.
- `exports_eu_to_uk`: An array of exports from the EU to the UK.
- `export_world_to_uk`: An array of exports from the world to the UK.
- `total_use`: An array of total use values.
- `services_export`: An array of services export values.
"""
struct InputOutput


    raw_data::DataFrame

    input_output_matrix::DataFrame
    industry_names::Array{String, 1}

    final_consumption::Array{Number, 1}
    gross_fixed_capital_formation::Array{Number, 1}

    delta_v_value_uk::Array{Number, 1}

    exports_eu_to_uk::Array{Number, 1}
    export_world_to_uk::Array{Number, 1}

    total_use::Array{Number, 1}

    services_export::Array{Number, 1}


    """
        InputOutput(raw_data::DataFrame, settings::Dict{String, Any})

    Constructs an `InputOutput` object from raw data and settings.

    # Arguments
    - `raw_data::DataFrame`: The raw datafrom the input files.
    - `settings::Dict{String, Any}`: A dictionary containing settings, including limits for processing the data.

    # Returns
    - `InputOutput`: An instance of the `InputOutput` struct containing processed data.
    """
    function InputOutput(raw_data::DataFrame, settings::Dict{String, Any})

        limits = settings["excel_limits"]["input_output"]

        input_output_matrix = raw_data[limits["row_range"][1]:limits["row_range"][2], limits["matrix_cols"][1]:limits["matrix_cols"][2]]
        temp_industry_names = Array(raw_data[limits["industry_names_row"], limits["matrix_cols"][1]:limits["matrix_cols"][2]])

        industry_names = Array{String, 1}(undef, length(temp_industry_names))
        for i in eachindex(temp_industry_names)
            industry_name = temp_industry_names[i]
            industry_name = split(industry_name, "CPA_")[2]
            industry_names[i] = industry_name

        end

        input_output_matrix = DataFrame(permutedims(input_output_matrix), industry_names)

        final_consumption = Array(raw_data[limits["row_range"][1]:limits["row_range"][2], limits["final_consumption_col"]])
        gross_fixed_capital_formation = Array(raw_data[limits["row_range"][1]:limits["row_range"][2], limits["gross_fixed_capital_formation_col"]])

        delta_v_value_uk = Array(raw_data[limits["row_range"][1]:limits["row_range"][2], limits["delta_v_value_uk_col_1"]]) + Array(raw_data[limits["row_range"][1]:limits["row_range"][2], limits["delta_v_value_uk_col_2"]])

        exports_eu_to_uk = Array(raw_data[limits["row_range"][1]:limits["row_range"][2], limits["exports_eu_to_uk_col"]])
        exports_world_to_uk = Array(raw_data[limits["row_range"][1]:limits["row_range"][2], limits["exports_world_to_uk_col"]])

        total_use = Array(raw_data[limits["row_range"][1]:limits["row_range"][2], limits["total_use_col"]])

        services_export = Array(raw_data[limits["row_range"][1]:limits["row_range"][2], limits["services_export_col"]])

        return new(raw_data, input_output_matrix, industry_names, final_consumption, gross_fixed_capital_formation, delta_v_value_uk, exports_eu_to_uk, exports_world_to_uk, total_use, services_export)

    end

end


"""
    struct Data

Stores all the data required for the model, including household, industry, input-output data, and other relevant information.

# Contains:

- `household`: An instance of `RawHouseHoldData` containing household income and hours worked data.
- `industry`: An instance of `RawIndustryData` containing industry capital, turnover, and inventory data.
- `input_output`: An instance of `InputOutput` containing input-output data and matrices.
- `imports`: An instance of `InputOutput` containing import data and matrices.
- `depreciation`: A DataFrame containing depreciation data.
- `risk_free_rate`: A DataFrame containing risk-free rate data.
- `assets`: A DataFrame containing asset data.
- `model_results`: A DataFrame containing model results.
- `merge_codes_105`: A DataFrame containing merge codes for 105 industries.
- `merge_codes_64`: A DataFrame containing merge codes for 64 industries.
- `others`: A DataFrame containing other relevant data such as taxes and subsidies.
"""
struct Data

    household::RawHouseHoldData
    industry::RawIndustryData

    input_output::InputOutput
    imports::InputOutput

    depreciation::DataFrame

    risk_free_rate::DataFrame
    assets::DataFrame
    model_results::DataFrame

    merge_codes_105::DataFrame
    merge_codes_64::DataFrame

    others::DataFrame

    # gdp::DataFrame

    """
        Data(data_struct::Dict{String, DataFrame}, settings::Dict{String, Any})

    Constructs a `Data` object from a dictionary of DataFrames and settings.

    # Arguments
    - `data_struct::Dict{String, DataFrame}`: A dictionary containing DataFrames for various data categories.
    - `settings::Dict{String, Any}`: A dictionary containing settings for processing the data.

    # Returns
    - `Data`: An instance of the `Data` struct containing all the processed data.
    """
    function Data(data_struct::Dict{String, DataFrame}, settings::Dict{String, Any})

        household = RawHouseHoldData(
            RawIncomeData(data_struct["hi_income"], data_struct["lo_income"]),
            RawIncomeData(data_struct["hi_hours"], data_struct["lo_hours"])
        )     

        industry = RawIndustryData(
            data_struct["capital"],
            data_struct["turnover"],
            data_struct["inventory"]
        )

        input_output = InputOutput(data_struct["input_output"], settings)
        imports = InputOutput(data_struct["imports"], settings)

        others = data_struct["others"]

        depreciation = data_struct["depreciation"]
        risk_free_rate = data_struct["risk_free_rate"]
        assets = data_struct["assets"]
        model_results = data_struct["model_results"]
        merge_codes_105 = data_struct["merge_codes_105"]
        merge_codes_64 = data_struct["merge_codes_64"]

        return new(household, industry, input_output, imports, depreciation, risk_free_rate, assets, model_results, merge_codes_105, merge_codes_64, others)
    end

end
