
struct IncomeData
    high::DataFrame
    low::DataFrame
end

struct HouseHoldData
    income::IncomeData
    hours::IncomeData
end

struct IndustryData
    capital::DataFrame
    turnover::DataFrame
    inventory::DataFrame
end


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

        input_output_matrix = DataFrame(input_output_matrix, industry_names)

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

struct Data

    household::HouseHoldData
    industry::IndustryData

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

    function Data(data_struct::Dict{String, DataFrame}, settings::Dict{String, Any})

        household = HouseHoldData(
            IncomeData(data_struct["hi_income"], data_struct["lo_income"]),
            IncomeData(data_struct["hi_hours"], data_struct["lo_hours"])
        )     

        industry = IndustryData(
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
