
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

struct InputOutputData

    m_value_uk::DataFrame
    f_value_uk::DataFrame
    i_value_uk::DataFrame

    delta_v_value_uk::DataFrame

    x1_value_uk::DataFrame
    x2_value_uk::DataFrame

end

struct InputOutput


    raw_data::DataFrame

    input_output_matrix::DataFrame
    industry_names::Array{String, 1}

    final_consumption::Array{Number, 1}
    gross_fixed_capital_formation::Array{Number, 1}

    delta_v_value_uk::Array{Number, 1}

    x1_value_uk::Array{Number, 1}
    x2_value_uk::Array{Number, 1}

    y_value_uk::Array{Number, 1}

    x_services::Array{Number, 1}


    function InputOutput(raw_data::DataFrame)


        input_output_matrix = raw_data[7:111, 3:107]
        industry_names = Array(raw_data[4, 3:107])

        final_consumption = Array(raw_data[7:111 , 109])
        gross_fixed_capital_formation = Array(raw_data[7:111 , 113])

        delta_v_value_uk = Array(raw_data[7:111 , 114]) + Array(raw_data[7:111 , 115])

        x1_value_uk = Array(raw_data[7:111 , 116])
        x2_value_uk = Array(raw_data[7:111 , 117])

        y_value_uk = Array(raw_data[7:111 , 119])

        x_services = Array(raw_data[7:111 , 118])

        return new(raw_data, input_output_matrix, industry_names, final_consumption, gross_fixed_capital_formation, delta_v_value_uk, x1_value_uk, x2_value_uk, y_value_uk, x_services)

    end

end

struct Data

    household::HouseHoldData
    industry::IndustryData

    input_output::InputOutput

    depreciation::DataFrame

    risk_free_rate::DataFrame
    assets::DataFrame
    model_results::DataFrame

    merge_codes_105::DataFrame

    # gdp::DataFrame

    function Data(data_struct::Dict{String, DataFrame})

        household = HouseHoldData(
            IncomeData(data_struct["hi_income"], data_struct["lo_income"]),
            IncomeData(data_struct["hi_hours"], data_struct["lo_hours"])
        )     

        industry = IndustryData(
            data_struct["capital"],
            data_struct["turnover"],
            data_struct["inventory"]
        )

        input_output = InputOutput(data_struct["input_output"])

        depreciation = data_struct["depreciation"]
        risk_free_rate = data_struct["risk_free_rate"]
        assets = data_struct["assets"]
        model_results = data_struct["model_results"]
        merge_codes_105 = data_struct["merge_codes_105"]

        return new(household, industry, input_output, depreciation, risk_free_rate, assets, model_results, merge_codes_105)
    end

end

function organise_data(data_struct::Dict{String, DataFrame})

    data = Data(data_struct)

    return data
    
end