using CSV
using DataFrames
using XLSX
using Logging

function read_excel(file::String; sheet::String="Sheet1", range::String="A1:Z1000")

    range_regex = r"^([A-Z]{1,3}[0-9]{1,7}(:[A-Z]{1,3}[0-9]{1,7})?|[A-Z]{1,3}:[A-Z]{1,3}|[0-9]{1,7}:[0-9]{1,7})$"
    if !occursin(range_regex, range)
        error("Invalid range: $range")
    end

    try
        XLSX.openxlsx(file, enable_cache=true) do f
            sheet_data = f[sheet]
            data = sheet_data[range]
            return DataFrames.DataFrame(data, :auto)
        end
    catch e

        if isa(e, ErrorException)
            error("Sheet not found: $sheet")
        else
            error("An error occurred: $e")
        end
    end

end


function read_csv(file::String)
    return CSV.read(file, DataFrames.DataFrame)
end


function read_data(filepaths::Dict{String, FilePath}, settings::Dict{String, Any})
    data = Dict{String, DataFrame}()

    for (key, filepath) in filepaths
        @debug "Reading data from $(filepath.path)"

        if endswith(filepath.path, ".csv")
            data[key] = read_csv(filepath.path)
        elseif endswith(filepath.path, ".xlsx")
            sheet = settings["excel_limits"][key]["sheet"]
            top_left = settings["excel_limits"][key]["top_left"]
            bottom_right = settings["excel_limits"][key]["bottom_right"]
            range = "$top_left:$bottom_right"
            data[key] = read_excel(filepath.path, sheet=sheet, range=range)
        else
            @warn "Invalid file format: $(filepath.path)"
        end
    end

    data = Data(data, settings)

    @info "Data read successfully."

    return data
end
