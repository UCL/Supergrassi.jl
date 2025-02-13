using CSV
using DataFrames
using XLSX

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


"""
    Reads a data file and returns a DataFrame.

    # Arguments
    - `file::String`: The path to the file.

    # Examples
    ```julia
    df = read_data("data.csv")
    ```
"""
function read_data(file::String)
    
    if !isfile(file)
        error("File not found: $file")
    end

    if occursin(r"\.csv$", file)
        return read_csv(file)
    else
        error("Invalid file format: $file")
    end

end

"""
    Reads a data file and returns a DataFrame.

    # Arguments
    - `file::String`: The path to the file.
    - `sheet::String`: The name of the sheet.
    - `range::String`: The range of cells to read.

    # Examples
    ```julia
    df = read_data("data.xlsx", "Sheet1", "A1:Z1000")
    ```
"""
function read_data(file::String, sheet::String, range::String)

    if !isfile(file)
        error("File not found: $file")
    end

    if !occursin(r"\.xlsx$", file)
        error("Invalid file format: $file")
    end

    return read_excel(file, sheet=sheet, range=range)

end

function read_data(filepaths::Dict{String, FilePath})
    data = Dict{String, DataFrame}()

    for (key, filepath) in filepaths
        println("Reading data from $(filepath.path)")

        if occursin(r"\.csv$", filepath.path)
            data[key] = read_csv(filepath.path)
        elseif occursin(r"\.xlsx$", filepath.path)
            @warn "Excel files not supported yet"
        else
            @warn "Invalid file format: $(filepath.path)"
        end
    end

    return data
end