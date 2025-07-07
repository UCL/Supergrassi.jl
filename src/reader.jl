using CSV
using DataFrames
using XLSX
using Logging

"""
    read_excel(file::String; sheet::String="Sheet1", range::String="A1:Z1000")

Reads an Excel file and returns a DataFrame.

# Arguments
- `file::String`: The path to the Excel file.
- `sheet::String`: The name of the sheet to read. Default is "Sheet1".
- `range::String`: The range of cells to read, in Excel format (e.g., "A1:Z1000"). Default is "A1:Z1000".

# Returns
- `DataFrame`: A DataFrame containing the data from the specified sheet and range.

# Examples
```julia
df = read_excel("data.xlsx", sheet="Sheet1", range="A1:Z1000")
```
"""
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


"""
    read_csv(file::String)

Reads a CSV file and returns a DataFrame.

# Arguments
- `file::String`: The path to the CSV file.

# Returns
- `DataFrame`: A DataFrame containing the data from the CSV file.

# Examples
```julia
df = read_csv("data.csv")
```
"""
function read_csv(file::String)
    return CSV.read(file, DataFrames.DataFrame)
end


"""
    read_data(file::String)

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

    if endswith(file, ".csv")
        return read_csv(file)
    else
        error("Invalid file format: $file")
    end

end

"""
    read_data(file::String, sheet::String, range::String)

Reads a data file and returns a DataFrame.

# Arguments
- `file::String`: The path to the file.
- `sheet::String`: The name of the sheet.
- `range::String`: The range of cells to read.

# Returns
- `DataFrame`: A DataFrame containing the data from the specified sheet and range.

# Examples
```julia
df = read_data("data.xlsx", "Sheet1", "A1:Z1000")
```
"""
function read_data(file::String, sheet::String, range::String)

    if !isfile(file)
        error("File not found: $file")
    end

    if !endswith(file, ".xlsx")
        error("Invalid file format: $file")
    end

    return read_excel(file, sheet=sheet, range=range)

end

"""
    read_data(filepaths::Dict{String, FilePath}, settings::Dict{String, Any})

Reads multiple data files and returns a `Data` object.

# Arguments
- `filepaths::Dict{String, FilePath}`: A dictionary where keys are identifiers and values are `FilePath` objects.
- `settings::Dict{String, Any}`: A dictionary containing settings, including limits for processing the data.

# Returns
- `Data`: An instance of the `Data` struct containing processed data from the files.

# Examples
```julia
data_files = Dict(
    "input_output" => FilePath("data/input_output.xlsx"),
    "assets_liabilities" => FilePath("data/assets_liabilities.csv")
)
```
"""
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
