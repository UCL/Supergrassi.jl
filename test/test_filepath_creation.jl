using Test
using Supergrassi

@testset "Filename Templating" begin


    substitutions = Dict(
        "name" => "data",
        "format" => "csv",
        "directory" => "dir",
        "suffix" => 1,
    )

    template = "{directory}/{name}{suffix}.{format}"

    expected_path = "dir/data1.csv"
    expected_filename = "data1.csv"
    expected_directory = "dir"
    expected_extension = ".csv"

    created_filepath = create_filepath_from_template(template, substitutions)

    @test created_filepath.path == expected_path
    @test created_filepath.file == expected_filename
    @test created_filepath.directory == expected_directory
    @test created_filepath.extension == expected_extension


end