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
    expected_filename = "data1.csv"
    expected_directory = "dir"

    created_filepath = create_filepath_from_template(template, substitutions)

    @test created_filepath.file == expected_filename
    @test created_filepath.directory == expected_directory


end