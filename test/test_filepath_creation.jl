using Test
using Supergrassi

@testset "Filename Templating" begin


    substitutions = Dict(
        "name" => "data",
        "format" => "csv"
    )

    template = "{name}.{format}"
    expected_filename = "data.csv"

    created_filepath = create_filepath_from_template(template, substitutions)

    println(created_filepath.file)

    @test created_filepath.file == expected_filename


end