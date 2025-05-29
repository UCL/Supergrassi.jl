using Documenter, Supergrassi

makedocs(
    modules = [Supergrassi],
    sitename = "Supergrassi",
    pages    = [
        "Introduction" => "index.md",
    ],
    format = Documenter.HTML(
        ;
        prettyurls = get(ENV, "CI", "false") == "true",
    ),
)

deploydocs(
    repo = "github.com/UCL/Supergrassi.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
    push_preview = true,
)
