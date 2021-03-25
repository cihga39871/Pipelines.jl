#!julia --color=yes
push!(LOAD_PATH,"../src/")

using Documenter, Pipelines

makedocs(
    sitename="Pipelines.jl",
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "command_program.md",
            "command_dependency.md",
        ],
        "API.md"
    ]
)
