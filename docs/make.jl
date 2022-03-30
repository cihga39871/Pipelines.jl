#!julia --color=yes
push!(LOAD_PATH,"../src/")

using Documenter, Pipelines

makedocs(
    sitename="Pipelines.jl",
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "command_program.md",
            "julia_program.md",
            "command_dependency.md",
        ],
        "API.md",
        "Change Log" => "changelog.md"
    ]
)

deploydocs(
    repo = "github.com/cihga39871/Pipelines.jl.git",
    devbranch = "main"
)
