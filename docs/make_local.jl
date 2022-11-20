#!julia --color=yes
push!(LOAD_PATH,"../src/")

using Documenter

include("../src/Pipelines.jl")
using .Pipelines

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
        "Tips & Troubleshoots" => "tips_and_troubleshoots.md",
        "Change Log" => "changelog.md"
    ]
)
