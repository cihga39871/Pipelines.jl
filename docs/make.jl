#!julia --color=yes

using Documenter, Pipelines

# include("../src/Pipelines.jl")
# using .Pipelines

makedocs(
    sitename="Pipelines.jl",
    authors = "Dr. Jiacheng Chuan, and contributors.",
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
    ],
    format = Documenter.HTML(
        sidebar_sitename=false,
        assets = ["assets/favicon.ico"]
    )
)

if haskey(ENV, "GITHUB_TOKEN")
    deploydocs(
        repo = "github.com/cihga39871/Pipelines.jl.git",
        devbranch = "main"
    )
end