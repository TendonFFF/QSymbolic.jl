using Documenter
using QSymbolic

makedocs(
    sitename = "QSymbolic.jl",
    modules = [QSymbolic],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://YOUR_USERNAME.github.io/QSymbolic.jl",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Guide" => [
            "Getting Started" => "guide/getting_started.md",
            "Basis Transforms" => "guide/transforms.md",
            "Composite Systems" => "guide/composite.md",
        ],
        "API Reference" => [
            "Spaces" => "api/spaces.md",
            "Bases" => "api/bases.md",
            "States" => "api/states.md",
            "Transforms" => "api/transforms.md",
        ],
    ],
    warnonly = [:missing_docs],
)

deploydocs(
    repo = "github.com/YOUR_USERNAME/QSymbolic.jl.git",
    devbranch = "main",
    push_preview = true,
)
