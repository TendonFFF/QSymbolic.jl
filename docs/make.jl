using Documenter
using QSymbolic

# Setup for doctests - makes QSymbolic exports available
DocMeta.setdocmeta!(QSymbolic, :DocTestSetup, :(using QSymbolic); recursive=true)

makedocs(
    sitename = "QSymbolic.jl",
    modules = [QSymbolic],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://TendonFFF.github.io/QSymbolic.jl",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Guide" => [
            "Getting Started" => "guide/getting_started.md",
            "Basis Transforms" => "guide/transforms.md",
            "Composite Systems" => "guide/composite.md",
            "Operators" => "guide/operators.md",
            "Symbolic Scalars" => "guide/symbolic.md",
        ],
        "API Reference" => [
            "Spaces" => "api/spaces.md",
            "Bases" => "api/bases.md",
            "States" => "api/states.md",
            "Transforms" => "api/transforms.md",
            "Operators" => "api/operators.md",
            "Symbolic Scalars" => "api/symbolic.md",
        ],
    ],
    warnonly = [:missing_docs, :doctest, :docs_block, :cross_references],
)

deploydocs(
    repo = "github.com/TendonFFF/QSymbolic.jl.git",
    devbranch = "main",
    push_preview = true,
)
