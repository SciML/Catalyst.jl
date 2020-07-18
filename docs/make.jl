using Documenter
using Catalyst, ModelingToolkit

makedocs(
    sitename = "Catalyst.jl",
    authors = "Samuel Isaacson",
    format = Documenter.HTML(prettyurls = (get(ENV, "CI", nothing) == "true")),
    modules = [Catalyst,ModelingToolkit],
    doctest = false,
    clean = true,
    pages = Any[
        "Home" => "index.md"#,
        "API" => Any[
            "api/catalyst_api.md"
        ]
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
