using Documenter
using Catalyst, ModelingToolkit


# supposed hack to get mhchem but doesn't seem to work yet...
# const katex_version = "0.11.1"
# function Documenter.Writers.HTMLWriter.RD.mathengine!(r::Documenter.Utilities.JSDependencies.RequireJS, engine::Documenter.KaTeX)
#     push!(r, Documenter.Utilities.JSDependencies.RemoteLibrary(
#         "katex",
#         "https://cdnjs.cloudflare.com/ajax/libs/KaTeX/$(katex_version)/katex.min.js"
#     ))
#     push!(r, Documenter.Utilities.JSDependencies.RemoteLibrary(
#         "mhchem",
#         "https://cdnjs.cloudflare.com/ajax/libs/KaTeX/$(katex_version)/contrib/mhchem.min.js"
#     ))
#     push!(r, Documenter.Utilities.JSDependencies.RemoteLibrary(
#         "katex-auto-render",
#         "https://cdnjs.cloudflare.com/ajax/libs/KaTeX/$(katex_version)/contrib/auto-render.min.js",
#         deps = ["katex"],
#     ))
#     push!(r, Documenter.Utilities.JSDependencies.Snippet(
#         ["jquery", "katex", "mhchem", "katex-auto-render"],
#         ["\$", "katex", "mhchem", "renderMathInElement"],
#         """
#         \$(document).ready(function() {
#           renderMathInElement(
#             document.body,
#             $(Documenter.Utilities.JSDependencies.json_jsescape(engine.config, 2))
#           );
#         })
#         """
#     ))
# end

makedocs(
    sitename = "Catalyst.jl",
    authors = "Samuel Isaacson",
    format = Documenter.HTML(mathengine=Documenter.Writers.HTMLWriter.MathJax(), prettyurls = (get(ENV, "CI", nothing) == "true")),
    modules = [Catalyst,ModelingToolkit],
    doctest = false,
    clean = true,
    pages = Any[
        "Home" => "index.md",
        "Tutorials" => Any[
            "tutorials/using_catalyst.md",
            "tutorials/basics.md",
            "tutorials/models.md",
            "tutorials/basic_examples.md",
            "tutorials/advanced.md",
            "tutorials/generated_systems.md",
            "tutorials/advanced_examples.md",
            "tutorials/bifurcation_diagram.md",
            "tutorials/parameter_estimation.md"
            "tutorials/generating_reactions_programmatically.md"
        ],
        "API" => Any[
            "api/catalyst_api.md"
        ]
    ]
)

deploydocs(
   repo = "github.com/SciML/Catalyst.jl.git";
   push_preview = true
)