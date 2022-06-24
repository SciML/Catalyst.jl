using Documenter
using Catalyst, ModelingToolkit

include("pages.jl")

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

# mathengine = MathJax3(Dict(
#     :loader => Dict("load" => ["[tex]/require"])) #,
#     #:tex => Dict("packages" => ["base", "ams", "autoload", "physics","[+]","require"] ))
#     )
# mathengine = MathJax3()

mathengine = MathJax3(Dict(:loader => Dict("load" => ["[tex]/require", "[tex]/mathtools"]),
                           :tex => Dict("inlineMath" => [["\$", "\$"], ["\\(", "\\)"]],
                                        "packages" => [
                                            "base",
                                            "ams",
                                            "autoload",
                                            "mathtools",
                                            "require",
                                        ])))

makedocs(sitename = "Catalyst.jl",
         authors = "Samuel Isaacson",
         format = Documenter.HTML(; analytics = "UA-90474609-3",
                                  mathengine,
                                  prettyurls = (get(ENV, "CI", nothing) == "true")),
         modules = [Catalyst, ModelingToolkit],
         doctest = false,
         clean = true,
         pages = pages)

deploydocs(repo = "github.com/SciML/Catalyst.jl.git";
           push_preview = true)
