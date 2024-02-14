using Documenter
using Catalyst, ModelingToolkit

docpath = Base.source_dir()
assetpath = joinpath(docpath, "src", "assets")
cp(joinpath(docpath, "Manifest.toml"), joinpath(assetpath, "Manifest.toml"), force = true)
cp(joinpath(docpath, "Project.toml"), joinpath(assetpath, "Project.toml"), force = true)

include("pages.jl")

# mathengine = MathJax3(Dict(:loader => Dict("load" => ["[tex]/mathtools", "[tex]/mhchem"]),
#                            :tex => Dict("inlineMath" => [["\$", "\$"], ["\\(", "\\)"]],
#                                         "packages" => [
#                                             "base",
#                                             "ams",
#                                             "autoload",
#                                             "mathtools",
#                                             "mhchem"
#                                         ])))

# makedocs(sitename = "Catalyst.jl",
#          authors = "Samuel Isaacson",
#          format = Documenter.HTML(; analytics = "UA-90474609-3",
#                                   mathengine,
#                                   prettyurls = (get(ENV, "CI", nothing) == "true"),
#                                   assets = ["assets/favicon.ico"],
#                                   canonical = "https://docs.sciml.ai/Catalyst/stable/"),
#          modules = [Catalyst, ModelingToolkit],
#          doctest = false,
#          clean = true,
#          pages = pages)

makedocs(sitename = "Catalyst.jl",
    authors = "Samuel Isaacson",
    format = Documenter.HTML(; analytics = "UA-90474609-3",
        prettyurls = (get(ENV, "CI", nothing) == "true"),
        assets = ["assets/favicon.ico"],
        canonical = "https://docs.sciml.ai/Catalyst/stable/"),
    modules = [Catalyst, ModelingToolkit],
    doctest = false,
    clean = true,
    pages = pages)

deploydocs(repo = "github.com/SciML/Catalyst.jl.git";
    push_preview = true)
