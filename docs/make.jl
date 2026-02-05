using Documenter
using Catalyst, ModelingToolkitBase
# Add packages for plotting
using GraphMakie, CairoMakie

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
#          modules = [Catalyst, ModelingToolkitBase],
#          doctest = false,
#          clean = true,
#          pages = pages)

makedocs(sitename = "Catalyst.jl",
    authors = "Samuel Isaacson",
    format = Documenter.HTML(; analytics = "UA-90474609-3",
        prettyurls = (get(ENV, "CI", nothing) == "true"),
        collapselevel = 1,
        assets = ["assets/favicon.ico"],
        canonical = "https://docs.sciml.ai/Catalyst/stable/"),
    modules = [Catalyst, ModelingToolkitBase,
        isdefined(Base, :get_extension) ?
        Base.get_extension(Catalyst, :CatalystGraphMakieExtension) :
        Catalyst.CatalystGraphMakieExtension],
    doctest = false,
    clean = true,
    pages = pages,
    pagesonly = true,
    warnonly = [:missing_docs, :cross_references]) # `:cross_references` here temporarily while getting docs to work on v16.

deploydocs(repo = "github.com/SciML/Catalyst.jl.git";
    push_preview = true)
