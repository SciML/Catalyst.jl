using Documenter
using Catalyst, ModelingToolkit
# Add packages for plotting
using GraphMakie, CairoMakie
#using MultiDocumenter

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
        collapselevel = 1,
        assets = ["assets/favicon.ico"],
        canonical = "https://docs.sciml.ai/Catalyst/stable/"),
    modules = [Catalyst, ModelingToolkit,
        isdefined(Base, :get_extension) ?
        Base.get_extension(Catalyst, :CatalystGraphMakieExtension) :
        Catalyst.CatalystGraphMakieExtension],
    doctest = false,
    clean = true,
    pages = pages,
    pagesonly = true,
    warnonly = [:missing_docs])

#clonedir = mktempdir()

#docs = [MultiDocumenter.MultiDocRef(
#        upstream = joinpath(@__DIR__, "build"),
#        path = "docs",
#        name = "Catalyst",
#        fix_canonical_url = false,
#        ),
#        MultiDocumenter.MultiDocRef(
#        upstream = joinpath(clonedir, "build"),
#        path = "NetworkAnalysis",
#        name = "CatalystNetworkAnalysis",
#        giturl = "https://github.com/SciML/CatalystNetworkAnalysis.jl",
#        )]

deploydocs(repo = "github.com/SciML/Catalyst.jl.git";
    push_preview = true)

#outpath = joinpath(@__DIR__, "build")
#MultiDocumenter.make(outpath, docs;
#                     assets_dir = "docs/src/assets",
#                     search_engine = MultiDocumenter.SearchConfig(index_versions = [
#                                                                      "stable",
#                                                                  ],
#                                                                  engine = MultiDocumenter.FlexSearch),
#                     brand_image = MultiDocumenter.BrandImage("https://sciml.ai",
#                                                              joinpath("assets",
#                                                                       "logo.png")))
#
