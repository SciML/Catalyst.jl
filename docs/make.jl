using Documenter
using Catalyst

makedocs(
    sitename = "Catalyst",
    format = Documenter.HTML(),
    modules = [Catalyst]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
