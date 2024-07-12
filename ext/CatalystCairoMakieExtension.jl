module CatalystCairoMakieExtension

# Fetch packages.
using Catalyst, CairoMakie

# Creates and exports hc_steady_states function.
include("CatalystCairoMakieExtension/cairo_makie_extension_spatial_modelling.jl")

end
