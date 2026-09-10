module CatalystCairoMakieExtension

# Fetch packages. The plotting functions used here are Makie's; `CairoMakie` is the
# backend that has to be loaded for them to render, which the extension trigger ensures.
using Catalyst: Catalyst, DiscreteSpaceReactionSystem, grid_size, spat_getu
using Makie: Makie, heatmap, record, scatterlines
import Catalyst: dspace_plot, dspace_animation, dspace_kymograph

# Creates and exports utilities for plotting discrete space simulations.
include("CatalystCairoMakieExtension/cairo_makie_extension_spatial_modelling.jl")

end
