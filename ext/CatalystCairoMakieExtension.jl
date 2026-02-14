module CatalystCairoMakieExtension

# Fetch packages.
using Catalyst, CairoMakie, SparseArrays
import Catalyst: dspace_plot, dspace_animation, dspace_kymograph,
                 demask_vals, extract_vals, extract_grid_axes

# Creates and exports utilities for plotting discrete space simulations.
include("CatalystCairoMakieExtension/cairo_makie_extension_spatial_modelling.jl")

end
