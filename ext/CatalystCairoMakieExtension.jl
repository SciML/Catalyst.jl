module CatalystCairoMakieExtension

# Fetch packages.
using Catalyst, CairoMakie, SparseArrays
import Catalyst: lattice_plot, lattice_animation, lattice_kymograph,
                 demask_vals, extract_vals, extract_grid_axes

# Creates and exports utilities for plotting lattice simulations.
include("CatalystCairoMakieExtension/cairo_makie_extension_spatial_modelling.jl")

end
