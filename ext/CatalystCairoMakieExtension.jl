module CatalystCairoMakieExtension

# Fetch packages.
using Catalyst, CairoMakie, SparseArrays
import Catalyst: lattice_plot, lattice_animation, lattice_kymograph

# Creates and exports hc_steady_states function.
include("CatalystCairoMakieExtension/cairo_makie_extension_spatial_modelling.jl")

end
