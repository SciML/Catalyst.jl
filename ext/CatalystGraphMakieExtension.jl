module CatalystGraphMakieExtension

# Fetch packages.
using Catalyst, GraphMakie
import Catalyst: lattice_plot, lattice_animation, extract_vals
import Graphs: AbstractGraph, SimpleGraph

# Creates and exports hc_steady_states function.
include("CatalystGraphMakieExtension/graph_makie_extension_spatial_modelling.jl")

end
