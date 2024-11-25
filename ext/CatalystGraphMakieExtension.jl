module CatalystGraphMakieExtension

# Fetch packages.
using Catalyst, GraphMakie, Graphs, Symbolics, SparseArrays
using Symbolics: get_variables!
import Catalyst: species_reaction_graph, incidencematgraph, lattice_plot, lattice_animation
import SparseArrays: rowvals, nonzeros, nzrange

# Creates and exports graph plotting functions.
include("CatalystGraphMakieExtension/graph_makie_extension_spatial_modelling.jl")
include("CatalystGraphMakieExtension/rn_graph_plot.jl")
end
