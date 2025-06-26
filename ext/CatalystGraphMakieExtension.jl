module CatalystGraphMakieExtension

# Fetch packages.
using Catalyst, GraphMakie, Graphs, Symbolics, SparseArrays, NetworkLayout, Makie
using Symbolics: get_variables!
using Latexify, LaTeXStrings
import Catalyst: species_reaction_graph, incidencematgraph, lattice_plot, lattice_animation

# Creates and exports graph plotting functions.
include("CatalystGraphMakieExtension/graph_makie_extension_spatial_modelling.jl")
include("CatalystGraphMakieExtension/rn_graph_plot.jl")
end
