module CatalystGraphMakieExtension

# Fetch packages.
using Catalyst: Catalyst, species_reaction_graph, incidencematgraph,
                lattice_plot, lattice_animation
using GraphMakie: GraphMakie
using Graphs: Graphs
using Symbolics: Symbolics, get_variables!
using SparseArrays: SparseArrays
using NetworkLayout: NetworkLayout
using Makie: Makie

# Creates and exports graph plotting functions.
include("CatalystGraphMakieExtension/graph_makie_extension_spatial_modelling.jl")
include("CatalystGraphMakieExtension/rn_graph_plot.jl")
end
