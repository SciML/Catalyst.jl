module CatalystGraphMakieExtension

# Fetch packages.
using Catalyst, GraphMakie, Graphs 
import Catalyst: lattice_plot, lattice_animation, extract_vals
import Graphs: AbstractGraph, SimpleGraph, SimpleDiGraph, SimpleEdge, src, dst, ne, nv 
import Catalyst: speciesreactiongraph, incidencematgraph

# Creates and exports hc_steady_states function.
include("CatalystGraphMakieExtension/graph_makie_extension_spatial_modelling.jl")
include("CatalystGraphMakieExtension/rn_graph_plot.jl")
export SRGraph, ComplexGraph
end
