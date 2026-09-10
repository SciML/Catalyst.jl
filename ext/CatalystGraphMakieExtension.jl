module CatalystGraphMakieExtension

# Fetch packages.
using Catalyst: Catalyst, DiscreteSpaceReactionSystem, ReactionSystem, incidencemat,
    prodstoichmat, reactioncomplexes, reactions, spat_getu, species, speciesmap,
    substoichmat
using GraphMakie: graphplot
using Graphs: Graphs, AbstractGraph, Edge, SimpleDiGraph, SimpleGraph, add_edge!, dst,
    edges, edgetype, has_edge, has_vertex, inneighbors, is_connected, is_directed, ne,
    nv, outneighbors, src, vertices
using Makie: Makie, DataAspect, hidedecorations!, hidespines!, record
using NetworkLayout: Stress
using SparseArrays: nonzeros, nzrange, rowvals
using SymbolicUtils: unwrap
using Symbolics: get_variables!, sorted_arguments
import Catalyst: species_reaction_graph, incidencematgraph, dspace_plot, dspace_animation
import SymbolicIndexingInterface: getname

# Creates and exports graph plotting functions.
include("CatalystGraphMakieExtension/graph_makie_extension_spatial_modelling.jl")
include("CatalystGraphMakieExtension/rn_graph_plot.jl")
end
