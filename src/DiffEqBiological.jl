__precompile__()

module DiffEqBiological

using Reexport
using DiffEqBase
using DiffEqJump
using SymEngine
using MacroTools
using DataStructures
using Parameters
@reexport using DiffEqBase, DiffEqJump
using Compat
@reexport using DynamicPolynomials
using HomotopyContinuation
using LinearAlgebra
using Plots

const ExprValues = Union{Expr,Symbol,Float64,Int}

include("reaction_network.jl")
include("maketype.jl")
include("network_properties.jl")
include("massaction_jump_utils.jl")
include("problem.jl")
include("equilibrate_utils.jl")

# reaction network macro
export @reaction_network, @reaction_func, @min_reaction_network

# functions to query network properties
export speciesmap, paramsmap, numspecies, numreactions, numparams
export odeexprs, jacobianexprs, noiseexprs, jumpexprs
export get_substrate_stoich, get_net_stoich
export rxtospecies_depgraph, speciestorx_depgraph, rxtorx_depgraph

# functions to add mathematical equations to the network
export addodes!, addsdes!, addjumps!, addequi!, manage_reaction_network!

# problems that can be solved from the network
export ODEProblem, SDEProblem, DiscreteProblem, JumpProblem, SteadyStateProblem

# tolls for finding equilibrium solutions and bifurcation diagrams.
export @add_constraint, @add_constraints, internal___add___constraint, fix_parameters, @make_hc_template, make_hc_template, steady_states, stability
export bifurcations, bifurcations_grid, bifurcations_grid_2d, bifurcations_diagram_grid
export bif_plot, bif_plot!, bif_scatter, bif_scatter!

end # module
