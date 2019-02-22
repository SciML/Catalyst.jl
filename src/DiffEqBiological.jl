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

include("reaction_network.jl")
include("maketype.jl")
include("massaction_jump_utils.jl")
include("problem.jl")
include("equilibrate_utils.jl")

export @reaction_network, @reaction_func, @min_reaction_network
export addodes!, addsdes!, addjumps!, addequi!
export ODEProblem, SDEProblem, DiscreteProblem, JumpProblem, SteadyStateProblem
export @fixed_concentration, internal___fix___concentrations, fix_parameters, @make_hc_template, make_hc_template, steady_states, stability
export bifurcations, bifurcations_grid, bifurcations_grid_2d, bifurcations_diagram_grid
export bif_plot, bif_plot!, bif_scatter, bif_scatter!

end # module
