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
using RecipesBase

import Base: (==)

const ExprValues = Union{Expr,Symbol,Float64,Int}

include("reaction_network.jl")
include("maketype.jl")
include("network_properties.jl")
include("massaction_jump_utils.jl")
include("problem.jl")
include("equilibrate_utils.jl")
include("plot_recipes.jl")

# reaction network macro
export @reaction_network, @reaction_func, @min_reaction_network, @empty_reaction_network

# functions to query network properties
export species, params, speciesmap, paramsmap, numspecies, numreactions, numparams
export oderhsfun, jacfun, paramjacfun, odefun, noisefun, sdefun, jumps, regularjumps
export odeexprs, jacobianexprs, noiseexprs, jumpexprs, rateexpr, oderatelawexpr, ssaratelawexpr
export substratestoich, productstoich, netstoich, ismassaction, dependants, dependents, substrates, products, substratesymstoich, productsymstoich
export rxtospecies_depgraph, speciestorx_depgraph, rxtorx_depgraph

# functions to extend empty_ and min_ reaction_networks
export addspecies!, addparam!, add_scale_noise_param!, addreaction!

# functions to add mathematical equations to the network
export addodes!, addsdes!, addjumps!, addequi!

# problems that can be solved from the network
export ODEProblem, SDEProblem, DiscreteProblem, JumpProblem, SteadyStateProblem

# tolls for finding equilibrium solutions and bifurcation diagrams.
export EquilibrateContent
export @add_constraint, @add_constraints, internal___add___constraint, fix_parameters, @make_hc_template, make_hc_template, @add_hc_template, add_hc_template
export steady_states, stability, HcBifurcationSolver1, HcBifurcationSolver3
export bifurcations, bifurcations_grid, bifurcations_grid_2d, bifurcations_diagram_grid

end # module
