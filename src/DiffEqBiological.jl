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

const ExprValues = Union{Expr,Symbol,Float64,Int}                   

include("reaction_network.jl")
include("maketype.jl")
include("network_properties.jl")
include("massaction_jump_utils.jl")
include("problem.jl")

# reaction network macro
export @reaction_network, @reaction_func, @min_reaction_network

# functions to query network properties
export speciesmap, paramsmap, get_substrate_stoich, get_net_stoich

# functions to add mathematical equations to the network
export addodes!, addsdes!, addjumps!

# problems that can be solved from the network
export ODEProblem, SDEProblem, DiscreteProblem, JumpProblem, SteadyStateProblem

end # module
