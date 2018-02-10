__precompile__()

module DiffEqBiological

using Reexport
using DiffEqBase
using DiffEqJump
using SymEngine
using MacroTools
using DataStructures
@reexport using DiffEqJump

using Compat
abstract type AbstractReaction end

include("ReactionNetwork.jl")
include("maketype.jl")
include("problem.jl")

export AbstractReaction

export @reaction_network, @reaction_func
export ODEProblem, SDEProblem, JumpProblem

Reaction(args...) = warn("""
The old Reaction DSL is deprecated for a new 
macro-based DSL which supports parameters, regulation,
etc. Please see the documentation for details
""")

export Reaction 

end # module
