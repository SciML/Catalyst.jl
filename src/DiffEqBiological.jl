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

include("reaction_network.jl")
include("maketype.jl")
include("problem.jl")

export @reaction_network, @reaction_func
export ODEProblem, SDEProblem, JumpProblem, SteadyStateProblem

Reaction(args...) = error("""
 The old Reaction DSL is deprecated for a new
 macro-based DSL which supports parameters, regulation,
 etc. Please see the documentation for details
""")

end # module
