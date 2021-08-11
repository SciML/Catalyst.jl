"""
$(DocStringExtensions.README)
"""
module Catalyst

using DocStringExtensions
using DiffEqBase, Reexport, ModelingToolkit, DiffEqJump
using ModelingToolkit: Symbolic, value, istree, get_states, get_ps, get_iv, get_systems, get_eqs, get_defaults, toparam
import ModelingToolkit: get_variables, namespace_expr
import ModelingToolkit: namespace_equation, get_variables!, modified_states!

# internal but needed ModelingToolkit functions
import ModelingToolkit: check_variables, check_parameters, _iszero, _merge

const DEFAULT_IV = (@parameters t)[1]
@reexport using ModelingToolkit
import MacroTools
import Base: (==), merge!, merge, hash, size, getindex, setindex, isless, Sort.defalg, length
using Symbolics
using Latexify, Requires
import AbstractAlgebra

# as used in Catlab
const USE_GV_JLL = Ref(false)
function __init__()
    @require Graphviz_jll="3c863552-8265-54e4-a6dc-903eb78fde85" begin
      USE_GV_JLL[] = true
      let cfg = joinpath(Graphviz_jll.artifact_dir, "lib", "graphviz", "config6")
        if !isfile(cfg)
          Graphviz_jll.dot(path -> run(`$path -c`))
        end
      end
    end
  end

# base system type and features
include("reactionsystem.jl")
export Reaction, ReactionSystem, ismassaction, oderatelaw, jumpratelaw
export ODEProblem, SDEProblem, JumpProblem, NonlinearProblem, DiscreteProblem, SteadyStateProblem

const ExprValues = Union{Expr,Symbol,Float64,Int}

include("expression_utils.jl")
include("reaction_network.jl")

# reaction network macro
export @reaction_network, @add_reactions

# registers CRN specific functions using Symbolics.jl
include("registered_functions.jl")

# functions to query network properties
include("networkapi.jl")
export species, params, reactions, speciesmap, paramsmap, numspecies, numreactions, numparams
export make_empty_network, addspecies!, addparam!, addreaction!
export dependants, dependents, substoichmat, prodstoichmat, netstoichmat
export conservationlaws, conservedquantities
export reaction_complexes, reaction_rates, complex_stoich_matrix, complex_incidence_matrix, complex_outgoing_matrix
    
# for Latex printing of ReactionSystems
include("latexify_recipes.jl")

# for making and saving graphs
import Base.Iterators: flatten
import DataStructures: OrderedDict, OrderedSet
import Parameters: @with_kw_noshow
include("graphs.jl")
export Graph, savegraph

end # module
