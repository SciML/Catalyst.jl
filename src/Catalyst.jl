"""
$(DocStringExtensions.README)
"""
module Catalyst

using DocStringExtensions
using SparseArrays, DiffEqBase, Reexport
using LaTeXStrings, Latexify, Requires
using JumpProcesses: JumpProcesses, JumpProblem, MassActionJump, ConstantRateJump,
                  VariableRateJump

# ModelingToolkit imports and convenience functions we use
using ModelingToolkit
const MT = ModelingToolkit
@reexport using ModelingToolkit
using Symbolics
using ModelingToolkit: Symbolic, value, istree, get_states, get_ps, get_iv, get_systems,
                       get_eqs, get_defaults, toparam, get_var_to_name, get_observed, getvar

import ModelingToolkit: get_variables, namespace_expr, namespace_equation, get_variables!,
                        modified_states!, validate, namespace_variables,
                        namespace_parameters,
                        rename, renamespace, getname, flatten
# internal but needed ModelingToolkit functions
import ModelingToolkit: check_variables, check_parameters, _iszero, _merge, check_units,
                        get_unit

import Base: (==), hash, size, getindex, setindex, isless, Sort.defalg, length, show
import MacroTools, Graphs
import DataStructures: OrderedDict, OrderedSet
import Parameters: @with_kw_noshow

# globals for the modulate
const DEFAULT_IV = (@parameters t)[1]

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
export ODEProblem, SDEProblem, JumpProblem, NonlinearProblem, DiscreteProblem,
       SteadyStateProblem
export get_constraints, has_constraints, get_combinatoric_ratelaws

# reaction_network macro
const ExprValues = Union{Expr, Symbol, Float64, Int}
include("expression_utils.jl")
include("reaction_network.jl")
export @reaction_network, @add_reactions, @reaction

# registers CRN specific functions using Symbolics.jl
include("registered_functions.jl")
export mm, mmr, hill, hillr, hillar

# functions to query network properties
include("networkapi.jl")
export species, reactionparams, reactions, speciesmap, paramsmap, reactionparamsmap
export numspecies, numreactions, numreactionparams, setdefaults!, symmap_to_varmap
export make_empty_network, addspecies!, addparam!, addreaction!
export dependants, dependents, substoichmat, prodstoichmat, netstoichmat
export conservationlaws, conservedquantities, conservedequations, conservationlaw_constants

# depreciated functions to remove in future releases
export params, numparams

# network analysis functions
export reactioncomplexmap, reactioncomplexes, incidencemat, reactionrates, complexstoichmat
export complexoutgoingmat, incidencematgraph, linkageclasses, deficiency, subnetworks
export linkagedeficiencies, isreversible, isweaklyreversible

# for Latex printing of ReactionSystems
include("latexify_recipes.jl")

# for making and saving graphs
include("graphs.jl")
export Graph, savegraph, complexgraph

end # module
