"""
$(DocStringExtensions.README)
"""
module Catalyst

using DocStringExtensions
using SparseArrays, DiffEqBase, Reexport, Setfield
using LaTeXStrings, Latexify, Requires
using JumpProcesses: JumpProcesses,
                     JumpProblem, MassActionJump, ConstantRateJump,
                     VariableRateJump

# ModelingToolkit imports and convenience functions we use
using ModelingToolkit
const MT = ModelingToolkit
using DynamicQuantities#, Unitful # Having Unitful here as well currently gives an error.


@reexport using ModelingToolkit
using Symbolics
using LinearAlgebra
using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

import Symbolics: BasicSymbolic
import SymbolicUtils
using ModelingToolkit: Symbolic, value, istree, get_unknowns, get_ps, get_iv, get_systems,
                       get_eqs, get_defaults, toparam, get_var_to_name, get_observed,
                       getvar

import ModelingToolkit: get_variables, namespace_expr, namespace_equation, get_variables!,
                        modified_unknowns!, validate, namespace_variables,
                        namespace_parameters, rename, renamespace, getname, flatten,
                        is_alg_equation, is_diff_equation

# internal but needed ModelingToolkit functions
import ModelingToolkit: check_variables,
                        check_parameters, _iszero, _merge, check_units,
                        get_unit, check_equations, iscomplete

import Base: (==), hash, size, getindex, setindex, isless, Sort.defalg, length, show
import MacroTools, Graphs
import Graphs: DiGraph, SimpleGraph, SimpleDiGraph, vertices, edges, add_vertices!, nv, ne
import DataStructures: OrderedDict, OrderedSet
import Parameters: @with_kw_noshow
import Symbolics: occursin, wrap

# globals for the modulate
function default_time_deriv()
    return ModelingToolkit.D_nounits
end
function default_t()
    return ModelingToolkit.t_nounits
end
const DEFAULT_IV = default_t()
const DEFAULT_IV_SYM = Symbol(DEFAULT_IV)
export default_t, default_time_deriv

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

### Package Constants ###

# Union type of types that can occur in expressions.
const ExprValues = Union{Expr, Symbol, Float64, Int, Bool}

# The symbol used for conserved quantities in conservation law eliminations.
const CONSERVED_CONSTANT_SYMBOL = :Γ

# Declares symbols which may neither be used as parameters nor unknowns.
const forbidden_symbols_skip = Set([:ℯ, :pi, :π, :t, :∅])
const forbidden_symbols_error = union(Set([:im, :nothing, CONSERVED_CONSTANT_SYMBOL]),
                                      forbidden_symbols_skip)
const forbidden_variables_error = let
    fvars = copy(forbidden_symbols_error)
    delete!(fvars, :t)
    fvars
end

### Package Main ###

# The `Reaction` structure and its functions.
include("reaction.jl")
export isspecies
export Reaction

# The `ReactionSystem` structure and its functions.
include("reactionsystem.jl")
export ReactionSystem, isspatial
export species, nonspecies, reactions, nonreactions, speciesmap, paramsmap
export numspecies, numreactions, setdefaults!
export make_empty_network
export dependants, dependents, substoichmat, prodstoichmat, netstoichmat
export isautonomous
export reactionrates
export isequivalent
export set_default_noise_scaling

# depreciated functions to remove in future releases
export params, numparams

# Conversions of the `ReactionSystem` structure.
include("reactionsystem_conversions.jl")
export ODEProblem,
       SDEProblem, JumpProblem, NonlinearProblem, DiscreteProblem,
       SteadyStateProblem
export ismassaction, oderatelaw, jumpratelaw
export symmap_to_varmap

# reaction_network macro
include("expression_utils.jl")
include("dsl.jl")
export @reaction_network, @network_component, @reaction, @species

# Network analysis functionality.
include("network_analysis.jl")
export reactioncomplexmap, reactioncomplexes, incidencemat
export complexstoichmat
export complexoutgoingmat, incidencematgraph, linkageclasses, deficiency, subnetworks
export linkagedeficiencies, isreversible, isweaklyreversible
export conservationlaws, conservedquantities, conservedequations, conservationlaw_constants

# registers CRN specific functions using Symbolics.jl
include("registered_functions.jl")
export mm, mmr, hill, hillr, hillar

# functions to query network properties

# for Latex printing of ReactionSystems
include("latexify_recipes.jl")

# for making and saving graphs
include("graphs.jl")
export Graph, savegraph, complexgraph

# for creating compounds
include("chemistry_functionality.jl")
export @compound, @compounds
export iscompound, components, coefficients, component_coefficients
export balance_reaction, balance_system

# Functionality for computing the stability of system steady states.
include("steady_state_stability.jl")
export steady_state_stability, steady_state_jac

### Extensions ###

# HomotopyContinuation
function hc_steady_states end
export hc_steady_states

# StructuralIdentifiability
function make_si_ode end
export make_si_ode

### Spatial Reaction Networks ###

# spatial reactions
include("spatial_reaction_systems/spatial_reactions.jl")
export TransportReaction, TransportReactions, @transport_reaction
export isedgeparameter

# lattice reaction systems
include("spatial_reaction_systems/lattice_reaction_systems.jl")
export LatticeReactionSystem
export spatial_species, vertex_parameters, edge_parameters

# variosu utility functions
include("spatial_reaction_systems/utility.jl")

# spatial lattice ode systems.
include("spatial_reaction_systems/spatial_ODE_systems.jl")

end # module
