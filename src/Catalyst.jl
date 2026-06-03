"""
$(DocStringExtensions.README)
"""
module Catalyst

using DocStringExtensions
using SparseArrays, DiffEqBase, Reexport, Setfield, EnumX
import SciMLBase
using LaTeXStrings, Latexify
using LinearAlgebra, Combinatorics
using JumpProcesses: JumpProcesses, JumpProblem,
                     MassActionJump, ConstantRateJump, VariableRateJump,
                     SpatialMassActionJump, CartesianGrid, CartesianGridRej

# ModelingToolkit imports and convenience functions we use
using ModelingToolkitBase
const MT = ModelingToolkitBase
using DynamicQuantities

@reexport using ModelingToolkitBase
using Symbolics
using LinearAlgebra
using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

import Symbolics: SymbolicT
using Symbolics: iscall, sorted_arguments, value
using ModelingToolkitBase: get_unknowns, get_ps, get_iv, get_systems,
                       get_eqs, toparam, get_var_to_name, get_observed,
                       getvar, has_iv, JumpType

import ModelingToolkitBase: get_variables, namespace_expr, namespace_equation,
                        modified_unknowns!, namespace_variables,
                        namespace_parameters, renamespace, flatten,
                        is_alg_equation, is_diff_equation, collect_vars!,
                        eqtype_supports_collect_vars

# Import from owner modules (not re-exporters) per ExplicitImports.jl audit
import Symbolics: get_variables!, rename
import SymbolicIndexingInterface
import SymbolicIndexingInterface: getname
import ModelingToolkitBase: SymmapT

# internal but needed ModelingToolkit functions
import ModelingToolkitBase: check_variables, check_parameters,
                        check_equations, iscomplete

# Import from owner module (SymbolicUtils) per ExplicitImports.jl audit
import SymbolicUtils: _iszero, unwrap

import Base: (==), hash, size, getindex, setindex, isless, Sort.defalg, length, show
import MacroTools, Graphs
using MacroTools: striplines
import Graphs: DiGraph, SimpleGraph, SimpleDiGraph, vertices, edges, add_vertices!, nv, ne
import DataStructures: OrderedDict, OrderedSet
import Parameters: @with_kw_noshow
# Note: occursin is from Base (not Symbolics), so we don't import it
import Symbolics: wrap
import Symbolics.RewriteHelpers: hasnode, replacenode
import SymbolicUtils: getmetadata, hasmetadata, setmetadata
import SciMLPublic: @public

# globals for the modulate
function default_time_deriv()
    return ModelingToolkitBase.D_nounits
end
function default_t()
    return ModelingToolkitBase.t_nounits
end
const DEFAULT_IV = default_t()
const DEFAULT_IV_SYM = Symbol(DEFAULT_IV)
export default_t, default_time_deriv

### Package Constants ###

# Union type of types that can occur in expressions.
const ExprValues = Union{Expr, Symbol, Float64, Int, Bool}

# The symbol used for conserved quantities in conservation law eliminations.
const CONSERVED_CONSTANT_SYMBOL = :Γ

# Declares symbols which may neither be used as parameters nor unknowns.
const forbidden_symbols_skip = Set([:ℯ, :pi, :π, :t, :∅, :Ø])
const forbidden_symbols_error = union(Set([:im, :nothing, CONSERVED_CONSTANT_SYMBOL]),
    forbidden_symbols_skip)

### Unit Helpers ###

# SymbolicDimensions-preserving unit inference (replaces MTKBase's `get_unit` for validation).
include("unit_helpers.jl")

### LaTeX Utilities ###

# Accessor functions for Symbolics' SymLatexWrapper metadata.
include("latex_utils.jl")

### Package Main ###

# The `Reaction` structure and its functions.
include("reaction.jl")
export isspecies
export Reaction, PhysicalScale

# Union type for `Reaction`s and `Equation`s.
const CatalystEqType = Union{Reaction, Equation}

# The `ReactionSystem` structure and its functions.
include("reactionsystem.jl")
export ReactionSystem, isspatial
export species, nonspecies, reactions, nonreactions, speciesmap, paramsmap
export numspecies, numreactions, numparams
export make_empty_network
export dependants, dependents, substoichmat, prodstoichmat, netstoichmat
export isautonomous
export reactionrates
export set_default_noise_scaling
export ode_model, sde_model, jump_model, ss_ode_model, hybrid_model

# Mark unit validation APIs as public without exporting them.
@public validate_units, assert_valid_units, unit_validation_report
@public UnitValidationError, UnitValidationIssue, UnitValidationReport

# System-level metadata key types and accessors.
include("reactionsystem_metadata.jl")
@public U0Map, ParameterMap
@public has_u0_map, get_u0_map, set_u0_map
@public has_parameter_map, get_parameter_map, set_parameter_map

# Conversions of the `ReactionSystem` structure.
include("reactionsystem_conversions.jl")
export ODEProblem, SDEProblem, JumpProblem, NonlinearProblem,
       SteadyStateProblem, HybridProblem
export ismassaction, oderatelaw, jumpratelaw

# reaction_network macro
include("expression_utils.jl")
include("dsl.jl")
export @reaction_network, @network_component, @reaction, @species

# Network analysis functionality.
include("network_analysis.jl")
export reactioncomplexmap, reactioncomplexes, incidencemat
export complexstoichmat, laplacianmat, fluxmat, massactionvector, complexoutgoingmat,
       adjacencymat
export incidencematgraph, linkageclasses, stronglinkageclasses,
       terminallinkageclasses, deficiency, subnetworks
export linkagedeficiencies, isreversible, isweaklyreversible
export conservationlaws, conservedquantities, conservedequations, conservationlaw_constants
export satisfiesdeficiencyone, satisfiesdeficiencyzero
export iscomplexbalanced, isdetailedbalanced, robustspecies

# Containes the `nullspace` function required for conservation law elimination.
include("mtk_nullspace_function.jl")

# registers CRN specific functions using Symbolics.jl
include("registered_functions.jl")
export mm, mmr, hill, hillr, hillar

# functions to query network properties

# for Latex printing of ReactionSystems
include("latexify_recipes.jl")

# for making and saving graphs/plots
include("plotting.jl")

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

# GraphMakie: functionality for plotting species-reaction graphs and complexes
function plot_network end
function plot_complexes end
export plot_network, plot_complexes

### Spatial Reaction Networks ###

# Spatial reactions.
include("spatial_reaction_systems/spatial_reactions.jl")
export TransportReaction, @transport_reaction
export isedgeparameter

# Lattice reaction systems.
include("spatial_reaction_systems/discrete_space_reaction_systems.jl")
export DiscreteSpaceReactionSystem
export spatial_species, vertex_parameters, edge_parameters
export CartesianGrid, CartesianGridReJ # (Implemented in JumpProcesses)
export has_cartesian_dspace, has_masked_dspace, has_grid_dspace, has_graph_dspace,
       grid_dims, grid_size
export make_edge_p_values, make_directed_edge_values

# Specific spatial problem types.
include("spatial_reaction_systems/spatial_ODE_systems.jl")
include("spatial_reaction_systems/discrete_space_jump_systems.jl")

# General spatial modelling utility functions.
include("spatial_reaction_systems/utility.jl")

# Methods for interfacing with from DiscreteSpaceReactionSystem derived problems, integrators, and solutions.
include("spatial_reaction_systems/discrete_space_sim_struct_interfacing.jl")
export spat_getp, spat_setp!, spat_getu, spat_setu!, rebuild_spat_internals!

# Functions for plotting of discrete space simulations (most of the code is in extensions, not here).
include("spatial_reaction_systems/discrete_space_simulation_plotting.jl")
export dspace_plot, dspace_animation, dspace_kymograph

### ReactionSystem Serialisation ###
# Has to be at the end (because it uses records of all metadata declared by Catalyst).
include("reactionsystem_serialisation/serialisation_support.jl")
include("reactionsystem_serialisation/serialise_fields.jl")
include("reactionsystem_serialisation/serialise_reactionsystem.jl")
export save_reactionsystem

end # module
