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

import Base: ==, hash, size, getindex, setindex, isless, Sort.defalg, length, show
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
"""
    default_time_deriv()

Return Catalyst's default time derivative operator.

This is the derivative with respect to [`default_t`](@ref), and is useful when
programmatically constructing equations involving time derivatives.

# Examples
```julia
t = default_t()
D = default_time_deriv()
@species X(t)
eq = D(X) ~ -X
```
"""
function default_time_deriv()
    return ModelingToolkitBase.D_nounits
end

"""
    default_t()

Return Catalyst's default independent time variable.

Use this when programmatically declaring species or constructing
[`ReactionSystem`](@ref)s so that the independent variable matches Catalyst's
DSL defaults.

# Examples
```julia
t = default_t()
@species X(t)
```
"""
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
"""
    hc_steady_states(rs::ReactionSystem, ps; filter_negative = true, neg_thres = -1e-15, u0 = [], kwargs...)

Find steady states of a [`ReactionSystem`](@ref) using HomotopyContinuation.jl.

This extension function is available after loading HomotopyContinuation.jl.
`ps` supplies parameter values, `u0` supplies initial conditions needed for
systems with conservation laws, and extra keyword arguments are passed to
HomotopyContinuation's solver.

# Examples
```julia
using Catalyst, HomotopyContinuation

rs = @reaction_network begin
    k1, Y --> 2X
    k2, 2X --> X + Y
    k3, X + Y --> Y
    k4, X --> 0
end
hc_steady_states(rs, [:k1 => 8.0, :k2 => 2.0, :k3 => 1.0, :k4 => 1.5])
```
"""
function hc_steady_states end
export hc_steady_states

# StructuralIdentifiability
"""
    make_si_ode(rs::ReactionSystem; measured_quantities = [], known_p = [],
        ignore_no_measured_warn = false, remove_conserved = true)

Convert a [`ReactionSystem`](@ref) to the ODE representation used by
StructuralIdentifiability.jl.

This extension function is available after loading StructuralIdentifiability.jl.
`measured_quantities` lists measured species or equations, `known_p` lists
parameters treated as known, and `remove_conserved` controls conservation-law
elimination before conversion.

# Examples
```julia
using Catalyst, StructuralIdentifiability

rs = @reaction_network begin
    (p, d), 0 <--> X
end
make_si_ode(rs; measured_quantities = [:X], known_p = [:p])
```
"""
function make_si_ode end
export make_si_ode

# GraphMakie: functionality for plotting species-reaction graphs and complexes
"""
    plot_network(rn::ReactionSystem; kwargs...)

Plot the species-reaction graph of `rn`.

This extension function is available after loading GraphMakie.jl and
NetworkLayout.jl. Keyword arguments are forwarded to the GraphMakie plotting
recipe.

# Examples
```julia
using Catalyst, GraphMakie, CairoMakie

rn = @reaction_network begin
    k, A --> B
end
plot_network(rn)
```
"""
function plot_network end

"""
    plot_complexes(rn::ReactionSystem; show_rate_labels = false, kwargs...)

Plot the reaction-complex graph of `rn`.

This extension function is available after loading GraphMakie.jl and
NetworkLayout.jl. Set `show_rate_labels = true` to label graph edges by their
reaction rates.

# Examples
```julia
using Catalyst, GraphMakie, CairoMakie

rn = @reaction_network begin
    k, A --> B
end
plot_complexes(rn)
```
"""
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
export CartesianGrid, CartesianGridRej # (Implemented in JumpProcesses)
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
