module CatalystHomotopyContinuationExtension

# Fetch packages.
using Catalyst: Catalyst, ReactionSystem, conservationlaw_constants, conservedequations,
    isautonomous, ss_ode_model
using ModelingToolkitBase: ModelingToolkitBase, complete, equations, get_iv, parameters,
    unknowns, @parameters
using SciMLBase: NonlinearProblem
using SymbolicUtils: Rewriters, simplify_fractions, substitute, unwrap, unwrap_const
using Symbolics: Symbolics, SymbolicT, iscall, wrap, @variables
using TermInterface: arguments, maketerm, metadata, operation, sorted_arguments
import DiffEqBase
import DynamicPolynomials
import ModelingToolkitBase as MT
import HomotopyContinuation as HC
using DataStructures: OrderedSet

# Workaround for Julia 1.10 precompilation bug (related to ModelingToolkit.jl#4211).
# On Julia 1.10, when Catalyst and HomotopyContinuation are loaded simultaneously,
# the precompiled CodeInstance for collect_vars! can be incorrectly optimized to a no-op.
# This forces fresh compilation after both packages are loaded.
function __init__()
    if VERSION < v"1.11"
        _warmup_collect_vars()
    end
    return nothing
end

function _warmup_collect_vars()
    @variables _dummy_t
    @parameters _dummy_p
    _us = OrderedSet{SymbolicT}()
    _ps = OrderedSet{SymbolicT}()
    expr = unwrap(_dummy_p + _dummy_t)
    MT.collect_vars!(_us, _ps, expr, unwrap(_dummy_t))
    return nothing
end

# Creates and exports hc_steady_states function.
include("CatalystHomotopyContinuationExtension/homotopy_continuation_extension.jl")

end
