module CatalystHomotopyContinuationExtension

# Fetch packages.
using Catalyst: Catalyst, ReactionSystem, NonlinearSystem, ModelingToolkit,
                isautonomous, get_iv, complete, conservationlaw_constants,
                conservedequations, species, unknowns, parameters, substitute,
                symmap_to_varmap, get_networkproperties, conservationlaw_errorcheck,
                hc_steady_states, Initial, equations, defaults, mm, hill
import DynamicPolynomials
import ModelingToolkit as MT
import HomotopyContinuation as HC
import Setfield: @set
import Symbolics: unwrap, wrap, Rewriters, symtype, issym, maketerm, BasicSymbolic, metadata
using Symbolics: iscall, simplify, simplify_fractions, solve, arguments, operation,
                 sorted_arguments

# Creates and exports hc_steady_states function.
include("CatalystHomotopyContinuationExtension/homotopy_continuation_extension.jl")

end
