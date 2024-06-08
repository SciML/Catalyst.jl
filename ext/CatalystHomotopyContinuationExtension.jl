module CatalystHomotopyContinuationExtension

# Fetch packages.
using Catalyst
import DynamicPolynomials
import ModelingToolkit as MT
import HomotopyContinuation as HC
import Setfield: @set
import Symbolics: unwrap, wrap, Rewriters, symtype, issym
using Symbolics: iscall

# Creates and exports hc_steady_states function.
include("CatalystHomotopyContinuationExtension/homotopy_continuation_extension.jl")

end
