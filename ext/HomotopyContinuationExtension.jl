module HomotopyContinuationExtension

# Fetch packages.
using Catalyst
import ModelingToolkit as MT
import HomotopyContinuation as HC
import Symbolics: unwrap, wrap, Rewriters, symtype, issym, istree

# Creates and exports hc_steady_states function.
include("CatalystHomotopyContinuationExtension/homotopy_continuation_extension.jl")

end