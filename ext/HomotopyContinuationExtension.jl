module HomotopyContinuationExtension

# Fetch packages.
using Catalyst
import ModelingToolkit as MT
import HomotopyContinuation as HC

# Creates and exports hc_steady_states function.
include("HomotopyContinuationExtension/homotopy_continuation_extension.jl")

end