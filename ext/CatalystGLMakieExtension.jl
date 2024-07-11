module CatalystGLMakieExtension

# Fetch packages.
using Catalyst, GLMakie

# Creates and exports hc_steady_states function.
include("CatalystGLMakieExtension/glmakie_extension_spatial_modelling.jl")

end
