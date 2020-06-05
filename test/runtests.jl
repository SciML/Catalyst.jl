### Fecth the require packages ###
using DiffEqBiological, ModelingToolkit
using OrdinaryDiffEq, StochasticDiffEq, DiffEqJump
using Test, SafeTestsets
using UnPack, LinearAlgebra, Random, Statistics
using Latexify, Plots


### Decalre constants realted to groups.
const GROUP = get(ENV, "GROUP", "All")
const is_APPVEYOR = Sys.iswindows() && haskey(ENV,"APPVEYOR")
const is_TRAVIS = haskey(ENV,"TRAVIS")


### Fetched the set of test networks ###
include("test_networks.jl")

### Ruyn the tests ###
@time begin

# Standard tests realting to the model.
if GROUP == "All" || GROUP == "Core"
  @time @safetestset "Model Macro" begin include("make_model.jl") end
  @time @safetestset "Custom Functions" begin include("custom_functions.jl") end
  @time @safetestset "Higher Order" begin include("higher_order_reactions.jl") end
end

# Tests relating to modification an already created model.
if GROUP == "All" || GROUP == "Model Modifcations"
  @time @safetestset "" begin include(".jl") end
  @time @safetestset "" begin include(".jl") end
end

# Tests related to solving Ordinary Differential Equations.
if GROUP == "All" || GROUP == "Deterministic Tests"
  @time @safetestset "" begin include("solve_ODEs.jl") end
  @time @safetestset "" begin include("make_jacobian.jl") end
end

# Tests related to solving Stochastic Differential Equations.
if GROUP == "All" || GROUP == "SDE Problems"
  @time @safetestset "" begin include("solve_SDEs.jl") end
end

# Tests related to solvingJump Systems.
if GROUP == "All" || GROUP == "Jump Problems"
  @time @safetestset "" begin include("latexify.jl") end
  @time @safetestset "" begin include("tests.jl") end
end

# Miscellaneous tests.
if GROUP == "All" || GROUP == "Miscellaneous"
  @time @safetestset "" begin include(".jl") end
  #@time @safetestset "" begin include("latexify.jl") end
end

end # @time
