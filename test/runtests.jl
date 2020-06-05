### Fecth the require packages ###
using DiffEqBiological, ModelingToolkit
using OrdinaryDiffEq, StochasticDiffEq, DiffEqJump
using Test, SafeTestsets
using LinearAlgebra, Random, Statistics
using UnPack, SparseArrays


### Decalre constants realted to groups.
const GROUP = get(ENV, "GROUP", "All")
const is_APPVEYOR = Sys.iswindows() && haskey(ENV,"APPVEYOR")
const is_TRAVIS = haskey(ENV,"TRAVIS")



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
if GROUP == "All" || GROUP == "ODE Problems"
  @time @safetestset "" begin include(".jl") end
  @time @safetestset "" begin include(".jl") end
end

# Tests related to solving Stochastic Differential Equations.
if GROUP == "All" || GROUP == "SDE Problems"
  @time @safetestset "" begin include(".jl") end
  @time @safetestset "" begin include(".jl") end
end

# Tests related to solvingJump Systems.
if GROUP == "All" || GROUP == "Jump Problems"
  @time @safetestset "" begin include(".jl") end
  @time @safetestset "" begin include(".jl") end
end

# Miscellaneous tests.
if GROUP == "All" || GROUP == "Miscellaneous"
  @time @safetestset "" begin include(".jl") end
  @time @safetestset "" begin include(".jl") end
end



if GROUP == "All" || GROUP == "Core"
  @time @safetestset "Model Macro" begin include("make_model_test.jl") end
  @time @safetestset "Gillespie Tests" begin include("gillespie.jl") end
  @time @safetestset "Test Solvers" begin include("solver_test.jl") end
  @time @safetestset "Higher Order" begin include("higher_order_reactions.jl") end
  @time @safetestset "Additional Functions" begin include("func_test.jl") end
  @time @safetestset "Steady State Solver" begin include("steady_state.jl") end
  @time @safetestset "Mass Action Jumps" begin include("mass_act_jump_tests.jl") end
  @time @safetestset "Other Tests" begin include("misc_tests.jl") end
  @time @safetestset "Network query tests" begin include("networkquery_test.jl") end
end

if !is_APPVEYOR && (GROUP == "All" || GROUP == "EquibrilateI")
  @time @safetestset "Equilibrate (1)" begin include("equilibrate_test_1.jl") end
  @time @safetestset "Equilibrate (2)" begin include("equilibrate_test_2.jl") end
end

if !is_APPVEYOR && (GROUP == "All" || GROUP == "EquibrilateII")
  @time @safetestset "Equilibrate (3)" begin include("equilibrate_test_3.jl") end
  @time @safetestset "Equilibrate (4)" begin include("equilibrate_test_4.jl") end
end

# min macro tests
if !is_APPVEYOR && (GROUP == "All" || GROUP == "Min")
  @time @safetestset "Model Macro (Min)" begin include("make_model_test_min.jl") end
  @time @safetestset "Gillespie Tests (Min)" begin include("gillespie_min.jl") end
  @time @safetestset "Test Solvers (Min)" begin include("solver_test_min.jl") end
  @time @safetestset "Higher Order (Min)" begin include("higher_order_reactions_min.jl") end
  @time @safetestset "Additional Functions (Min)" begin include("func_test_min.jl") end
  @time @safetestset "Steady State Solver (Min)" begin include("steady_state_min.jl") end
  @time @safetestset "Equilibrate (Min)" begin include("equilibrate_test_min.jl") end
  @time @safetestset "Mass Action Jumps (Min)" begin include("mass_act_jump_tests_min.jl") end
end

# tests that handle both macros
if GROUP == "All" || GROUP == "Misc"
  @time @safetestset "Discrete Problem" begin include("discreteproblem_test.jl") end
  @time @safetestset "Add Reactions API" begin include("addreactions_test.jl") end
  @time @safetestset "Latexify recipe" begin include("latexify_test.jl") end
end

end # @time
