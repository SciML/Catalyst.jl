using Test, Statistics, Random, LinearAlgebra, SparseArrays, Plots
using DiffEqBiological, DiffEqBase
using OrdinaryDiffEq, StochasticDiffEq, DiffEqJump, SteadyStateDiffEq

# fetches a set of networks which can be used for tests.
include("test_networks.jl")

# full macro tests
@time begin
  @time @testset "Model Macro" begin include("make_model_test.jl") end
  @time @testset "Gillespie Tests" begin include("gillespie.jl") end
  @time @testset "Test Solvers" begin include("solver_test.jl") end
  @time @testset "Higher Order" begin include("higher_order_reactions.jl") end
  @time @testset "Additional Functions" begin include("func_test.jl") end
  @time @testset "Steady State Solver" begin include("steady_state.jl") end
  @time @testset "Mass Action Jumps" begin include("mass_act_jump_tests.jl") end
  @time @testset "Equilibrate (1)" begin include("equilibrate_test_1.jl") end
  @time @testset "Equilibrate (2)" begin include("equilibrate_test_2.jl") end
  @time @testset "Equilibrate (3)" begin include("equilibrate_test_3.jl") end
  @time @testset "Equilibrate (4)" begin include("equilibrate_test_4.jl") end
  @time @testset "Other Tests" begin include("misc_tests.jl") end
  @time @testset "Network query tests" begin include("networkquery_test.jl") end
end

# min macro tests
@time begin
  @time @testset "Model Macro (Min)" begin include("make_model_test_min.jl") end
  @time @testset "Gillespie Tests (Min)" begin include("gillespie_min.jl") end
  @time @testset "Test Solvers (Min)" begin include("solver_test_min.jl") end
  @time @testset "Higher Order (Min)" begin include("higher_order_reactions_min.jl") end
  @time @testset "Additional Functions (Min)" begin include("func_test_min.jl") end
  @time @testset "Steady State Solver (Min)" begin include("steady_state_min.jl") end
  @time @testset "Equilibrate (Min)" begin include("equilibrate_test_min.jl") end
  @time @testset "Mass Action Jumps (Min)" begin include("mass_act_jump_tests_min.jl") end
end

# tests that handle both macros
@time begin
  @time @testset "Discrete Problem" begin include("discreteproblem_test.jl") end
  @time @testset "Add Reactions API" begin include("addreactions_test.jl") end
end
