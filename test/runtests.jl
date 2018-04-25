using DiffEqBiological
using Base.Test

tic()
@time @testset "Model Macro" begin include("make_model_test.jl") end
@time @testset "Gillespie Tests" begin include("gillespie.jl") end
@time @testset "Test Solvers" begin include("solver_test.jl") end
@time @testset "Higher Order" begin include("higher_order_reactions.jl") end
@time @testset "Additional Functions" begin include("func_test.jl") end
@time @testset "Steady state solver" begin include("steady_state.jl") end
@time @testset "Mass Action Jumps" begin include("mass_act_jump_tests.jl") end
toc()
