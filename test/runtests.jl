using DiffEqBiological
using Base.Test

tic()
#@time @testset "Gillespie Tests" begin include("gillespie.jl") end
#@time @testset "Variable Rate Reaction Tests" begin include("variable_rate_reactions.jl") end
#@time @testset "@reaction_network Tests" begin include("reaction_network.jl") end
#@time @testset "Higher order reaction Tests" begin include("higher_order_reactions.jl") end

@time @testset "Test the model creation macro" begin include("make_model.jl") end
@time @testset "Tests the solver methods when run on the reaction networks" begin include("solver_test.jl") end
@time @testset "Tests the jump reaction simulations when run on the reaction networks Tests" begin include("jump_reaction_test.jl") end
@time @testset "Tests some of the additional functionalities" begin include("func_test.jl") end
toc()
