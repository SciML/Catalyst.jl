using DiffEqBiological
using Base.Test

tic()
@time @testset "Gillespie Tests" begin include("gillespie.jl") end
@time @testset "Variable Rate Reaction Tests" begin include("variable_rate_reactions.jl") end
toc()
