using DiffEqBiological
using Base.Test

tic()
@time @testset "Gillespie Tests" begin include("gillespie.jl") end
toc()
