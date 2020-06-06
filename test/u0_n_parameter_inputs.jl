### Fetch required packages and reaction networks ###
using DiffEqBiological, DiffEqJump, OrdinaryDiffEq, Random, Statistics, Test
include("test_networks.jl")

### Tests various ways to input u0 and p for various functions ###

# Tests for reaction_networks_standard[7]
test_network = reaction_networks_standard[7]
@variables X1 X2 X3 X4 X5 X6 X
@parameters p1 p2 p3 k1 k2 k3 v1 K1 d1 d2 d3 d4 d5

for factor in [1e0, 1e1]
    u0_1 = factor*rand(length(test_network.states))
    u0_2 = [X1=>u0_1[1], X2=>u0_1[2], X3=>u0_1[3], X4=>u0_1[4], X5=>u0_1[5]]
    u0_3 = [X2=>u0_1[2], X5=>u0_1[5], X4=>u0_1[4], X3=>u0_1[3], X1=>u0_1[1]]
    p_1 = 0.01 .+ factor*rand(length(test_network.ps))
    p_2 = [p1=>p_1[1], p2=>p_1[2], p3=>p_1[3], k1=>p_1[4], k2=>p_1[5], k3=>p_1[6], v1=>p_1[7], K1=>p_1[8], d1=>p_1[9], d2=>p_1[10], d3=>p_1[11], d4=>p_1[12], d5=>p_1[13]]
    p_3 = [k2=>p_1[5], k3=>p_1[6], v1=>p_1[7], d5=>p_1[13], p2=>p_1[2], p1=>p_1[1], d2=>p_1[10], K1=>p_1[8], d1=>p_1[9], d4=>p_1[12], d3=>p_1[11], p3=>p_1[3], k1=>p_1[4]]

    sols = []
    push!(sols,solve(ODEProblem(test_network,u0_1,(0.,10.),p_1),Rosenbrock23()))
    push!(sols,solve(ODEProblem(test_network,u0_1,(0.,10.),p_2),Rosenbrock23()))
    push!(sols,solve(ODEProblem(test_network,u0_1,(0.,10.),p_3),Rosenbrock23()))
    push!(sols,solve(ODEProblem(test_network,u0_2,(0.,10.),p_1),Rosenbrock23()))
    push!(sols,solve(ODEProblem(test_network,u0_2,(0.,10.),p_2),Rosenbrock23()))
    push!(sols,solve(ODEProblem(test_network,u0_2,(0.,10.),p_3),Rosenbrock23()))
    push!(sols,solve(ODEProblem(test_network,u0_3,(0.,10.),p_1),Rosenbrock23()))
    push!(sols,solve(ODEProblem(test_network,u0_3,(0.,10.),p_2),Rosenbrock23()))
    push!(sols,solve(ODEProblem(test_network,u0_3,(0.,10.),p_3),Rosenbrock23()))

    ends = map(sol -> sol.u[end],sols)
    for i in 1:length(u0_1)
        @test abs(maximum(getindex.(ends,1))-minimum(first.(getindex.(ends,1)))) < 1e-5
    end

    u0_1 = rand(1:Int64(factor*100),length(test_network.states))
    u0_2 = [X1=>u0_1[1], X2=>u0_1[2], X3=>u0_1[3], X4=>u0_1[4], X5=>u0_1[5]]
    u0_3 = [X2=>u0_1[2], X5=>u0_1[5], X4=>u0_1[4], X3=>u0_1[3], X1=>u0_1[1]]

    sols = []
    push!(sols,solve(JumpProblem(test_network,DiscreteProblem(test_network,u0_1,(0.,10000000.),p_1),Direct()),SSAStepper()))
    push!(sols,solve(JumpProblem(test_network,DiscreteProblem(test_network,u0_1,(0.,10000000.),p_2),Direct()),SSAStepper()))
    push!(sols,solve(JumpProblem(test_network,DiscreteProblem(test_network,u0_1,(0.,10000000.),p_3),Direct()),SSAStepper()))
    push!(sols,solve(JumpProblem(test_network,DiscreteProblem(test_network,u0_2,(0.,10000000.),p_1),Direct()),SSAStepper()))
    push!(sols,solve(JumpProblem(test_network,DiscreteProblem(test_network,u0_2,(0.,10000000.),p_2),Direct()),SSAStepper()))
    push!(sols,solve(JumpProblem(test_network,DiscreteProblem(test_network,u0_2,(0.,10000000.),p_3),Direct()),SSAStepper()))
    push!(sols,solve(JumpProblem(test_network,DiscreteProblem(test_network,u0_3,(0.,10000000.),p_1),Direct()),SSAStepper()))
    push!(sols,solve(JumpProblem(test_network,DiscreteProblem(test_network,u0_3,(0.,10000000.),p_2),Direct()),SSAStepper()))
    push!(sols,solve(JumpProblem(test_network,DiscreteProblem(test_network,u0_3,(0.,10000000.),p_3),Direct()),SSAStepper()))

    for i in 1:length(u0_1)
        means = map(sol -> mean(getindex.(sol.u,i)), sols)
        stds = map(sol -> std(getindex.(sol.u,i)), sols)
        (minimum(means)>0.001) && @test maximum(means)/minimum(means) < 1.5
        (minimum(stds)>0.001) && @test maximum(stds)/minimum(stds) < 1.5
    end
end
