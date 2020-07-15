### Not currently run by runtests 
### Fetch required packages and reaction networks ###
using DiffEqBiological, DiffEqJump, OrdinaryDiffEq, Plots, Random, Test
include("test_networks.jl")


### Tests the plot() function on a few basic simulation solutions, checks that there are no errors ###

plotting_test_networks = [reaction_networks_standard[6], reaction_networks_constraint[6],reaction_networks_real[3]]
for network in plotting_test_networks, factor in [1e-1, 1e0, 1e1]
    u0 = factor*rand(length(network.states))
    p = factor*rand(length(network.ps))
    sol = solve(ODEProblem(network,u0,(0.,1),p),Rosenbrock23())
    plot(sol)

    u0 = rand(1:Int64(factor*100),length(network.states))
    sol = solve(JumpProblem(network,DiscreteProblem(network,u0,(0.,1),p),Direct()),SSAStepper())
    plot(sol)
end
