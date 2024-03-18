### Not currently run by runtests

### Fetch Packages and Reaction Networks ###
using Catalyst, JumpProcesses, OrdinaryDiffEq, Plots, Random, Test
using ModelingToolkit: get_unknowns, get_ps

# Sets rnd number.
using StableRNGs
rng = StableRNG(12345)

# Fetch test networks.
include("../test_networks.jl")

### Run Tests ###

# Tests the plot() function on a few basic simulation solutions, checks that there are no errors.
let
    plotting_test_networks = [
        reaction_networks_standard[6],
        reaction_networks_constraint[6],
        reaction_networks_real[3],
    ]
    for network in plotting_test_networks, factor in [1e-1, 1e0, 1e1]
        u0 = factor * rand(rng, length(get_unknowns(network)))
        p = factor * rand(rng, length(get_ps(network)))
        sol = solve(ODEProblem(network, u0, (0.0, 1), p), Rosenbrock23())
        plot(sol)

        u0 = rand(rng, 1:Int64(factor * 100), length(get_unknowns(network)))
        sol = solve(JumpProblem(network, DiscreteProblem(network, u0, (0.0, 1), p),
                                Direct()),
                    SSAStepper())
        plot(sol)
    end
end
