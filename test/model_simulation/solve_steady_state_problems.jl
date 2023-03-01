### Fetch Packages and Reaction Networks ###

# Fetch packages.
using Catalyst, OrdinaryDiffEq, Random, SteadyStateDiffEq, Test

# Sets rnd number.
using StableRNGs
rng = StableRNG(12345)

# Fetch test networks.
include("../test_networks.jl")

### Compares Solution to Known Steady States ###

# Simple network.
let
    steady_state_network_1 = @reaction_network begin
        (k1, k2), ∅ ↔ X1
        (k3, k4), ∅ ↔ 3X2
        (k5, k6), ∅ ↔ X3 + X4
    end

    for factor in [1e-1, 1e0, 1e1]
        u0 = factor * rand(rng, length(states(steady_state_network_1)))
        u0[4] = u0[3]
        p = 0.01 .+ factor * rand(rng, length(parameters(steady_state_network_1)))
        prob = SteadyStateProblem(steady_state_network_1, u0, p)
        sol = solve(prob, SSRootfind()).u
        (minimum(sol[1:1]) > 1e-2) && (@test abs.(sol[1] - p[1] / p[2]) < 1e-8)
        (minimum(sol[2:2]) > 1e-2) &&
            (@test abs.(sol[2]^3 / factorial(3) - p[3] / p[4]) < 1e-8)
        (minimum(sol[3:4]) > 1e-2) && (@test abs.(sol[3] * sol[4] - p[5] / p[6]) < 1e-8)
    end
end

# These are disabled due to problem in SteadyStateProblem solution exactness. We do not recommend this method for finding steady states.

# steady_state_network_2 = @reaction_network begin
#     v / 10 + hill(X, v, K, n), ∅ → X
#     d, X → ∅
# end

# for factor in [1e-1, 1e1, 1e1], repeat in 1:3
#   u0_small = factor * rand(rng, length(get_states(steady_state_network_2))) / 100
#     u0_large = factor * rand(rng, length(get_states(steady_state_network_2))) * 100
#     p = factor * rand(rng, length(get_ps(steady_state_network_2)))
#     p[3] = round(p[3]) + 1
#     sol1 = solve(SteadyStateProblem(steady_state_network_2, u0_small, p), SSRootfind()).u[1]
#     sol2 = solve(SteadyStateProblem(steady_state_network_2, u0_large, p), SSRootfind()).u[1]
#     diff1 = abs(p[1] / 10 + p[1] * (sol1^p[3]) / (sol1^p[3] + p[2]^p[3]) - p[4] * sol1)
#     diff2 = abs(p[1] / 10 + p[1] * (sol2^p[3]) / (sol2^p[3] + p[2]^p[3]) - p[4] * sol2)
#     @test (diff1 < 1e-8) || (diff2 < 1e-8)
# end

### For a couple of networks, test that the steady state solution is identical to the long term ODE solution. ###

# steady_state_test_networks = [
#     reaction_networks_standard[8],
#     reaction_networks_standard[10],
#     reaction_networks_weird[1],
# ]
# for network in steady_state_test_networks, factor in [1e-1, 1e0, 1e1]
#     u0 = factor * rand(rng, length(get_states(network)))
#     p = factor * rand(rng, length(get_ps(network)))
#     sol_ode = solve(ODEProblem(network, u0, (0.0, 1000000), p), Rosenbrock23())
#     sol_ss = solve(SteadyStateProblem(network, u0, p), SSRootfind())
# 
#     @test all(abs.(sol_ode.u[end] .- sol_ss.u) .< 1e-4)
# end

### No parameter test ###

# no_param_network = @reaction_network begin (0.6, 3.2), ∅ ↔ X end
# for factor in [1e0, 1e1, 1e2]
#     u0 = factor * rand(rng, length(get_states(no_param_network)))
#     sol_ss = solve(SteadyStateProblem(no_param_network, u0), SSRootfind(), abstol = 1e-11)
#     @test abs.(sol_ss.u[1] - 0.6 / 3.2) < 1e-8
# end
