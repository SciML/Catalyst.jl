### Fetch required packages and reaction networks ###
using Catalyst, OrdinaryDiffEq, Random, SteadyStateDiffEq, Test
include("test_networks.jl")


### Compares to netowork with know steady state ###
steady_state_network_1 = @reaction_network begin
    (k1,k2), ∅ ↔ X1
    (k3,k4), ∅ ↔ 3X2
    (k5,k6), ∅ ↔ X3 + X4
end k1 k2 k3 k4 k5 k6

for factor in [1e-1, 1e0, 1e1], repeat = 1:3
    u0 = factor*rand(length(steady_state_network_1.states)); u0[4] = u0[3];
    p = 0.01 .+ factor*rand(length(steady_state_network_1.ps))
    prob = SteadyStateProblem(steady_state_network_1,u0,p)
    sol = solve(prob,SSRootfind()).u
    (minimum(sol[1:1]) > 1e-2) && (@test abs.(sol[1] - p[1]/p[2]) < 0.01)
    (minimum(sol[2:2]) > 1e-2) && (@test abs.(sol[2]^3/factorial(3) - p[3]/p[4]) < 0.01)
    (minimum(sol[3:4]) > 1e-2) && (@test abs.(sol[3]*sol[4] - p[5]/p[6]) < 0.01)
end

steady_state_network_2 = @reaction_network begin
    v/10+hill(X,v,K,n), ∅ → X
    d, X → ∅
end v K n d

for factor in [1e-1, 1e1, 1e1], repeat = 1:3
    u0_small = factor*rand(length(steady_state_network_2.states))/100
    u0_large = factor*rand(length(steady_state_network_2.states))*100
    p = factor*rand(length(steady_state_network_2.ps))
    p[3] = round(p[3])+1
    sol1 = solve(SteadyStateProblem(steady_state_network_2,u0_small,p),SSRootfind()).u[1]
    sol2 = solve(SteadyStateProblem(steady_state_network_2,u0_large,p),SSRootfind()).u[1]
    diff1 = abs(p[1]/10 + p[1]*(sol1^p[3])/(sol1^p[3]+p[2]^p[3]) - p[4]*sol1)
    diff2 = abs(p[1]/10 + p[1]*(sol2^p[3])/(sol2^p[3]+p[2]^p[3]) - p[4]*sol2)
    @test (diff1 < 1e-8) || (diff2 < 1e-8)
end


### For a couple of networks, test that the steady state solution is identical to the long term ODE solution. ###
steady_state_test_networks = [reaction_networks_standard[8], reaction_networks_standard[10], reaction_networks_weird[1]]
for network in steady_state_test_networks, factor in [1e-1, 1e0, 1e1]
    u0 = factor*rand(length(network.states))
    p = factor*rand(length(network.ps))
    sol_ode = solve(ODEProblem(network,u0,(0.,1000000),p),Rosenbrock23())
    sol_ss = solve(SteadyStateProblem(network,u0,p),SSRootfind())

    @test all(abs.(sol_ode.u[end] .- sol_ss.u) .< 1e-4)
end
