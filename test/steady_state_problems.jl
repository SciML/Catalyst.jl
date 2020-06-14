### Fetch required packages and reaction networks ###
using DiffEqBiological, OrdinaryDiffEq, Random, SteadyStateDiffEq, Test
include("test_networks.jl")


### Compares to netowork with know steady state ###
steady_state_network_1 = @reaction_network begin
    (k1,k2), X1 ↔ X2
    (k3,k4), X3 + X4 ↔ X5
    (k5,k6), 2X6 ↔ 3X7
    (k7,k8), ∅ ↔ X8
end k1 k2 k3 k4 k5 k6 k7 k8

for factor in [1e-1, 1e0, 1e1], repeat = 1:5
    u0 = factor*rand(length(steady_state_network_1.states))
    p = 0.01 .+ factor*rand(length(steady_state_network_1.ps))
    prob = SteadyStateProblem(steady_state_network_1,u0,p)
    sol = solve(prob,SSRootfind()).u
    (sol[1] > 1e-7) && (@test abs.(sol[1]/sol[2] - p[2]/p[1]) < 1e-10)
    (sol[3] > 1e-7) && (@test abs.(sol[3]*sol[4]/sol[5] - p[4]/p[3]) < 1e-6)
    (sol[6] > 1e-7) && (@test abs.((sol[6]^2/factorial(2))/(sol[7]^3/factorial(3))- p[6]/p[5]) < 1e-2)
    (sol[8] > 1e-7) && (@test abs.(sol[8] - p[7]/p[8]) < 1e-10)
end

steady_state_network_2 = @reaction_network begin
    v/10+hill(X,v,K,n), ∅ → X
    d, X → ∅
end v K n d

for factor in [1e-1, 1e1, 1e1], repeat = 1:3
    u0 = factor*rand(length(steady_state_network_2.states))
    p = factor*rand(length(steady_state_network_2.ps))
    p[3] = round(p[3])+1
    prob = SteadyStateProblem(steady_state_network_2,u0,p)
    sol = solve(prob,SSRootfind()).u[1]
    @test_broken abs(p[1]/10 + p[1]*(sol^p[3])/(sol^p[3]+p[2]^p[3]) - p[4]*sol) < 1e-8
end


### For a couple of networks, test that the steady state solution is identical to the long term ODE solution. ###
steady_state_test_networks = [reaction_networks_standard[8], reaction_networks_real[1], reaction_networks_weird[1]]
for network in steady_state_test_networks, factor in [1e-1, 1e0, 1e1]
    u0 = factor*rand(length(network.states))
    p = factor*rand(length(network.ps))
    sol_ode = solve(ODEProblem(network,u0,(0.,1000000),p),Rosenbrock23())
    sol_ss = solve(SteadyStateProblem(network,u0,p),SSRootfind())

    @test all(abs.(sol_ode.u[end] .- sol_ss.u) .< 1e-4)
end
