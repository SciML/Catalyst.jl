#! format: off

### Fetch Packages and Reaction Networks ###

# Fetch packages.
using Catalyst, OrdinaryDiffEq, JumpProcesses, Random, Test
using ModelingToolkit: get_states, get_ps
t = default_t()

# Sets rnd number.
using StableRNGs
rng = StableRNG(12345)

# Fetch test networks.
include("../test_networks.jl")

### Run Tests ###

# Tests various ways to input u0 and p for various functions.
let
    test_network = reaction_networks_standard[7]
    test_osys = complete(convert(ODESystem, test_network))
    @parameters p1 p2 p3 k1 k2 k3 v1 K1 d1 d2 d3 d4 d5
    @species X1(t) X2(t) X3(t) X4(t) X5(t) X6(t) X(t)

    for factor = [1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3]
        u0_1 = factor*rand(rng,length(unknowns(test_network)))
        u0_2 = [X1=>u0_1[1], X2=>u0_1[2], X3=>u0_1[3], X4=>u0_1[4], X5=>u0_1[5]]
        u0_3 = [X2=>u0_1[2], X5=>u0_1[5], X4=>u0_1[4], X3=>u0_1[3], X1=>u0_1[1]]
        p_1 = 0.01 .+ factor*rand(rng,length(parameters(test_network)))
        p_2 = [p1=>p_1[1], p2=>p_1[2], p3=>p_1[3], k1=>p_1[4], k2=>p_1[5], k3=>p_1[6],
            v1=>p_1[7], K1=>p_1[8], d1=>p_1[9], d2=>p_1[10], d3=>p_1[11], d4=>p_1[12],
            d5=>p_1[13]]
        p_3 = [k2=>p_1[5], k3=>p_1[6], v1=>p_1[7], d5=>p_1[13], p2=>p_1[2], p1=>p_1[1],
            d2=>p_1[10], K1=>p_1[8], d1=>p_1[9], d4=>p_1[12], d3=>p_1[11], p3=>p_1[3],
            k1=>p_1[4]]

        sols = []
        push!(sols, solve(ODEProblem(test_osys,u0_1, (0.,10.), p_1), Rosenbrock23()))
        push!(sols, solve(ODEProblem(test_osys,u0_1, (0.,10.), p_2), Rosenbrock23()))
        push!(sols, solve(ODEProblem(test_osys,u0_1, (0.,10.), p_3), Rosenbrock23()))
        push!(sols, solve(ODEProblem(test_osys,u0_2, (0.,10.), p_1), Rosenbrock23()))
        push!(sols, solve(ODEProblem(test_osys,u0_2, (0.,10.), p_2), Rosenbrock23()))
        push!(sols, solve(ODEProblem(test_osys,u0_2, (0.,10.), p_3), Rosenbrock23()))
        push!(sols, solve(ODEProblem(test_osys,u0_3, (0.,10.), p_1), Rosenbrock23()))
        push!(sols, solve(ODEProblem(test_osys,u0_3, (0.,10.), p_2), Rosenbrock23()))
        push!(sols, solve(ODEProblem(test_osys,u0_3, (0.,10.), p_3), Rosenbrock23()))

        ends = map(sol -> sol.u[end],sols)
        for i in 1:length(u0_1)
            @test abs(maximum(getindex.(ends,1)) - minimum(first.(getindex.(ends,1)))) < 1e-5
        end

        u0_1 = rand(rng, 1:Int64(factor*100), length(unknowns(test_network)))
        u0_2 = [X1=>u0_1[1], X2=>u0_1[2], X3=>u0_1[3], X4=>u0_1[4], X5=>u0_1[5]]
        u0_3 = [X2=>u0_1[2], X5=>u0_1[5], X4=>u0_1[4], X3=>u0_1[3], X1=>u0_1[1]]

        discrete_probs = []
        push!(discrete_probs, DiscreteProblem(test_network, u0_1, (0.,1.), p_1))
        push!(discrete_probs, DiscreteProblem(test_network, u0_1, (0.,1.), p_2))
        push!(discrete_probs, DiscreteProblem(test_network, u0_1, (0.,1.), p_3))
        push!(discrete_probs, DiscreteProblem(test_network, u0_2, (0.,1.), p_1))
        push!(discrete_probs, DiscreteProblem(test_network, u0_2, (0.,1.), p_2))
        push!(discrete_probs, DiscreteProblem(test_network, u0_2, (0.,1.), p_3))
        push!(discrete_probs, DiscreteProblem(test_network, u0_3, (0.,1.), p_1))
        push!(discrete_probs, DiscreteProblem(test_network, u0_3, (0.,1.), p_2))
        push!(discrete_probs, DiscreteProblem(test_network, u0_3, (0.,1.), p_3))

        for i in 2:9
            @test discrete_probs[1].p == discrete_probs[i].p
            @test discrete_probs[1].u0 == discrete_probs[i].u0
        end
    end
end

# Tests solving for various inputs types across various problem types.
let 
    model = complete(@reaction_network begin
        (kp,kd), 0 <--> X
        (k1,k2), X <--> Y
    end)
    @unpack X, Y, kp, kd, k1, k2 = model
    
    tspan = (0.0, 10.0)

    u0_vals_1 = [X => 1, Y => 0]
    u0_vals_2 = [model.X => 1, model.Y => 0]
    u0_vals_3 = [:X => 1, :Y => 0]
    p_vals_1 = [kp => 1.0, kd => 0.2, k1 => 1.0, k2 => 2.0]
    p_vals_2 = [model.kp => 1.0, model.kd => 0.2, model.k1 => 1.0, model.k2 => 2.0]
    p_vals_3 = [:kp => 1.0, :kd => 0.2, :k1 => 1.0, :k2 => 2.0]

    oprob_1 = ODEProblem(model, u0_vals_1, tspan, p_vals_1)
    oprob_2 = ODEProblem(model, u0_vals_2, tspan, p_vals_2)
    oprob_3 = ODEProblem(model, u0_vals_3, tspan, p_vals_3)
    sprob_1 = SDEProblem(model,u0_vals_1, tspan, p_vals_1)
    sprob_2 = SDEProblem(model,u0_vals_2, tspan, p_vals_2)
    sprob_3 = SDEProblem(model,u0_vals_3, tspan, p_vals_3)
    dprob_1 = DiscreteProblem(model, u0_vals_1, tspan, p_vals_1)
    dprob_2 = DiscreteProblem(model, u0_vals_2, tspan, p_vals_2)
    dprob_3 = DiscreteProblem(model, u0_vals_3, tspan, p_vals_3)
    jprob_1 = JumpProblem(model, dprob_1, Direct())
    jprob_2 = JumpProblem(model, dprob_2, Direct())
    jprob_3 = JumpProblem(model, dprob_3, Direct())
    nprob_1 = NonlinearProblem(model, u0_vals_1, p_vals_1)
    nprob_2 = NonlinearProblem(model, u0_vals_2, p_vals_2)
    nprob_3 = NonlinearProblem(model, u0_vals_3, p_vals_3)
    
    @test solve(oprob_1, Tsit5()) == solve(oprob_2, Tsit5()) == solve(oprob_3, Tsit5())
    @test solve(sprob_1, ImplicitEM(); seed=1234) == solve(sprob_2, ImplicitEM(); seed=1234) == solve(sprob_3, ImplicitEM(); seed=1234)
    @test solve(jprob_1, SSAStepper(); seed=1234) == solve(jprob_2, SSAStepper(); seed=1234) == solve(jprob_3, SSAStepper(); seed=1234)
    @test solve(nprob_1, NewtonRaphson()) == solve(nprob_2, NewtonRaphson()) == solve(nprob_3, NewtonRaphson())
end