### Fetch Packages ###

# Fetch packages.
using Catalyst, NonlinearSolve, OrdinaryDiffEq, SteadyStateDiffEq
using Random, Test

# Sets rnd number.
using StableRNGs
rng = StableRNG(12345)

# Fetch test networks and functions.
include("../test_functions.jl")

### Run Tests ###

# Creates a simple problem and find steady states just different approaches. Compares to analytic solution.
let 
    # Model with easily computable steady states.
    steady_state_network_1 = @reaction_network begin
        (k1, k2), ∅ ↔ X1
        (k3, k4), ∅ ↔ 3X2
        (k5, k6*X1), ∅ ↔ X3
    end
    
    # Creates NonlinearProblem.
    u0 = rnd_u0(steady_state_network_1, rng)
    ps = rnd_ps(steady_state_network_1, rng)
    nlprob = NonlinearProblem(steady_state_network_1, u0, ps)
    
    # Solves it using standard algorithm and simulation based algorithm.
    sol1 = solve(nlprob; abstol=1e-12, reltol=1e-12)
    sol2 = solve(nlprob, DynamicSS(Rosenbrock23()); abstol=1e-12, reltol=1e-12)
    
    # Tests solutions are correct.
    @test sol1[:X1] ≈ nlprob.ps[:k1] / nlprob.ps[:k2] atol=1e-10
    @test sol1[:X2]^3 / factorial(3) ≈ nlprob.ps[:k3] / nlprob.ps[:k4] atol=1e-10
    @test sol1[:X3] ≈ (nlprob.ps[:k5] * nlprob.ps[:k2]) / (nlprob.ps[:k6] * nlprob.ps[:k1]) atol=1e-10
    @test sol2[:X1] ≈ nlprob.ps[:k1] / nlprob.ps[:k2] atol=1e-10
    @test sol2[:X2]^3 / factorial(3) ≈ nlprob.ps[:k3] / nlprob.ps[:k4] atol=1e-10
    @test sol2[:X3] ≈ (nlprob.ps[:k5] * nlprob.ps[:k2]) / (nlprob.ps[:k6] * nlprob.ps[:k1]) atol=1e-10
end

# Creates a system with multiple steady states.
# Checks that corresponding ODEFunction return 0.0 in both cases. Checks for manually computed function as well.
let 
    # Creates steady state network.
    steady_state_network_2 = @reaction_network begin
        v / 10 + hill(X, v, K, n), ∅ → X
        d, X → ∅
    end

    # Creates NonlinearProblem.
    u0 = rnd_u0(steady_state_network_2, rng; min = 1.0)
    ps = rnd_ps(steady_state_network_2, rng)
    nlprob = NonlinearProblem(steady_state_network_2, u0, ps)
    
    # Solves it using standard algorithm and simulation based algorithm.
    sol1 = solve(nlprob; abstol=1e-12, reltol=1e-12)
    sol2 = solve(nlprob, DynamicSS(Rosenbrock23()); abstol=1e-12, reltol=1e-12)
    
    # Computes NonlinearFunction (manually and automatically).
    nlprob_eval_1 = remake(nlprob; u0 = [:X => sol1[1]])
    nlprob_eval_2 = remake(nlprob; u0 = [:X => sol2[1]])
    function nf_manual(u,p)
        X = u[:X]
        v, K, n, d = p
        return v/10 + v * X^n/(X^n + K^n) - d*X
    end

    # Tests solutions are correct.
    @test nlprob_eval_1.f(nlprob_eval_1.u0, nlprob_eval_1.p)[1] ≈ 0.0 atol=1e-10
    @test nlprob_eval_2.f(nlprob_eval_2.u0, nlprob_eval_2.p)[1] ≈ 0.0 atol=1e-10
    @test nf_manual(sol1, last.(ps)) ≈ 0.0 atol=1e-10
    @test nf_manual(sol2, last.(ps)) ≈ 0.0 atol=1e-10
end

# Checks for system with conservation laws.
# Checks using interfacing with output solution.
let 
    # Conservation laws currently broken (you get stuck in an infinite loop in MTK or something).
    return (@test_broken false)
    # Creates steady state network, unpack the parameter values.
    steady_state_network_3 = @reaction_network begin
        (p,d), 0 <--> X
        (k1, k2), 2Y <--> Y2
        (k3, k4), X + Y2 <--> XY2
    end
    @unpack X, Y, Y2, XY2 = steady_state_network_3

    # Creates NonlinearProblem.
    u0 = [steady_state_network_3.X => rand(), steady_state_network_3.Y => rand() + 1.0, steady_state_network_3.Y2 => rand() + 3.0, steady_state_network_3.XY2 => 0.0]
    p = [:p => rand()+1.0, :d => 0.5, :k1 => 1.0, :k2 => 2.0, :k3 => 3.0, :k4 => 4.0]
    nl_prob_1 = NonlinearProblem(steady_state_network_3, u0, p; remove_conserved = true)
    nl_prob_2 = NonlinearProblem(steady_state_network_3, u0, p)

    # Solves it using standard algorithm and simulation based algorithm.
    sol1 = solve(nl_prob_1; abstol=1e-12, reltol=1e-12)
    sol2 = solve(nl_prob_2, DynamicSS(Rosenbrock23()); abstol=1e-12, reltol=1e-12)

    # Checks output using the ODE's drift function
    nfunc = NonlinearFunction(convert(NonlinearSystem, steady_state_network_3))   
    @test f_eval([:X => sol1[X], :Y => sol1[Y], :Y2 => sol1[Y2], :XY2 => sol1[XY2]], p, 0.0) ≈ [0.0, 0.0, 0.0, 0.0] atol=1e-10
    @test f_eval([:X => sol2[X], :Y => sol2[Y], :Y2 => sol2[Y2], :XY2 => sol2[XY2]], p, 0.0) ≈ [0.0, 0.0, 0.0, 0.0] atol=1e-10
end