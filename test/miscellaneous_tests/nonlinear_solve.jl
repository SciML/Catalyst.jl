### Fetch Packages ###

# Fetch packages.
using Catalyst, NonlinearSolve, OrdinaryDiffEq, SteadyStateDiffEq
using Random, Test

# Sets rnd number.
using StableRNGs
rng = StableRNG(12345)

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
    u0 = rand(rng, length(states(steady_state_network_1)))
    p = rand(rng, length(parameters(steady_state_network_1)))
    nl_prob = NonlinearProblem(steady_state_network_1, u0, p)
    
    # Solves it using standard algorithm and simulation based algorithm.
    sol1 = solve(nl_prob; abstol=1e-12, reltol=1e-12).u
    sol2 = solve(nl_prob, DynamicSS(Rosenbrock23()); abstol=1e-12, reltol=1e-12).u
    
    # Tests solutions are correct.
    @test isapprox(sol1[1], p[1] / p[2]; atol=1e-10)
    @test isapprox(sol1[2]^3 / factorial(3), p[3] / p[4]; atol=1e-10)
    @test isapprox(sol1[3], (p[5] * p[2]) / (p[6] * p[1]); atol=1e-10)
    @test isapprox(sol2[1], p[1] / p[2]; atol=1e-10)
    @test isapprox(sol2[2]^3 / factorial(3), p[3] / p[4]; atol=1e-10)
    @test isapprox(sol2[3], (p[5] * p[2]) / (p[6] * p[1]); atol=1e-10)
end

# Creates a system with multiple steady states.
# Checks that corresponding ODEFunction return 0.0 in both cases. Checks for manually computed function as well.
# Checks with SYmbol `u0` and Symbolic `p`.
let 
    # Creates steady state network, unpack the parameter values.
    steady_state_network_2 = @reaction_network begin
        v / 10 + hill(X, v, K, n), ∅ → X
        d, X → ∅
    end
    @unpack v, K, n, d = steady_state_network_2

    # Creates NonlinearProblem.
    u0 = [:X => 1.0]
    p = [v => 1.0 + rand(rng), K => 0.8, n => 3, d => 1.0]
    nl_prob = NonlinearProblem(steady_state_network_2, u0, p)
    
    # Solves it using standard algorithm and simulation based algorithm.
    sol1 = solve(nl_prob; abstol=1e-12, reltol=1e-12).u
    sol2 = solve(nl_prob, DynamicSS(Rosenbrock23(); abstol=1e-12, reltol=1e-12); abstol=1e-12, reltol=1e-12).u
    
    # Computes NonlinearFunction (manually and automatically).
    nfunc = NonlinearFunction(convert(NonlinearSystem, steady_state_network_2))    
    function nf_manual(u,p)
        X = u[1]
        v, K, n, d = p
        return v/10 + v * X^n/(X^n + K^n) - d*X
    end

    # Tests solutions are correct.
    @test isapprox(nfunc(sol1, last.(p))[1], 0.0; atol=1e-10)
    @test isapprox(nfunc(sol2, last.(p))[1], 0.0; atol=1e-10)
    @test isapprox(nf_manual(sol1, last.(p)), 0.0; atol=1e-10)
    @test isapprox(nf_manual(sol2, last.(p)), 0.0; atol=1e-10)
end

# Checks for system with conservation laws.
# Checks using interfacing with output solution.
let 
    # Creates steady state network, unpack the parameter values.
    steady_state_network_3 = complete(@reaction_network begin
        (p,d), 0 <--> X
        (k1, k2), 2Y <--> Y2
        (k3, k4), X + Y2 <--> XY2
    end)
    @unpack X, Y, Y2, XY2 = steady_state_network_3

    # Creates NonlinearProblem.
    u0 = [steady_state_network_3.X => rand(), steady_state_network_3.Y => rand() + 1.0, steady_state_network_3.Y2 => rand() + 3.0, steady_state_network_3.XY2 => 0.0]
    p = [:p => rand()+1.0, :d => 0.5, :k1 => 1.0, :k2 => 2.0, :k3 => 3.0, :k4 => 4.0]
    nl_prob_1 = NonlinearProblem(steady_state_network_3, u0, p; remove_conserved = true)
    nl_prob_2 = NonlinearProblem(steady_state_network_3, u0, p)

    # Solves it using standard algorithm and simulation based algorithm.
    sol1 = solve(nl_prob_1; abstol=1e-12, reltol=1e-12)
    sol2 = solve(nl_prob_2, DynamicSS(Rosenbrock23(); abstol=1e-12, reltol=1e-12); abstol=1e-12, reltol=1e-12)

    # Checks output using NonlinearFunction.
    nfunc = NonlinearFunction(convert(NonlinearSystem, steady_state_network_3))   
    @test isapprox(nfunc([sol1[X], sol1[Y], sol1[Y2], sol1[XY2]], last.(p)), [0.0, 0.0, 0.0, 0.0]; atol=1e-10)
    @test isapprox(nfunc([sol2[X], sol2[Y], sol2[Y2], sol2[XY2]], last.(p)), [0.0, 0.0, 0.0, 0.0]; atol=1e-10)
end
