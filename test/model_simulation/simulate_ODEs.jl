### Fetch Packages and Reaction Networks ###

# Fetch packages.
using Catalyst, OrdinaryDiffEq, Random, Test
using ModelingToolkit: get_unknowns, get_ps

# Sets rnd number.
using StableRNGs
rng = StableRNG(12345)

# Fetch test networks and functions.
include("../test_networks.jl")
include("../test_functions.jl")

### Compares to Known Solution ###

# Exponential decay, should be identical to the (known) analytical solution.
let
    exponential_decay = @reaction_network begin d, X → ∅ end

    for factor in [1e-2, 1e-1, 1e0, 1e1, 1e2]
        u0 = rnd_u0(exponential_decay, rng; factor)    
        t_stops = range(0.0, 100 / factor, length = 101)
        p = rnd_ps(exponential_decay, rng; factor)
        prob = ODEProblem(exponential_decay, u0, (0.0, t_stops[end]), p)

        sol = solve(prob, Rosenbrock23(), saveat = t_stops, abstol = 1e-14, reltol = 1e-14)
        analytic_sol = [u0[1][2] * exp(-p[1][2] * t) for t in t_stops]
        @test sol[:X] ≈ analytic_sol
    end
end

# Networks with know equilibrium.
let
    known_equilibrium = @reaction_network begin
        (k1, k2), X1 ↔ X2
        (k3, k4), X3 + X4 ↔ X5
        (k5, k6), 2X6 ↔ 3X7
        (k7, k8), ∅ ↔ X8
    end

    for factor in [1e-1, 1e0, 1e1]
        u0 = rnd_u0(known_equilibrium, rng; factor)    
        p = rnd_ps(known_equilibrium, rng; factor, min = 0.1)
        prob = ODEProblem(known_equilibrium, u0, (0.0, 100000.0), p)
        sol = solve(prob, Vern7(); abstol = 1e-12, reltol = 1e-12)

        @test sol[:X1][end] / sol[:X2][end] ≈ prob.ps[:k2] / prob.ps[:k1] atol=1e-8
        @test sol[:X3][end] * sol[:X4][end] / sol[:X5][end] ≈ prob.ps[:k4] / prob.ps[:k3] atol=1e-8
        @test (sol[:X6][end]^2 / factorial(2)) / (sol[:X7][end]^3 / factorial(3)) ≈ prob.ps[:k6] / prob.ps[:k5] atol=1e-8
        @test sol[:X8][end] ≈ prob.ps[:k7] / prob.ps[:k8] atol=1e-8
    end
end

### Compares to Known ODE Function ###

let
    identical_networks_1 = Vector{Pair}()
    
    function real_functions_1(du, u, p, t)
        X1, X2, X3 = u
        p1, p2, p3, k1, k2, k3, k4, d1, d2, d3 = p
        du[1] = p1 + k1 * X2 - k2 * X1 * X3^2 / factorial(2) - k3 * X1 + k4 * X3 - d1 * X1
        du[2] = p2 - k1 * X2 + k2 * X1 * X3^2 / factorial(2) - d2 * X2
        du[3] = p3 + 2 * k1 * X2 - 2 * k2 * X1 * X3^2 / factorial(2) + k3 * X1 - k4 * X3 -
                d3 * X3
    end
    push!(identical_networks_1, reaction_networks_standard[1] => real_functions_1)

    function real_functions_2(du, u, p, t)
        X1, X2 = u
        v1, K1, v2, K2, d = p
        du[1] = v1 * K1 / (K1 + X2) - d * X1 * X2
        du[2] = v2 * X1 / (K2 + X1) - d * X1 * X2
    end
    push!(identical_networks_1, reaction_networks_standard[2] => real_functions_2)

    function real_functions_3(du, u, p, t)
        X1, X2, X3 = u
        v1, v2, v3, K1, K2, K3, n1, n2, n3, d1, d2, d3 = p
        du[1] = v1 * K1^n1 / (K1^n1 + X3^n1) - d1 * X1
        du[2] = v2 * K2^n2 / (K2^n2 + X1^n2) - d2 * X2
        du[3] = v3 * K3^n3 / (K3^n3 + X2^n3) - d3 * X3
    end
    push!(identical_networks_1, reaction_networks_hill[2] => real_functions_3)

    function real_functions_4(du, u, p, t)
        X1, X2, X3 = u
        k1, k2, k3, k4, k5, k6 = p
        du[1] = -k1 * X1 + k2 * X2 + k5 * X3 - k6 * X1
        du[2] = -k3 * X2 + k4 * X3 + k1 * X1 - k2 * X2
        du[3] = -k5 * X3 + k6 * X1 + k3 * X2 - k4 * X3
    end
    push!(identical_networks_1, reaction_networks_constraint[1] => real_functions_4)

    function real_functions_5(du, u, p, t)
        X, Y, Z = u
        k1, k2, k3, k4 = p
        du[1] = k1 - k2 * log(12 + X) * X
        du[2] = k2 * log(12 + X) * X - k3 * log(3 + Y) * Y
        du[3] = k3 * log(3 + Y) * Y - log(5, 6 + k4) * Z
    end
    push!(identical_networks_1, reaction_networks_weird[2] => real_functions_5)

    for (i, networks) in enumerate(identical_networks_1)
        for factor in [1e-2, 1e-1, 1e0, 1e1]
            u0_1 = rnd_u0(networks[1], rng; factor)
            p_1 = rnd_ps(networks[1], rng; factor)
            (i == 3) && (p_1 =rnd_ps_Int64(networks[1], rng))
            u0_2 = last.(u0_1)
            p_2 = last.(p_1)
                                
            prob1 = ODEProblem(networks[1], u0_1, (0.0, 10000.0), p_1)
            prob2 = ODEProblem(networks[2], u0_2, (0.0, 10000.0), p_2)
            sol1 = solve(prob1, Rosenbrock23(); abstol = 1e-10, reltol = 1e-10, saveat = 1.0)
            sol2 = solve(prob2, Rosenbrock23(); abstol = 1e-10, reltol = 1e-10, saveat = 1.0)
            @test sol1.u ≈ sol2.u
        end
    end
end

### Checks Simulations Don't Error ###
let
    for (i, network) in enumerate(reaction_networks_all)
        for factor in [1e-1, 1e0, 1e1]
            u0 = rnd_u0(network, rng; factor)
            #If parameter in exponent, want to avoid possibility of (-small u)^(decimal). Also avoid large exponents.
            if in(i, [[11:20...]..., 34, 37, 42])
                p = rnd_ps(network, rng)
            else
                p = rnd_ps(network, rng; factor)
            end
            prob = ODEProblem(network, u0, (0.0, 1.0), p)
            @test SciMLBase.successful_retcode(solve(prob, Rosenbrock23()))
        end
    end
end

### Other Tests ###

# No parameter test.
let
    no_param_network = @reaction_network begin (1.5, 2), ∅ ↔ X end
    for factor in [1e0, 1e1, 1e2]
        u0 = rnd_u0(no_param_network, rng; factor)
        prob = ODEProblem(no_param_network, u0, (0.0, 1000.0))
        sol = solve(prob, Rosenbrock23())
        @test sol[:X][end] ≈ 1.5 / 2.0
    end
end

# Test solving with floating point stoichiometry.
let
    # Prepare model.
    function oderhs(du, u, p, t)
        du[1] = -2.5 * p[1] * u[1]^2.5
        du[2] = 3 * p[1] * u[1]^2.5
        nothing
    end
    rn = @reaction_network begin 
        k, 2.5 * A --> 3 * B 
    end
    u = rnd_u0(rn, rng)
    tspan = (0.0, 1.0)
    p = [:k => 1.0]
    
    # Check equivalence.
    du1 = du2 = zeros(2) 
    oprob = ODEProblem(rn, u, tspan, p; combinatoric_ratelaws = false)
    oprob.f(du1, oprob.u0, oprob.p, 90.0)
    oderhs(du2, last.(u), last.(p), 0.0)
    @test du1 ≈ du2
end