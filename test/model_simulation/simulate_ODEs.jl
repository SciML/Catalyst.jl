### Fetch Packages and Reaction Networks ###

# Fetch packages.
using Catalyst, OrdinaryDiffEq, Random, Test
using ModelingToolkit: get_states, get_ps

# Sets rnd number.
using StableRNGs
rng = StableRNG(12345)

# Fetch test networks.
include("../test_networks.jl")

### Compares to Known Solution ###

# Exponential decay, should be identical to the (known) analytical solution.
let
    exponential_decay = @reaction_network begin d, X → ∅ end

    for factor in [1e-2, 1e-1, 1e0, 1e1, 1e2]
        u0 = factor * rand(rng, length(get_states(exponential_decay)))
        p = factor * rand(rng, length(get_ps(exponential_decay)))
        prob = ODEProblem(exponential_decay, u0, (0.0, 100 / factor), p)
        sol = solve(prob, Rosenbrock23(), saveat = range(0.0, 100 / factor, length = 101))
        analytic_sol = map(t -> u0[1] * exp(-p[1] * t),
                           range(0.0, 100 / factor, length = 101))
        @test all(abs.(first.(sol.u) .- analytic_sol) .< 0.1)
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

    for factor in [1e-1, 1e0, 1e1, 1e2, 1e3]
        u0 = factor * rand(rng, length(get_states(known_equilibrium)))
        p = 0.01 .+ factor * rand(rng, length(get_ps(known_equilibrium)))
        prob = ODEProblem(known_equilibrium, u0, (0.0, 1000000.0), p)
        sol = solve(prob, Rosenbrock23())
        @test abs.(sol.u[end][1] / sol.u[end][2] - p[2] / p[1]) < 10000 * eps()
        @test abs.(sol.u[end][3] * sol.u[end][4] / sol.u[end][5] - p[4] / p[3]) <
              10000 * eps()
        @test abs.((sol.u[end][6]^2 / factorial(2)) / (sol.u[end][7]^3 / factorial(3)) -
                   p[6] / p[5]) < 1e-8
        @test abs.(sol.u[end][8] - p[7] / p[8]) < 10000 * eps()
    end
end

### Compares to Known ODE Function ###

identical_networks_1 = Vector{Pair}()

let
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
            u0 = factor * rand(rng, length(get_states(networks[1])))
            p = factor * rand(rng, length(get_ps(networks[1])))
            (i == 3) && (p = min.(round.(p) .+ 1, 10))                      #If parameter in exponent, want to avoid possibility of (-small u)^(decimal). Also avoid large exponents.
            prob1 = ODEProblem(networks[1], u0, (0.0, 10000.0), p)
            sol1 = solve(prob1, Rosenbrock23(), saveat = 1.0)
            prob2 = ODEProblem(networks[2], u0, (0.0, 10000.0), p)
            sol2 = solve(prob2, Rosenbrock23(), saveat = 1.0)
            @test all(abs.(hcat((sol1.u .- sol2.u)...)) .< 1e-7)
        end
    end
end

### Checks Simulations Don't Error ###

let
    for (i, network) in enumerate(reaction_networks_all)
        (i % 5 == 0) &&
            println("Iteration " * string(i) * " at line 104 in file solve_ODEs.jl")
        for factor in [1e-1, 1e0, 1e1]
            u0 = factor * rand(rng, length(get_states(network)))
            p = factor * rand(rng, length(get_ps(network)))
            in(i, [[11:20...]..., 34, 37, 42]) && (p = min.(round.(p) .+ 1, 10))  #If parameter in exponent, want to avoid possibility of (-small u)^(decimal). Also avoid large exponents.
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
        u0 = factor * rand(rng, length(get_states(no_param_network)))
        prob = ODEProblem(no_param_network, u0, (0.0, 1000.0))
        sol = solve(prob, Rosenbrock23())
        @test abs.(sol.u[end][1] - 1.5 / 2) < 1e-8
    end
end

# Test solving with floating point stoichiometry.
let
    function oderhs(du, u, p, t)
        du[1] = -2.5 * p[1] * u[1]^2.5
        du[2] = 3 * p[1] * u[1]^2.5
        nothing
    end
    rn = @reaction_network begin k, 2.5 * A --> 3 * B end
    u0 = [:A => 1.0, :B => 0.0]
    tspan = (0.0, 1.0)
    p = [:k => 1.0]
    oprob = ODEProblem(rn, u0, tspan, p; combinatoric_ratelaws = false)
    du1 = du2 = zeros(2)
    u = rand(2)
    oprob.f(du1, u, [1.0], 0.0)
    oderhs(du2, u, [1.0], 0.0)
    @test isapprox(du1, du2, rtol = 1e3 * eps())
end
