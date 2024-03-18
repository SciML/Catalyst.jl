
### Fetch Packages and Set Global Variables ###
using DiffEqBase, Catalyst, Random, Symbolics, Test
using ModelingToolkit: get_unknowns, get_ps

using StableRNGs
rng = StableRNG(12345)

### Tests Custom Functions ###
let
    new_hill(x, v, k, n) = v * x^n / (k^n + x^n)
    new_poly(x, p1, p2) = 0.5 * p1 * x^2
    new_exp(x, p) = exp(-p * x)

    custom_function_network_1 = @reaction_network begin
        hill(X1, v1, K1, 2), X1 + Y1 --> Z1
        mm(X2, v2, K2), X2 + Y2 --> Z2
        p1 * X3^2 + p2, X3 + Y3 --> Z3
        exp(-p3 * Y4), X4 + Y4 --> Z4
        hillr(X5, v3, K3, 2), X5 + Y5 --> Z5
        mmr(X6, v4, K4), X6 + Y6 --> Z6
        hillar(X7, Y7, v5, K5, 2), X7 + Y7 --> Z7
    end

    custom_function_network_2 = @reaction_network begin
        new_hill(X1, v1, K1, 2), X1 + Y1 --> Z1
        v2 * X2 / (X2 + K2), X2 + Y2 --> Z2
        2 * new_poly(X3, p1, p2) + p2, X3 + Y3 --> Z3
        new_exp(Y4, p3), X4 + Y4 --> Z4
        v3 * (K3^2) / (K3^2 + X5^2), X5 + Y5 --> Z5
        v4 * K4 / (X6 + K4), X6 + Y6 --> Z6
        v5 * (X7^2) / (K5^2 + X7^2 + Y7^2), X7 + Y7 --> Z7
    end

    function permute_ps(pvals, rn1, rn2)
        ps1 = parameters(rn1)
        ps2 = parameters(rn2)
        pvals2 = similar(pvals)
        for (i, p) in enumerate(ps2)
            pidx = findfirst(isequal(p), ps1)
            pvals2[i] = pvals[pidx]
        end
        pvals2
    end

    f1 = ODEFunction(convert(ODESystem, custom_function_network_1), jac = true)
    f2 = ODEFunction(convert(ODESystem, custom_function_network_2), jac = true)
    g1 = SDEFunction(convert(SDESystem, custom_function_network_1))
    g2 = SDEFunction(convert(SDESystem, custom_function_network_2))
    for factor in [1e-2, 1e-1, 1e0, 1e1, 1e2]
        u0 = factor * rand(rng, length(get_unknowns(custom_function_network_1)))
        p = factor * rand(rng, length(get_ps(custom_function_network_2)))

        # needed as this code assumes an ordering of the parameters and species...
        p2 = permute_ps(p, custom_function_network_1, custom_function_network_2)

        t = rand(rng)
        @test all(abs.(f1(u0, p, t) .- f2(u0, p2, t)) .< 10e-10)
        @test all(abs.(f1.jac(u0, p, t) .- f2.jac(u0, p2, t)) .< 10e-10)
        @test all(abs.(g1(u0, p, t) .- g2(u0, p2, t)) .< 10e-10)
    end
end

### Tests that the various notations gives identical results ###

# Michaelis-Menten function.
let
    mm_network = @reaction_network begin
        (1.0, 1.0), 0 ↔ X
        mm(X, v, K), 0 --> X1
        mm(X, v, K), 0 --> X2
        mm(X, v, K), 0 --> X3
    end
    f_mm = ODEFunction(convert(ODESystem, mm_network), jac = true)

    u0 = 10 * rand(rng, length(get_unknowns(mm_network)))
    p = 10 * rand(rng, length(get_ps(mm_network)))
    t = 10 * rand(rng)

    f_mm_output = f_mm(u0, p, t)[2:end]
    f_mm_jac_output = f_mm.jac(u0, p, t)[2:end, 1]
    @test (maximum(f_mm_output) - minimum(f_mm_output)) .< 100 * eps()
    @test (maximum(f_mm_jac_output) - minimum(f_mm_jac_output)) .< 100 * eps()
end

# Repressing Michaelis-Menten function.
let
    mmr_network = @reaction_network begin
        (1.0, 1.0), 0 ↔ X
        mmr(X, v, K), 0 --> X1
        mmr(X, v, K), 0 --> X2
        mmr(X, v, K), 0 --> X3
    end
    f_mmr = ODEFunction(convert(ODESystem, mmr_network), jac = true)

    u0 = 10 * rand(rng, length(get_unknowns(mmr_network)))
    p = 10 * rand(rng, length(get_ps(mmr_network)))
    t = 10 * rand(rng)

    f_mmr_output = f_mmr(u0, p, t)[2:end]
    f_mmr_jac_output = f_mmr.jac(u0, p, t)[2:end, 1]
    @test (maximum(f_mmr_output) - minimum(f_mmr_output)) .< 100 * eps()
    @test (maximum(f_mmr_jac_output) - minimum(f_mmr_jac_output)) .< 100 * eps()
end

# Hill function.
let
    hill_network = @reaction_network begin
        (1.0, 1.0), 0 ↔ X
        hill(X, v, K, 2), 0 --> X1
        hill(X, v, K, 2), 0 --> X2
    end
    f_hill = ODEFunction(convert(ODESystem, hill_network), jac = true)

    u0 = 10 * rand(rng, length(get_unknowns(hill_network)))
    p = 10 * rand(rng, length(get_ps(hill_network)))
    t = 10 * rand(rng)

    f_hill_output = f_hill(u0, p, t)[2:end]
    f_hill_jac_output = f_hill.jac(u0, p, t)[2:end, 1]
    @test (maximum(f_hill_output) - minimum(f_hill_output)) .< 100 * eps()
    @test (maximum(f_hill_jac_output) - minimum(f_hill_jac_output)) .< 100 * eps()
end

# Repressing Hill function.
let
    hillr_network = @reaction_network begin
        (1.0, 1.0), 0 ↔ X
        hillr(X, v, K, 2), 0 --> X1
        hillr(X, v, K, 2), 0 --> X2
    end
    f_hillr = ODEFunction(convert(ODESystem, hillr_network), jac = true)

    u0 = 10 * rand(rng, length(get_unknowns(hillr_network)))
    p = 10 * rand(rng, length(get_ps(hillr_network)))
    t = 10 * rand(rng)

    f_hillr_output = f_hillr(u0, p, t)[2:end]
    f_hillr_jac_output = f_hillr.jac(u0, p, t)[2:end, 1]
    @test (maximum(f_hillr_output) - minimum(f_hillr_output)) .< 100 * eps()
    @test (maximum(f_hillr_jac_output) - minimum(f_hillr_jac_output)) .< 100 * eps()
end

# Activation/repressing Hill function.
let
    hillar_network = @reaction_network begin
        (1.0, 1.0), 0 ↔ (X, Y)
        hillar(X, Y, v, K, 2), 0 --> X1
        hillar(X, Y, v, K, 2), 0 --> X2
        hillar(X, Y, v, K, 2), 0 --> X3
        hillar(X, Y, v, K, 2), 0 --> X4
    end
    f_hillar = ODEFunction(convert(ODESystem, hillar_network), jac = true)

    u0 = 10 * rand(rng, length(get_unknowns(hillar_network)))
    p = 10 * rand(rng, length(get_ps(hillar_network)))
    t = 10 * rand(rng)

    f_hillar_output = f_hillar(u0, p, t)[3:end]
    f_hillar_jac_output = f_hillar.jac(u0, p, t)[3:end, 1]
    @test (maximum(f_hillar_output) - minimum(f_hillar_output)) .< 100 * eps()
    @test (maximum(f_hillar_jac_output) - minimum(f_hillar_jac_output)) .< 100 * eps()
end

### Test Symbolic Derivatives ###

let
    @variables X Y
    @parameters v K n

    @test isequal(Symbolics.derivative(Catalyst.mm(X, v, K), X), v * K / (K + X)^2)
    @test isequal(Symbolics.derivative(Catalyst.mm(X, v, K), v), X / (K + X))
    @test isequal(Symbolics.derivative(Catalyst.mm(X, v, K), K), -v * X / (K + X)^2)

    @test isequal(Symbolics.derivative(Catalyst.mmr(X, v, K), X), -v * K / (K + X)^2)
    @test isequal(Symbolics.derivative(Catalyst.mmr(X, v, K), v), K / (K + X))
    @test isequal(Symbolics.derivative(Catalyst.mmr(X, v, K), K), v * X / (K + X)^2)

    @test isequal(Symbolics.derivative(Catalyst.hill(X, v, K, n), X),
                  n * v * (K^n) * (X^(n - 1)) / (K^n + X^n)^2)
    @test isequal(Symbolics.derivative(Catalyst.hill(X, v, K, n), v), X^n / (K^n + X^n))
    @test isequal(Symbolics.derivative(Catalyst.hill(X, v, K, n), K),
                  -n * v * (K^(n - 1)) * (X^n) / (K^n + X^n)^2)
    @test isequal(Symbolics.derivative(Catalyst.hill(X, v, K, n), n),
                  v * (X^n) * (K^n) * (log(X) - log(K)) / (K^n + X^n)^2)

    @test isequal(Symbolics.derivative(Catalyst.hillr(X, v, K, n), X),
                  -n * v * (K^n) * (X^(n - 1)) / (K^n + X^n)^2)
    @test isequal(Symbolics.derivative(Catalyst.hillr(X, v, K, n), v), K^n / (K^n + X^n))
    @test isequal(Symbolics.derivative(Catalyst.hillr(X, v, K, n), K),
                  n * v * (K^(n - 1)) * (X^n) / (K^n + X^n)^2)
    @test isequal(Symbolics.derivative(Catalyst.hillr(X, v, K, n), n),
                  v * (X^n) * (K^n) * (log(K) - log(X)) / (K^n + X^n)^2)

    @test isequal(Symbolics.derivative(Catalyst.hillar(X, Y, v, K, n), X),
                  n * v * (K^n + Y^n) * (X^(n - 1)) / (K^n + X^n + Y^n)^2)
    @test isequal(Symbolics.derivative(Catalyst.hillar(X, Y, v, K, n), Y),
                  -n * v * (Y^(n - 1)) * (X^n) / (K^n + X^n + Y^n)^2)
    @test isequal(Symbolics.derivative(Catalyst.hillar(X, Y, v, K, n), v),
                  X^n / (K^n + X^n + Y^n))
    @test isequal(Symbolics.derivative(Catalyst.hillar(X, Y, v, K, n), K),
                  -n * v * (v^(n - 1)) * (X^n) / (K^n + X^n + Y^n)^2)
    @test isequal(Symbolics.derivative(Catalyst.hillar(X, Y, v, K, n), n),
                  v * (X^n) * ((K^n + Y^n) * log(X) - (K^n) * log(K) - (Y^n) * log(Y)) /
                  (K^n + X^n + Y^n)^2)
end

### Tests Current Function Expansion ###

# Tests `ReactionSystem`s.
let
    @variables t 
    @species x(t) y(t)
    @parameters k v n 
    rs1 = @reaction_network rs begin
        mm(x, v, k), 0 --> x
        mmr(x, v, k), 0 --> x
        hill(x, v, k, n), 0 --> x
        hillr(x, v, k, n), 0 --> x
        hillar(x, y, v, k, n), 0 --> x    
    end
    rs2 = @reaction_network rs begin
        v * x / (x + k), 0 --> x
        v * k / (x + k), 0 --> x
        v * (x^n) / (x^n + k^n), 0 --> x
        v * (k^n) / (x^n + k^n), 0 --> x
        v * (x^n) / (x^n + y^n + k^n), 0 --> x    
    end

    @test Catalyst.expand_registered_functions(rs1) == rs2
end

# Tests `Reaction`s.
let
    @variables t 
    @species x(t) y(t)
    @parameters k v n 
    
    r1 = @reaction mm(x, v, k), 0 --> x
    r2 = @reaction mmr(x, v, k), 0 --> x
    r3 = @reaction hill(x, v, k, n), 0 --> x
    r4 = @reaction hillr(x, v, k, n), 0 --> x
    r5 = @reaction hillar(x, y, v, k, n), 0 --> x + y
    
    @test isequal(Catalyst.expand_registered_functions(r1).rate, v * x / (x + k))
    @test isequal(Catalyst.expand_registered_functions(r2).rate, v * k / (x + k))
    @test isequal(Catalyst.expand_registered_functions(r3).rate, v * (x^n) / (x^n + k^n))
    @test isequal(Catalyst.expand_registered_functions(r4).rate, v * (k^n) / (x^n + k^n))
    @test isequal(Catalyst.expand_registered_functions(r5).rate, v * (x^n) / (x^n + y^n + k^n))
end

# Tests `Equation`s.
let
    @variables T X(T) Y(T)
    @parameters K V N 
    
    eq1 = 0 ~ mm(X, V, K)
    eq2 = 0 ~ mmr(X, V, K)
    eq3 = 0 ~ hill(X, V, K, N)
    eq4 = 0 ~ hillr(X, V, K, N)
    eq5 = 0 ~ hillar(X, Y, V, K, N)
    
    @test isequal(Catalyst.expand_registered_functions(eq1), 0 ~ V * X / (X + K))
    @test isequal(Catalyst.expand_registered_functions(eq2), 0 ~ V * K / (X + K))
    @test isequal(Catalyst.expand_registered_functions(eq3), 0 ~ V * (X^N) / (X^N + K^N))
    @test isequal(Catalyst.expand_registered_functions(eq4), 0 ~ V * (K^N) / (X^N + K^N))
    @test isequal(Catalyst.expand_registered_functions(eq5), 0 ~ V * (X^N) / (X^N + Y^N + K^N))
end