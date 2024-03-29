### Fetch Packages and Set Global Variables ###

### Prepares Tests ###

# Fetch packages.
using DiffEqBase, Catalyst, JumpProcesses, Random, Statistics, Test
using ModelingToolkit: get_unknowns, get_ps

# Sets stable rng number.
using StableRNGs
rng = StableRNG(12345)

# Fetch test functions.
include("../test_functions.jl")

# Declares a network used throughout all tests.
higher_order_network_1 = @reaction_network begin
    p, ∅ ⟼ X1
    r1, 2X1 ⟼ 3X2
    mm(X1, r2, K), 3X2 ⟼ X3 + 2X4
    r3, X3 + 2X4 ⟼ 3X5 + 3X6
    r4 * X2, 3X5 + 3X6 ⟼ 3X5 + 2X7 + 4X8
    r5, 3X5 + 2X7 + 4X8 ⟼ 10X9
    r6, 10X9 ⟼ X10
    d, 2X10 ⟼ ∅
end

### Basic Tests ###

# Tests that ODE and SDE functions are correct (by comparing to network with manually written higher order rates). 
let
    higher_order_network_2 = @reaction_network begin
        p, ∅ ⟾ X1
        r1 * X1^2 / factorial(2), 2X1 ⟾ 3X2
        mm(X1, r2, K) * X2^3 / factorial(3), 3X2 ⟾ X3 + 2X4
        r3 * X3 * X4^2 / factorial(2), X3 + 2X4 ⟾ 3X5 + 3X6
        r4 * X2 * X5^3 * X6^3 / (factorial(3) * factorial(3)), 3X5 + 3X6 ⟾ 3X5 + 2X7 + 4X8
        r5 * X5^3 * X7^2 * X8^4 / (factorial(3) * factorial(2) * factorial(4)),
        3X5 + 2X7 + 4X8 ⟾ 10X9
        r6 * X9^10 / factorial(10), 10X9 ⟾ X10
        d * X10^2 / factorial(2), 2X10 ⟾ ∅
    end

    for factor in [1e-1, 1e0, 1e1, 1e2]
        u0 = rnd_u0(higher_order_network_1, rng; factor)
        ps = rnd_ps(higher_order_network_1, rng; factor)
        t = rand(rng)

        @test f_eval(higher_order_network_1, u0, ps, t) == f_eval(higher_order_network_2, u0, ps, t)
        @test jac_eval(higher_order_network_1, u0, ps, t) == jac_eval(higher_order_network_2, u0, ps, t)
        @test g_eval(higher_order_network_1, u0, ps, t) == g_eval(higher_order_network_2, u0, ps, t)
    end
end

# Tests that Jump Systems are correct (by comparing to network with manually written higher order rates). 
# Currently fails because binomial only takes Int input (and X is Float64).
# I see several solutions, but depends on whether we allow species to be specified as Int64.
# I have marked this one as broken for now.
@test_broken let 
    higher_order_network_3 = @reaction_network begin
        p, ∅ ⟼ X1
        r1 * binomial(X1, 2), 2X1 ⟾ 3X2
        mm(X1, r2, K) * binomial(X2, 3), 3X2 ⟾ X3 + 2X4
        r3 * binomial(X3, 1) * binomial(X4, 2), X3 + 2X4 ⟾ 3X5 + 3X6
        r4 * X2 * binomial(X5, 3) * binomial(X6, 3), 3X5 + 3X6 ⟾ 3X5 + 2X7 + 4X8
        r5 * binomial(X5, 3) * binomial(X7, 2) * binomial(X8, 4), 3X5 + 2X7 + 4X8 ⟾ 10X9
        r6 * binomial(X9, 10), 10X9 ⟾ X10
        d * binomial(X10, 2), 2X10 ⟾ ∅
    end

    for factor in [1e-1, 1e0]
        u0 = rand(rng, 1:Int64(factor * 100), length(get_unknowns(higher_order_network_1)))
        p = factor * rand(rng, length(get_ps(higher_order_network_3)))
        prob1 = JumpProblem(higher_order_network_1,
                            DiscreteProblem(higher_order_network_1, u0, (0.0, 1000.0), p),
                            Direct(); rng)
        sol1 = solve(prob1, SSAStepper())
        prob2 = JumpProblem(higher_order_network_3,
                            DiscreteProblem(higher_order_network_3, u0, (0.0, 1000.0), p),
                            Direct(); rng)
        sol2 = solve(prob2, SSAStepper())
        for i in 1:length(u0)
            vals1 = getindex.(sol1.u, i)
            vals2 = getindex.(sol1.u, i)
            (mean(vals2) > 0.001) && @test 0.8 < mean(vals1) / mean(vals2) < 1.25
            (std(vals2) > 0.001) && @test 0.8 < std(vals1) / std(vals2) < 1.25
        end
    end
end
