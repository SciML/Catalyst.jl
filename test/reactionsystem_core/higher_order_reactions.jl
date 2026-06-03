### Prepares Tests ###

# Fetch packages.
using DiffEqBase, Catalyst, JumpProcesses, Statistics, Test

# Sets stable rng number.
using StableRNGs
rng = StableRNG(12345)
seed = rand(rng, 1:100)

# Fetch test functions.
include("../test_functions.jl")

### Basic Tests ###

# Declares a base network to you for comparisons.
base_higher_order_network = @reaction_network begin
    p, ∅ ⟼ X1
    r1, 2X1 ⟼ 3X2
    mm(X1, r2, K), 3X2 ⟼ X3 + 2X4
    r3, X3 + 2X4 ⟼ 3X5 + 3X6
    r4 * X2, 3X5 + 3X6 ⟼ 3X5 + 2X7 + 4X8
    r5, 3X5 + 2X7 + 4X8 ⟼ 10X9
    r6, 10X9 ⟼ X10
    d, 2X10 ⟼ ∅
end

# Tests that ODE and SDE functions are correct (by comparing to network with manually written higher order rates).
let
    higher_order_network_alt1 = @reaction_network begin
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
        u0 = rnd_u0(base_higher_order_network, rng; factor)
        ps = rnd_ps(base_higher_order_network, rng; factor)
        t = rand(rng)

        @test f_eval(base_higher_order_network, u0, ps, t) ≈ f_eval(higher_order_network_alt1, u0, ps, t)
        @test jac_eval(base_higher_order_network, u0, ps, t) ≈ jac_eval(higher_order_network_alt1, u0, ps, t)
        @test g_eval(base_higher_order_network, u0, ps, t) ≈ g_eval(higher_order_network_alt1, u0, ps, t)
    end
end

# Tests that Jump Systems are correct (by comparing to network with manually written higher order rates).
let
    # Declares the reactions using Catalyst, but defines the propensities manually.
    higher_order_network_alt1 = @reaction_network begin
        p, ∅ ⟼ X1
        r1 * binomial(X1, 2), 2X1 ⟾ 3X2
        mm(X1, r2, K) * binomial(X2, 3), 3X2 ⟾ X3 + 2X4
        r3 * binomial(X3, 1) * binomial(X4, 2), X3 + 2X4 ⟾ 3X5 + 3X6
        r4 * X2 * binomial(X5, 3) * binomial(X6, 3), 3X5 + 3X6 ⟾ 3X5 + 2X7 + 4X8
        r5 * binomial(X5, 3) * binomial(X7, 2) * binomial(X8, 4), 3X5 + 2X7 + 4X8 ⟾ 10X9
        r6 * binomial(X9, 10), 10X9 ⟾ X10
        d * binomial(X10, 2), 2X10 ⟾ ∅
    end

    # Declares both the reactions and the propensities manually.
    rate1(u, p, t) = p[1]
    rate2(u, p, t) = p[2] * binomial(u[1], 2)
    rate3(u, p, t) = mm(u[1], p[3], p[4]) * binomial(u[2], 3)
    rate4(u, p, t) = p[5] * binomial(u[3], 1) * binomial(u[4], 2)
    rate5(u, p, t) = p[6] * u[2] * binomial(u[5], 3) * binomial(u[6], 3)
    rate6(u, p, t) = p[7] * binomial(u[5], 3) * binomial(u[7], 2) * binomial(u[8], 4)
    rate7(u, p, t) = p[8] * binomial(u[9], 10)
    rate8(u, p, t) = p[9] * binomial(u[10], 2)

    affect1!(int) = (int.u[1] += 1)
    affect2!(int) = (int.u[1] -= 2; int.u[2] += 3;)
    affect3!(int) = (int.u[2] -= 3; int.u[3] += 1; int.u[4] += 2;)
    affect4!(int) = (int.u[3] -= 1; int.u[4] -= 2; int.u[5] += 3; int.u[6] += 3;)
    affect5!(int) = (int.u[5] -= 3; int.u[6] -= 3; int.u[5] += 3; int.u[7] += 2; int.u[8] += 4;)
    affect6!(int) = (int.u[5] -= 3; int.u[7] -= 2; int.u[8] -= 4; int.u[9] += 10;)
    affect7!(int) = (int.u[9] -= 10; int.u[10] += 1;)
    affect8!(int) = (int.u[10] -= 2;)

    higher_order_network_alt2 = ConstantRateJump.([rate1, rate2, rate3, rate4, rate5, rate6, rate7, rate8],
                                [affect1!, affect2!, affect3!, affect4!, affect5!, affect6!, affect7!, affect8!])

    # Prepares JumpProblem via Catalyst.
    u0_base = rnd_u0_Int64(base_higher_order_network, rng)
    ps_base = rnd_ps(base_higher_order_network, rng)
    jprob_base = JumpProblem(base_higher_order_network, u0_base, (0.0, 1000.0), ps_base; rng = StableRNG(1234))

    # Prepares JumpProblem partially declared manually.
    jprob_alt1 = JumpProblem(higher_order_network_alt1, u0_base, (0.0, 1000.0), ps_base; rng = StableRNG(1234))

    # Prepares JumpProblem via manually declared system.
    u0_alt2 = map_to_vec(u0_base, [:X1, :X2, :X3, :X4, :X5, :X6, :X7, :X8, :X9, :X10])
    ps_alt2 = map_to_vec(ps_base, [:p, :r1, :r2, :K, :r3, :r4, :r5, :r6, :d])
    dprob_alt2 = DiscreteProblem(u0_alt2, (0.0, 1000.0), ps_alt2)
    jprob_alt2 = JumpProblem(dprob_alt2, Direct(), higher_order_network_alt2...; rng = StableRNG(1234))

    # Simulates the models.
    sol_base = solve(jprob_base, SSAStepper(); seed, saveat = 1.0)
    sol_alt1 = solve(jprob_alt1, SSAStepper(); seed, saveat = 1.0)
    sol_alt2 = solve(jprob_alt2, SSAStepper(); seed, saveat = 1.0)

    # Checks that species means in the simulations are similar
    @test mean(sol_base[:X10]) ≈ mean(sol_alt1[:X10]) atol = 1e-1 rtol = 1e-1
    @test mean(sol_alt1[:X10]) ≈ mean(sol_alt2[10,:]) atol = 1e-1 rtol = 1e-1
end
