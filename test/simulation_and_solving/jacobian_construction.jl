### Prepares Tests ###

# Fetch packages.
using Catalyst, DiffEqBase, Test

# Sets stable rng number.
using StableRNGs
rng = StableRNG(12345)

# Fetch test functions.
include("../test_functions.jl")

### Basic Tests ###

# Checks that the jacobian is correct for networks without parameters.
let
    jacobian_network_1 = @reaction_network begin
        (2.0, 1.0), ∅ ↔ X
        (3.0, 1.0), ∅ ↔ Y
        (5.0, 2.0), X + Y ↔ XY
    end
    
    test_jac = jac_eval(jacobian_network_1, [:X => 1.0, :Y => 1.0, :XY => 1.0], [], 0.0)
    real_jac = [-6.0 -5.0 2.0; -5.0 -6.0 2.0; 5.0 5.0 -2.0]
    @test test_jac == real_jac
end

# Checks that the jacobian is correct for networks with parameters.
let
    jacobian_network_2 = @reaction_network begin
        (p1, 1.0), ∅ ↔ X
        (p2, 1.0), ∅ ↔ Y
        (p3 * X, 1.0), X + Y ↔ XY
    end
    @unpack X, Y, XY, p1, p2, p3 = jacobian_network_2
    
    for factor in [1e-2, 1e-1, 1e0, 1e1, 1e2], repeat in 1:10
        u = Dict(rnd_u0(jacobian_network_2, rng; factor))
        p = Dict(rnd_ps(jacobian_network_2, rng; factor))
        test_jac = jac_eval(jacobian_network_2, u, p, 0.0)
        real_jac = [-1-2 * p[p3] * u[Y] * u[X] -p[p3]*u[X]*u[X] 1.0;
                    -2*p[p3]*u[Y]*u[X] -1-p[p3] * u[X] * u[X] 1;
                    2*p[p3]*u[Y]*u[X] p[p3]*u[X]*u[X] -1.0]
        @test test_jac ≈ real_jac
    end
end

# Checks for a more complicated network, with non-unitary stoichiometries and a hill function.
let
    jacobian_network_3 = @reaction_network begin
        k1, 2A → B
        k2, B → 2A
        k3, A + B → C
        k4, C → A + B
        k5, 3C → 3A
        k6, 0 → 2B
        hill(A, k7, k8, 2), ∅ → B
    end
    @unpack A, B, C, k1, k2, k3, k4, k5, k6, k7, k8 = jacobian_network_3

    for factor in [1e-2, 1e-1, 1e0, 1e1, 1e2], repeat in 1:10
        u = Dict(rnd_u0(jacobian_network_3, rng; factor))
        p = Dict(rnd_ps(jacobian_network_3, rng; factor))
        test_jac = jac_eval(jacobian_network_3, u, p, 0.0)
        real_jac = [-2 * p[k1] * u[A]-p[k3] * u[B] 2 * p[k2]-p[k3] * u[A] p[k4]+3 * p[k5] * u[C]^2 / 2;
                    p[k1] * u[A] - p[k3] * u[B]+2 * p[k7] * p[k8]^2 * u[A] / (p[k8]^2 + u[A]^2)^2 -p[k2]-p[k3] * u[A] p[k4];
                    p[k3]*u[B] p[k3]*u[A] -p[k4]-3 * p[k5] * u[C]^2 / 2]
        @test test_jac ≈ real_jac
    end
end