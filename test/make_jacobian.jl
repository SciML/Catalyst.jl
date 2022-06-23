### Fetch required packages and reaction networks ###
using Catalyst, DiffEqBase, Random, Test

using StableRNGs
rng = StableRNG(12345)

### Checks that the jacobian is correct for networks without parameters ###
jacobian_network_1 = @reaction_network begin
    (2.0, 1.0), ∅ ↔ X
    (3.0, 1.0), ∅ ↔ Y
    (5.0, 2.0), X + Y ↔ XY
end
test_jac = ODEFunction(convert(ODESystem, jacobian_network_1); jac = true).jac([1.0, 1.0, 1.0], [],
                                                                               0.0)
real_jac = [-6.0 -5.0 2.0; -5.0 -6.0 2.0; 5.0 5.0 -2.0]
@test test_jac == real_jac

### Checks that the jacobian is correct for networks with parameters ###
jacobian_network_2 = @reaction_network begin
    (p1, 1.0), ∅ ↔ X
    (p2, 1.0), ∅ ↔ Y
    (p3 * X, 1.0), X + Y ↔ XY
end p1 p2 p3
for factor in [1e-2, 1e-1, 1e0, 1e1, 1e2], repeat in 1:10
    u = factor * rand(rng, 3)
    p = factor * rand(rng, 3)
    local test_jac = ODEFunction(convert(ODESystem, jacobian_network_2); jac = true).jac(u, p, 0.0)
    local real_jac = [-1-2 * p[3] * u[2] * u[1] -p[3]*u[1]*u[1] 1.0;
                      -2*p[3]*u[2]*u[1] -1-p[3] * u[1] * u[1] 1;
                      2*p[3]*u[2]*u[1] p[3]*u[1]*u[1] -1.0]
    @test all(abs.(test_jac .- real_jac) .< 1e-9)
end

jacobian_network_3 = @reaction_network begin
    k1, 2A → B
    k2, B → 2A
    k3, A + B → C
    k4, C → A + B
    k5, 3C → 3A
    k6, 0 → 2B
    hill(A, k7, k8, 2), ∅ → B
end k1 k2 k3 k4 k5 k6 k7 k8
for factor in [1e-2, 1e-1, 1e0, 1e1, 1e2], repeat in 1:10
    u = factor * rand(rng, 3)
    p = factor * rand(rng, 8)
    A, B, C = u
    k1, k2, k3, k4, k5, k6, k7, k8 = p
    local test_jac = ODEFunction(convert(ODESystem, jacobian_network_3); jac = true).jac(u, p, 0.0)
    local real_jac = [-2 * k1 * A-k3 * B 2 * k2-k3 * A k4+3 * k5 * C^2 / 2;
                      k1 * A - k3 * B+2 * k7 * k8^2 * A / (k8^2 + A^2)^2 -k2-k3 * A k4;
                      k3*B k3*A -k4-3 * k5 * C^2 / 2]
    @test all(abs.(test_jac .- real_jac) .< 10e-9)
end
