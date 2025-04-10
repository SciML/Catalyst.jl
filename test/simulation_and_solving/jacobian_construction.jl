### Prepares Tests ###

# Fetch packages.
using Catalyst, DiffEqBase, OrdinaryDiffEqRosenbrock, Test

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

# Checks that the Jacobians (dense and sparse) are identical for different system types.
let
    # Creates model (vaguely messy model without conserved quantities).
    rn = @reaction_network begin
        (p,d), 0 <--> (X,Y)
        (k1,k2), X + Y <--> XY
        (k1,k2), X + 2Y <--> XY2
        (k1,k2), XY + XY2 <--> X2Y3
        d, (XY2,X2Y3) --> 0
        mm(X2Y3,v,K), 0 --> Z
        (k3,k4), 3Z <--> Z3
        1.0, X3 --> 0
    end

    # Performs tests for different randomised values (to be thoroughly sure).
    for factor in [0.1, 1.0, 10.0]
        # Creates randomised species and parameter values. Generates jacobians (dense/sparse).
        u0 = rnd_u0(rn, rng; factor)
        ps = rnd_ps(rn, rng; factor)
        oprob_jac = ODEProblem(rn, u0, 1.0, ps; jac = true, sparse = false)
        oprob_sjac = ODEProblem(rn, u0, 1.0, ps; jac = true, sparse = true)
        sprob_jac = SDEProblem(rn, u0, 1.0, ps; jac = true, sparse = false)
        sprob_sjac = SDEProblem(rn, u0, 1.0, ps; jac = true, sparse = true)
        nlprob_jac = NonlinearProblem(rn, u0, ps; jac = true, sparse = false)
        nlprob_sjac = NonlinearProblem(rn, u0, ps; jac = true, sparse = true)

        # Checks that Jacobians ar identical.
        # Approx is due to https://github.com/SciML/ModelingToolkit.jl/issues/3554.
        function eval_jac(prob, sparse)
            J = sparse ? deepcopy(prob.f.jac_prototype) : zeros(length(prob.u0), length(prob.u0))
            ModelingToolkit.is_time_dependent(prob) ? prob.f.jac(J, prob.u0, prob.p, 0.0) : prob.f.jac(J, prob.u0, prob.p)
            return J
        end
        @test eval_jac(oprob_jac, false) == eval_jac(sprob_jac, false) == eval_jac(nlprob_jac, false)
        @test eval_jac(oprob_sjac, true) ≈ eval_jac(sprob_sjac, true) atol = 1e-14 rtol = 1e-14
        @test eval_jac(oprob_sjac, true) ≈ eval_jac(nlprob_sjac, true) atol = 1e-14 rtol = 1e-14
    end
end

### Sparse Jacobian Tests ###

# Checks that generated dense/sparse Jacobians are identical.
let
    # Creates model (vaguely messy model without conserved quantities).
    # Model includes a time-dependent reaction.
    rn = @reaction_network begin
        (p,d), 0 <--> (X,Y,Z)
        k1, X + Y --> XY
        k2, X + 2Z --> XZ2
        k3, Y3 +X2 --> Y3Z2
        k4, X + Y + Z --> XYZ
        k5, XZ2 + Y3Z2 --> XY3Z4
        k6, XYZ + XYZ --> X2Y2Z2
        d, (XY3Z4, X2Y2Z2) --> 0
        X + Y, V --> 0
        k7/(1 + t), 2V --> V2
        Z, V2 --> 0
    end

    # Performs tests for different randomised values (to be thoroughly sure).
    for factor in [0.1, 1.0, 10.0]
        # Creates randomised species and parameter values. Generates jacobians (dense/sparse).
        u0 = rnd_u0(rn, rng; factor)
        t_val = factor*rand()
        ps = rnd_ps(rn, rng; factor)
        jac = jac_eval(rn, u0, ps, t_val; sparse = false)
        jac_sparse = jac_eval(rn, u0, ps, t_val; sparse = true)

        # Check correctness (both by converting to sparse jac to dense, and through multiplication with other matrix).
        # Approx is due to https://github.com/SciML/ModelingToolkit.jl/issues/3554.
        @test Matrix(jac_sparse) ≈ jac atol = 1e-14 rtol = 1e-14
        mat = factor*rand(rng, length(u0), length(u0))
        @test jac_sparse * mat ≈ jac * mat
    end
end

# Tests that simulations with different Jacobian and sparsity options are identical.
let
    # Creates model (vaguely messy model without conserved quantities).
    rn = @reaction_network begin
        (v0 + mm(X,v,K),d), 0 <--> X + 2Y
        (k1,k2), X + Y <--> XY
        (k1,k2), X + Y2 <--> XY2
        (k3,k4), XY + XY2 <--> X2Y3
        1.0, (XY,XY2,X2Y3) --> 0
        mm(X2Y3,v,K), 0 --> Z
        (k3*X,k4*Y), 3Z <--> Z3
        d, Z --> 0
    end

    # Generates initial conditions and parameter values. Creates problems with/o (sparse/dense) jacobian.
    u0 = rnd_u0(rn, rng)
    ps = rnd_ps(rn, rng)
    oprob = ODEProblem(rn, u0, 1.0, ps; jac = false, sparse = false)
    oprob_j = ODEProblem(rn, u0, 1.0, ps; jac = true, sparse = false)
    oprob_s = ODEProblem(rn, u0, 1.0, ps; jac = false, sparse = true)
    oprob_js = ODEProblem(rn, u0, 1.0, ps; jac = true, sparse = true)

    # Simulates system with implicit solver. Checks that all solutions are identical.
    sol = solve(oprob, Rosenbrock23(), saveat = 0.1, abstol = 1e-8, reltol = 1e-8)
    sol_j = solve(oprob_j, Rosenbrock23(), saveat = 0.1, abstol = 1e-8, reltol = 1e-8)
    sol_s = solve(oprob_s, Rosenbrock23(), saveat = 0.1, abstol = 1e-8, reltol = 1e-8)
    sol_js = solve(oprob_js, Rosenbrock23(), saveat = 0.1, abstol = 1e-8, reltol = 1e-8)
    @test  sol ≈ sol_j ≈ sol_s ≈ sol_js
end
