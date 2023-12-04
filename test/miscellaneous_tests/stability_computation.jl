### Fetch Packages ###

# Fetch packages.
using Catalyst, OrdinaryDiffEq
using Random, Test
import HomotopyContinuation

# Sets rnd number.
using StableRNGs
rng = StableRNG(12345)

### Run Tests ###

# Tests that stability is correctly assessed (using simulation) in multi stable system.
# Tests that `steady_state_jac` function works.
# Tests with and without sparsity.
# tests using symbolic input.
let
    # System which may have between 1 and 7 fixed points.
    rn = @reaction_network begin
        v/20.0 + hillar(X,Y,v,K,n), 0 --> X
        v/20.0 + hillar(Y,X,v,K,n), 0 --> Y
        d, (X,Y) --> 0
    end
    ss_jac = steady_state_jac(rn)
    ss_jac_sparse = steady_state_jac(rn; sparse=true)

    # Repeats several times, most cases should be encountered several times.
    for i = 1:50
        # Generates random parameter values (which can generate all steady states cases).
        p = [:v => 1.0 + 3*rand(rng), :K => 0.5 + 2*rand(rng), :n => rand(rng,[1,2,3,4]), :d => 0.5 + rand(rng)]

        # Computes stability using various jacobian options.
        sss = hc_steady_states(rn, p)
        stabs_1 = steady_state_stability(sss, rn, p)
        stabs_2 = steady_state_stability(sss, rn, p; sparse=true)
        stabs_3 = steady_state_stability(sss, rn, p; ss_jac=ss_jac)
        stabs_4 = steady_state_stability(sss, rn, p; ss_jac=ss_jac_sparse)

        # Confirms stability using simulations.
        for (idx,ss) in enumerate(sss)
            oprob = ODEProblem(rn, [1.001, 0.999] .* ss, (0.0,1000.0), p)
            sol_end = solve(oprob, Rosenbrock23())[end]
            stabs_5 = ss ≈ sol_end
            @test stabs_1[idx] == stabs_2[idx] == stabs_3[idx] == stabs_4[idx] == stabs_5

            # Checks stability when steady state is given on a pair form ([:X => x_val, :Y => y_val]).
            stabs_6 = steady_state_stability(Pair.(states(rn),ss), rn, p)
            @test stabs_5 == stabs_6
        end
    end
end

# Checks stability for system with known stability structure.
# Tests for system with conservation laws.
# Tests for various input forms of u0 and ps.
let
    # Creates model.
    rn = complete(@reaction_network begin
        k1+Z, Y --> 2X
        k2, 2X --> X + Y
        k3, X + Y --> Y
        k4, X --> 0
        (kD1+X, kD2), 2Z <--> Z2
    end)

    # Creates various forms of input.
    @unpack k1, k2, k3, k4, kD1, kD2, X, Y, Z, Z2 = rn
    u0_1 = [X => 1.0, Y => 1.0, Z => 1.0, Z2 => 1.0]
    u0_2 = [:X => 1.0, :Y => 1.0, :Z => 1.0, :Z2 => 1.0]
    u0_3 = [rn.X => 1.0, rn.Y => 1.0, rn.Z => 1.0, rn.Z2 => 1.0]
    u0_4 = [1.0, 1.0, 1.0, 1.0]
    ps_1 = [k1 => 8.0, k2 => 2.0, k3 => 1.0, k4 => 1.5, kD1 => 0.5, kD2 => 2.0]
    ps_2 = [:k1 => 8.0, :k2 => 2.0, :k3 => 1.0, :k4 => 1.5, :kD1 => 0.5, :kD2 => 2.0]
    ps_3 = [rn.k1 => 8.0, rn.k2 => 2.0, rn.k3 => 1.0, rn.k4 => 1.5, rn.kD1 => 0.5, rn.kD2 => 4.0]
    ps_4 = [8.0, 2.0, 1.0, 1.5, 0.5, 4.0]
    
    # Computes stability using various input forms, and checks that the output is correct. 
    sss = hc_steady_states(rn, ps_1; u0=u0_1)
    for u0 in [u0_1, u0_2, u0_3, u0_4], ps in [ps_1, ps_2, ps_3, ps_4]
        stab_1 = steady_state_stability(sss, rn, ps)
        @test length(stab_1) == 3
        @test count(stab_1) == 2
    
        ss_jac = steady_state_jac(rn; u0=u0)
        stab_2 = steady_state_stability(sss, rn, ps; ss_jac=ss_jac)
        @test length(stab_2) == 3
        @test count(stab_2) == 2
    end

    # Confirms error when computing Jacobian with wrong length of u0.
    @test_throws Exception steady_state_jac(rn; u0=[1.0, 1.0])
end