### Fetch Packages ###

# Fetch packages.
using Catalyst, OrdinaryDiffEq, SteadyStateDiffEq, Test
import HomotopyContinuation

# Sets rnd number.
using StableRNGs
rng = StableRNG(12345)

### Basic Tests ###

# Tests that stability is correctly assessed (using simulation) in multi stable system.
# Tests that `steady_state_jac` function works.
# Tests using symbolic input.
let
    # System which may have between 1 and 7 fixed points.
    rn = @reaction_network begin
        v/20.0 + hillar(X,Y,v,K,n), 0 --> X
        v/20.0 + hillar(Y,X,v,K,n), 0 --> Y
        d, (X,Y) --> 0
    end
    ss_jac = steady_state_jac(rn)

    # Repeats several times, most steady state stability cases should be encountered several times.
    for repeat = 1:20
        # Generates random parameter values (which can generate all steady states cases).
        ps = (:v => 1.0 + 3*rand(rng), :K => 0.5 + 2*rand(rng), :n => rand(rng,[1,2,3,4]), 
              :d => 0.5 + rand(rng))

        # Computes stability using various jacobian options.
        sss = hc_steady_states(rn, ps; show_progress = false)
        stabs_1 = [steady_state_stability(ss, rn, ps) for ss in sss]
        stabs_2 = [steady_state_stability(ss, rn, ps; ss_jac = ss_jac) for ss in sss]

        # Confirms stability using simulations.
        for (idx, ss) in enumerate(sss)
            ssprob = SteadyStateProblem(rn, [1.001, 0.999] .* ss, ps)
            sol = solve(ssprob, DynamicSS(Vern7()); abstol = 1e-8, reltol = 1e-8)
            stabs_3 = isapprox(ss, sol.u; atol = 1e-6)
            @test stabs_1[idx] == stabs_2[idx] == stabs_3

            # Checks stability when steady state is given on a pair form ([X => x_val, Y => y_val]).
            stabs_4 = steady_state_stability(Pair.(unknowns(rn), ss), rn, ps)
            @test stabs_3 == stabs_4
        end
    end
end

# Checks stability for system with known stability structure.
# Tests for system with conservation laws.
# Tests for various input forms of u0 and ps.
let
    # Creates model.
    rn = @reaction_network begin
        k1+Z, Y --> 2X
        k2, 2X --> X + Y
        k3, X + Y --> Y
        k4, X --> 0
        (kD1+X, kD2), 2Z <--> Z2
    end

    # Creates various forms of input.
    @unpack k1, k2, k3, k4, kD1, kD2, X, Y, Z, Z2 = rn
    u0_1 = [X => 1.0, Y => 1.0, Z => 1.0, Z2 => 1.0]
    u0_2 = [:X => 1.0, :Y => 1.0, :Z => 1.0, :Z2 => 1.0]
    u0_3 = [rn.X => 1.0, rn.Y => 1.0, rn.Z => 1.0, rn.Z2 => 1.0]
    u0_4 = [1.0, 1.0, 1.0, 1.0]
    ps_1 = [k1 => 8.0, k2 => 2.0, k3 => 1.0, k4 => 1.5, kD1 => 0.5, kD2 => 2.0]
    ps_2 = [:k1 => 8.0, :k2 => 2.0, :k3 => 1.0, :k4 => 1.5, :kD1 => 0.5, :kD2 => 2.0]
    ps_3 = [rn.k1 => 8.0, rn.k2 => 2.0, rn.k3 => 1.0, rn.k4 => 1.5, rn.kD1 => 0.5, rn.kD2 => 4.0]
    
    # Computes stability using various input forms, and checks that the output is correct. 
    sss = hc_steady_states(rn, ps_1; u0 = u0_1, show_progress = false)
    for u0 in [u0_1, u0_2, u0_3, u0_4], ps in [ps_1, ps_2, ps_3]
        stab_1 =  [steady_state_stability(ss, rn, ps) for ss in sss]
        ss_jac = steady_state_jac(rn; u0 = u0)
        stab_2 =  [steady_state_stability(ss, rn, ps; ss_jac = ss_jac) for ss in sss]
        @test length(stab_1) == length(stab_2) == 3
        @test count(stab_1) == count(stab_2) == 2
    end

    # Confirms error when computing Jacobian with wrong length of u0.
    @test_throws Exception steady_state_jac(rn; u0 = [1.0, 1.0])
end

### Other Tests ###

# Tests `tol` option. In this case, the maximum eigenvalue real part is -1.0, which generates
# and error with `tol = 100`.
let
    rn = @reaction_network begin
        (p,d), 0 <--> X
    end
    p = [:p => 1.0, :d => 1.0]
    u = [1.0]
    @test_throws Exception steady_state_stability(u, rn, p; tol = 1e2)
end
