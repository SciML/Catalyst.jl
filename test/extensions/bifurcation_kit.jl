### Fetch Packages ###
using Bifurcationkit, Catalyst, Test

# Sets rnd number.
using StableRNGs
rng = StableRNG(12345)

### Run Tests ###

# Brusselator extended with conserved species.
# Runs full computation, checks values corresponds to known values.
# Checks that teh correct bifurcation point is found at the correct position.
# Checks that bifurcation diagrams can be computed for systems with conservation laws.
# Checks that bifurcation diagrams can be computed for systems with default values.
# Checks that bifurcation diagrams can be computed for systems with non-constant rate.
# Checks that not providing conserved species throws and appropriate error.
let 
    # Create model
    extended_brusselator = @reaction_network begin
        @species W(t) = 2.0
        @parameters k2 = 0.5
        A, ∅ → X
        1, 2X + Y → 3X
        B, X → Y
        1, X → ∅
        (k1*Y, k2), V <--> W
    end
    @unpack A, B, k1
    u0_guess = [:X => 1.0, :Y => 1.0, :V => 0.0, :W => 0.0]
    p_start = [A => 1.0, B => 4.0, k1 => 0.1]
    
    # Computes bifurcation diagram.
    BifurcationProblem(extended_brusselator, u0_guess, p_start, :B; plot_var=:V, u0 = [:V => 1.0])
    p_span = (0.1, 6.0)
    opt_newton = NewtonPar(tol = 1e-9, max_iterations = 100)
    opts_br = ContinuationPar(dsmin = 0.0001, dsmax = 0.001, ds = 0.0001,
        max_steps = 200000, nev = 2, newton_options = opt_newton,
        p_min = p_span[1], p_max = p_span[2],
        detect_bifurcation = 3, n_inversion = 4, tol_bisection_eigenvalue = 1e-8, dsmin_bisection = 1e-9)
    bif_dia = bifurcationdiagram(bprob, PALC(), 2, (args...) -> opts_br; bothside=true)
    
    # Checks computed V values are correct (Formula: V = k2*(V0+W0)/(k1*Y+k2), where Y=2*B.)
    B_vals = getfield.(bif_dia.γ.branch, :param)
    V_vals = getfield.(bif_dia.γ.branch, :x)
    @test all(V_vals .≈ 0.5*(1.0+2.0) ./ (0.1 .* 2*B_vals .+ 0.5))
    
    # Checks that the bifurcation point is correct.
    @test length(bif_dia.γ.specialpoint) == 3 # Includes start and end point.
    hopf_bif_point = filter(sp -> sp.type == :hopf, bif_dia.γ.specialpoint)[1]
    @test isapprox(hopf_bif_point.param, 1.5, atol=1e-5)

    # Tests that an error is thrown if information of conserved species is not fully provided.
    @test_throws Exception BifurcationProblem(extended_brusselator, u0_guess, p_start, :B; plot_var=:V, u0 = [:X => 1.0])
end

# Bistable switch.
# Checks that the same bifurcation problem is created as for BifurcationKit.
# Checks with Symbolics as bifurcation and plot vars.
# Tries setting `jac=false`.
let 
    # Creates BifurcationProblem via Catalyst.
    bistable_switch = @reaction_network begin
        0.1 + hill(X,v,K,n), 0 --> X
        d, X --> 0
    end
    @unpack x, v, K, n, d = rn
    u0_guess = [x => 1.0]
    p_start = [v => 5.0, K => 2.5, n => 3, d => 1.0]
    bprob = BifurcationProblem(bistable_switch, u0_guess, p_start, K; jac=false; plot_var=x)
    
    # Creates BifurcationProblem via BifurcationKit.
    function bistable_switch_BK(u, p)
        X, = u
        v, K, n, d = p
        return [0.1 + v*(X^n)/(X^n + K^n) - d*X]
    end
    bprob_BK = BifurcationProblem(bistable_switch_BK, [1.0], [5.0, 2.5, 3, 1.0], (@lens _[1]); record_from_solution = (x, p) -> x[1])
    
    # Check the same function have been generated.
    bprob.u0 == bprob_BK.u0
    bprob.params == bprob_BK.params
    for repeat = 1:20
        u0 = rand(rng, 1)
        p = rand(rng, 4)
        @test bprob_BK.VF.F(u0, p) == bprob.VF.F(u0, p)
    end
end
