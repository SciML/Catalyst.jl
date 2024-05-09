### Prepares Tests ###

# Fetch packages.
using BifurcationKit, Catalyst, Test

# Sets stable rng number.
using StableRNGs
rng = StableRNG(12345)


### Basic Tests ###

# Brusselator extended with conserved species.
# Runs full computation, checks values corresponds to known values.
# Checks that the correct bifurcation point is found at the correct position.
# Checks that bifurcation diagrams can be computed for systems with conservation laws.
# Checks that bifurcation diagrams can be computed for systems with default values.
# Checks that bifurcation diagrams can be computed for systems with non-constant rate.
# Checks that not providing conserved species throws and appropriate error.
let 
    # Create model.
    extended_brusselator = @reaction_network begin
        @species W(t) = 2.0
        @parameters k2 = 0.5
        A, ∅ → X
        1, 2X + Y → 3X
        B, X → Y
        1, X → ∅
        (k1*Y, k2), V <--> W
    end
    @unpack A, B, k1 = extended_brusselator
    u0_guess = [:X => 1.0, :Y => 1.0, :V => 0.0, :W => 0.0]
    p_start = [A => 1.0, B => 4.0, k1 => 0.1]
    
    # Computes bifurcation diagram.
    bprob = BifurcationProblem(extended_brusselator, u0_guess, p_start, :B; plot_var = :V, u0 = [:V => 1.0])
    p_span = (0.1, 6.0)
    opts_br = ContinuationPar(dsmin = 0.0001, dsmax = 0.001, ds = 0.0001, max_steps = 10000, p_min = p_span[1], p_max = p_span[2], n_inversion = 4)
    bif_dia = bifurcationdiagram(bprob, PALC(), 2, (args...) -> opts_br; bothside = true)
    
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
    @unpack X, v, K, n, d = bistable_switch
    u0_guess = [X => 1.0]
    p_start = [v => 5.0, K => 2.5, n => 3, d => 1.0]
    bprob = BifurcationProblem(bistable_switch, u0_guess, p_start, K; jac=false, plot_var=X)
    
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

# Creates a system where rn is composed of 4, somewhat nested, networks.
# Tests with defaults within nested networks.
let
    rn1 = @network_component rn1 begin
        @parameters p=1.0
        (p, d), 0 <--> X
    end
    rn2 = @network_component rn2 begin
        @parameters p=2.0
        (p, d), 0 <--> X
    end
    rn3 = @network_component rn3 begin
        @parameters p=3.0
        (p, d), 0 <--> X
    end
    rn4 = @network_component rn4 begin
        @parameters p=4.0
        (p, d), 0 <--> X
    end
    @named rn3 = compose(rn3, [rn4])
    @named rn = compose(rn1, [rn2, rn3])
    rn = complete(rn)

    # Declares parameter values and initial u guess.
    @unpack X, d = rn
    u0_guess = [X => 1.0, rn2.X => 1.0, rn3.X => 1.0, rn3.rn4.X => 1.0]
    p_start = [d => 1.0, rn2.d => 1.0, rn3.d => 1.0, rn3.rn4.d => 1.0]


    # Computes bifurcation diagrams and check them (the final point at [p_value = 6] should be =[p/d]).
    p_span = (0.1, 6.0)

    let
        # Checks top layer.
        bprob = BifurcationProblem(rn, u0_guess, p_start, d; plot_var = X)
        opts_br = ContinuationPar(dsmin = 0.0001, dsmax = 0.001, ds = 0.0001, max_steps = 10000, p_min = p_span[1], p_max = p_span[2], n_inversion = 4)
        bif_dia = bifurcationdiagram(bprob, PALC(), 2, (args...) -> opts_br; bothside = true)
        @test bif_dia.γ.branch[end].x ≈ 1.0/6

        # Checks second layer (1).
        bprob = BifurcationProblem(rn, u0_guess, p_start, rn2.d; plot_var = rn2.X)
        opts_br = ContinuationPar(dsmin = 0.0001, dsmax = 0.001, ds = 0.0001, max_steps = 10000, p_min = p_span[1], p_max = p_span[2], n_inversion = 4)
        bif_dia = bifurcationdiagram(bprob, PALC(), 2, (args...) -> opts_br; bothside = true)
        @test bif_dia.γ.branch[end].x ≈ 2.0/6

        # Checks second layer (2).
        bprob = BifurcationProblem(rn, u0_guess, p_start, rn3.d; plot_var = rn3.X)
        opts_br = ContinuationPar(dsmin = 0.0001, dsmax = 0.001, ds = 0.0001, max_steps = 10000, p_min = p_span[1], p_max = p_span[2], n_inversion = 4)
        bif_dia = bifurcationdiagram(bprob, PALC(), 2, (args...) -> opts_br; bothside = true)
        @test bif_dia.γ.branch[end].x ≈ 3.0/6

        # Checks third layer.
        bprob = BifurcationProblem(rn, u0_guess, p_start, rn3.rn4.d; plot_var = rn3.rn4.X)
        opts_br = ContinuationPar(dsmin = 0.0001, dsmax = 0.001, ds = 0.0001, max_steps = 10000, p_min = p_span[1], p_max = p_span[2], n_inversion = 4)
        bif_dia = bifurcationdiagram(bprob, PALC(), 2, (args...) -> opts_br; bothside = true)
        @test bif_dia.γ.branch[end].x ≈ 4.0/6
    end
end

# Tests for nested model with conservation laws.
let
    # Creates model.
    rn1 = @network_component rn1 begin
        (k1, k2), X1 <--> X2
    end
    rn2 = @network_component rn2 begin
        (l1, l2), Y1 <--> Y2
    end
    @named rn = compose(rn1, [rn2])
    rn = complete(rn)

    # Creates input parameter and species vectors.
    @unpack X1, X2, k1, k2 = rn1
    u_guess = [X1 => 1.0, X2 => 1.0, rn2.Y1 => 1.0, rn2.Y2 => 1.0]
    p_start = [k1 => 1.0, k2 => 1.0, rn2.l1 => 1.0, rn2.l2 => 1.0]
    u0 = [X1 => 1.0, X2 => 0.0, rn2.Y1 => 1.0, rn2.Y2 => 0.0]

    # Computes bifurcation diagram.
    p_span = (0.2, 5.0)
    bprob = BifurcationProblem(rn, u_guess, p_start, k1; plot_var = X1, u0=u0)
    opts_br = ContinuationPar(dsmin = 0.0001, dsmax = 0.001, ds = 0.0001, max_steps = 10000, p_min = p_span[1], p_max = p_span[2], n_inversion = 4)
    bif_dia = bifurcationdiagram(bprob, PALC(), 2, (args...) -> opts_br; bothside = true)

    # Checks that the bifurcation diagram is correct.
    xs = getfield.(bif_dia.γ.branch, :x)
    k1s = getfield.(bif_dia.γ.branch, :param)
    @test all(1 ./ k1s .* (1 .- xs) .≈ xs)

    # Checks that there is an error if information for conserved quantities computation is not provided.
    @test_throws Exception bprob = BifurcationProblem(rn, u_guess, p_start, k1; plot_var = X1)
end


### Other Tests ###

# Checks that bifurcation diagrams can be computed for coupled CRN/DAE systems.
let
    # Prepares the model (production/degradation of X, with equations for volume and X concentration).
    rs = @reaction_network begin
        @parameters k
        @variables C(t)
        @equations begin
            D(V) ~ k*X - V
            C ~ X/V
        end
        (p/V,d/V), 0 <--> X
    end
    u0_guess = [:X => 1.0, :V => 1.0, :C => 1.0]
    p_start = [:p => 2.0, :d => 1.0, :k => 5.0]

    # Computes bifurcation diagram.
    bprob = BifurcationProblem(rs, u0_guess, p_start, :p; plot_var = :C)
    p_span = (0.1, 6.0)
    opts_br = ContinuationPar(dsmin = 0.0001, dsmax = 0.001, ds = 0.0001, max_steps = 10000, p_min = p_span[1], p_max = p_span[2], n_inversion = 4)
    bif_dia = bifurcationdiagram(bprob, PALC(), 2, (args...) -> opts_br; bothside = true)
    @test all(getfield.(bif_dia.γ.branch, :x) .≈ 0.2)
end

# Checks that `BifurcationProblem`s cannot be generated from non-complete `ReactionSystems`s.
let 
    # Create model.
    incomplete_network = @network_component begin
        (p, d), 0 <--> X
    end
    u0_guess = [:X => 1.0]
    p_start = [:p => 1.0, :d => 0.2]
    
    # Computes bifurcation diagram.
    @test_throws Exception BifurcationProblem(incomplete_network, u0_guess, p_start, :p)
end