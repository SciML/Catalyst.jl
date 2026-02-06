# Fetch packages.
using Catalyst, NonlinearSolve, OrdinaryDiffEqVerner, OrdinaryDiffEqTsit5, OrdinaryDiffEqRosenbrock, Statistics, SteadyStateDiffEq, StochasticDiffEq, Test
using ModelingToolkitBase: getdefault, getdescription, getdefault
using Symbolics: unwrap

# Sets stable rng number.
using StableRNGs
rng = StableRNG(12345)
seed = rand(rng, 1:100)

# Sets the default `t` and `D` to use.
t = default_t()
D = default_time_deriv()


### Basic Coupled Differential Equations Tests ###

# Tests coupled CRN/ODE. Checks that known steady state is reached.
# Check that steady state can be found using NonlinearSolve and SteadyStateDiffEq.
let
    # Creates coupled reactions system.
    @parameters p d k
    @species X(t)
    @variables A(t)
    eqs = [
        Reaction(p, nothing, [X]),
        Reaction(d, [X], nothing),
        D(A) ~ p*X - k*A
    ]
    @named coupled_rs = ReactionSystem(eqs, t)
    coupled_rs = complete(coupled_rs)

    # Basic model checks.
    @test issetequal(parameters(coupled_rs), [p, d, k])
    @test issetequal(species(coupled_rs), unknowns(coupled_rs)[1:1])
    @test issetequal(unknowns(coupled_rs)[1:1], [X])
    @test issetequal(unknowns(coupled_rs)[2:2], [A])
    @test issetequal(reactions(coupled_rs), equations(coupled_rs)[1:2])
    @test issetequal(equations(coupled_rs)[1:2], eqs[1:2])
    @test issetequal(equations(coupled_rs)[3:3], eqs[3:3])

    # Set simulation inputs.
    u0 = [X => 0.1, A => 10.0]
    tspan = (0.0, 1000.0)
    ps = [p => 1.0, d => 0.5, k => 2.0]

    # Checks that the correct steady state is found through ODEProblem.
    oprob = ODEProblem(coupled_rs, u0, tspan, ps)
    osol = solve(oprob, Vern7(); abstol = 1e-8, reltol = 1e-8)
    @test osol[[X,A]][end] ≈ [2.0, 1.0]

    # Checks that the correct steady state is found through NonlinearProblem.
    nlprob = NonlinearProblem(coupled_rs, u0, ps)
    nlsol = solve(nlprob; abstol = 1e-8, reltol = 1e-8)
    @test nlsol[[X,A]] ≈ [2.0, 1.0]

    # Checks that the correct steady state is found through SteadyStateProblem.
    ssprob = SteadyStateProblem(coupled_rs, u0, ps)
    sssol = solve(ssprob, DynamicSS(Rosenbrock23()); abstol = 1e-8, reltol = 1e-8)
    @test sssol[[X,A]] ≈ [2.0, 1.0]
end

# Checks that coupled systems created via the DSL, extension, and programmatically are identical.
# Check that these contain the correct stuff (in the correct order).
# Checks that that these systems yield identical simulations reaching the known (correct) steady state.
# Checks interpolation of variables between the reaction system and ODE.
# Checks that reactions and equations are reordered, even if given in the wrong order to programmatic
# model creation.
let
    # Creates the model fully programmatically
    @parameters k1 k2 a b
    @species X1(t) X2(t)
    @variables A(t) B(t)
    eqs_prog = [
        D(A) ~ X1 + a - A,
        D(B) ~ X2 + b - B,
        Reaction(k1*A, [X1], [X2]),
        Reaction(k2*B, [X2], [X1])
    ]
    coupled_rs_prog = ReactionSystem(eqs_prog, t, [A, B, X1, X2], [k1, k2, a, b]; name = :coupled_rs)
    coupled_rs_prog = complete(coupled_rs_prog)

    # Creates model by extending a `ReactionSystem` with another ReactionSystem containing ODEs.
    rn_extended = @network_component begin
        ($k1*$A, $k2*$B), X1 <--> X2
    end
    eqs_extended = [
        D(A) ~ X1 + a - A
        D(B) ~ X2 + b - B
    ]
    @named rs_extended = ReactionSystem(eqs_extended, t)
    coupled_rs_extended = complete(extend(rs_extended, rn_extended; name = :coupled_rs))

    # Creates the model through the DSL.
    coupled_rs_dsl = @reaction_network coupled_rs begin
        @parameters k1 k2 a b
        @equations begin
            D(A) ~ X1 + a - A
            D(B) ~ X2 + b - B
        end
        (k1*A, k2*B), X1 <--> X2
    end

    # Checks that models are equivalent and contain the correct stuff.
    @test Catalyst.isequivalent(coupled_rs_prog, coupled_rs_extended)
    @test Catalyst.isequivalent(coupled_rs_extended, coupled_rs_dsl)
    @test issetequal(parameters(coupled_rs_extended), [a, b, k1, k2])
    @test issetequal(species(coupled_rs_extended), [X1, X2])
    @test issetequal(unknowns(coupled_rs_extended)[1:2], [X1, X2])
    @test issetequal(unknowns(coupled_rs_extended)[3:4], [A, B])
    @test issetequal(equations(coupled_rs_extended)[3:4], eqs_extended)

    # Simulates the three models, checking that they all yield the correct end point.
    u0 = [A => 1.0, B => 1.0, X1 => 10.0, X2 => 10.0]
    tspan = (0.0, 100.)
    ps = [a => 1.0, b => 1.0, k1 => 1.0, k2 => 1.0]
    for coupled_rs in [coupled_rs_prog, coupled_rs_extended, coupled_rs_dsl]
        oprob = ODEProblem(coupled_rs, u0, tspan, ps)
        osol = solve(oprob, Vern7(); abstol = 1e-8, reltol = 1e-8)
        osol[[A,B,X1,X2]][end] ≈ [10.0, 10.0, 11.0, 11.0]
    end
end


### Basic Coupled Algebraic Equations Tests ###

# Tests coupled CRN/algebraic equation. Checks that known steady state is reached using ODE solve.
# Check that steady state can be found using NonlinearSolve and SteadyStateDiffEq.
# Checks that errors are given if `structural_simplify = true` argument is not given.
let
    # Creates a simple coupled model with an algebraic equation.
    @parameters p d a b
    @species X(t)
    @variables A(t)
    eqs = [
        Reaction(p, nothing, [X]),
        Reaction(d, [X], nothing),
        a*A^2 ~ X + b
    ]
    @named coupled_rs = ReactionSystem(eqs, t)
    coupled_rs = complete(coupled_rs)

    # Check model content.
    @test issetequal(parameters(coupled_rs), [p, d, a, b])
    @test issetequal(species(coupled_rs), unknowns(coupled_rs)[1:1])
    @test issetequal(unknowns(coupled_rs)[1:1], [X])
    @test issetequal(unknowns(coupled_rs)[2:2], [A])
    @test issetequal(reactions(coupled_rs), equations(coupled_rs)[1:2])
    @test issetequal(equations(coupled_rs)[1:2], eqs[1:2])
    @test issetequal(equations(coupled_rs)[3:3], eqs[3:3])

    # Set simulation inputs.
    u0 = [X => 0.1]
    tspan = (0.0, 1000.0)
    ps = [p => 1.0, d => 0.5, a => 2.0, b => 16.0]

    # Checks not using `structural_simplify` argument yields an error.
    @test_throws Exception ODEProblem(coupled_rs, u0, tspan, ps)
    @test_throws Exception SteadyStateProblem(coupled_rs, u0, ps)

    # Checks that the correct steady state is found through ODEProblem.
    oprob = ODEProblem(coupled_rs, u0, tspan, ps; structural_simplify = true,
        guesses = [A => 1.0])
    osol = solve(oprob, Rosenbrock23(); abstol = 1e-8, reltol = 1e-8)
    @test osol[[X,A]][end] ≈ [2.0, 3.0]

    # Checks that the correct steady state is found through NonlinearProblem.
    u0 = [X => 0.1, A => 16.1/2]
    nlprob = NonlinearProblem(coupled_rs, u0, ps, structural_simplify = true)
    nlsol = solve(nlprob)
    @test nlsol[[X,A]] ≈ [2.0, 3.0]

    # Checks that the correct steady state is found through SteadyStateProblem.
    u0 = [X => 0.1, A => 1.0]
    ssprob = SteadyStateProblem(coupled_rs, u0, ps; structural_simplify = true)
    sssol = solve(ssprob, DynamicSS(Rosenbrock23()); abstol = 1e-8, reltol = 1e-8)
    @test_broken sssol[[X,A]] ≈ [2.0, 3.0] # The previous lines fails to solve. Issue at: https://github.com/SciML/ModelingToolkit.jl/issues/4174. This currently also yields a warning.
end


### Basic Combined Coupled Algebraic/Differential Equations Tests ###

# Checks that a combined reaction/differential/algebraic coupled system can be created.
# Checks that it can its ODE, SteadyState, and Nonlinear problems all can be solved.
# Checks that Tuple u0/ps input, and non-default independent variables works.
# The system is mostly made up to be non-trivial, but reliably solvable.
let
    @parameters p d a b c
    @parameters τ
    @variables A(τ) B(τ) C(τ)
    @species X(τ)
    Δ = Differential(τ)
    eqs = [
        Δ(A) ~ b + X - A,
        Δ(B) ~ sqrt(A + X + b) - B,
        Reaction(p, nothing, [X], nothing, [2]),
        Reaction(d, [X], nothing),
        (X + C)*B ~ A
    ]
    @named coupled_rs = ReactionSystem(eqs, τ)
    coupled_rs = complete(coupled_rs)

    # Set simulation inputs.
    u0 = (X => 2.0, A => 4.0, B => 1.0, C => 2.0)
    ps = (p => 1.0, d => 2.0, b => 4.0)

    # Creates and solves a ODE, SteadyState, and Nonlinear problems.
    # Success is tested by checking that the same steady state solution is found.
    oprob = ODEProblem(coupled_rs, u0, (0.0, 1000.0), ps; structural_simplify = true,
        warn_initialize_determined = false)
    ssprob = SteadyStateProblem(coupled_rs, u0, ps; structural_simplify = true,
        warn_initialize_determined = false)
    nlprob = NonlinearProblem(coupled_rs, u0, ps; structural_simplify = true)
    osol = solve(oprob, Rosenbrock23(); abstol = 1e-8, reltol = 1e-8)
    sssol = solve(ssprob, DynamicSS(Rosenbrock23()); abstol = 1e-8, reltol = 1e-8)
    nlsol = solve(nlprob; abstol = 1e-8, reltol = 1e-8)
    @test osol[[A, B, C, X]][end] ≈ sssol[[A, B, C, X]] ≈ nlsol[[A, B, C, X]]
end


### Accessor Tests ###

# Checks basic accessor functions for a basic coupled CRN/equation model.
let
    # Creates a reaction system.
    t = default_t()
    D = default_time_deriv()
    @parameters p d v
    @species X(t)
    @variables V(t) W(t)

    eqs = [
        Reaction(p, [], [X]),
        Reaction(d, [X], []),
        Reaction(d, [X], nothing, [2], nothing),
        D(V) ~ X - v*V,
        W^2 ~ log(V) + X
    ]
    @named coupled_rs = ReactionSystem(eqs, t)

    # Check unknowns-related accessors.
    @test Catalyst.has_species(coupled_rs)
    @test issetequal(Catalyst.get_species(coupled_rs), [X])
    @test issetequal(species(coupled_rs), [X])
    @test issetequal(ModelingToolkitBase.get_unknowns(coupled_rs), [X, V, W])
    @test issetequal(unknowns(coupled_rs), [X, V, W])
    @test issetequal(nonspecies(coupled_rs), [V, W])
    @test numspecies(coupled_rs) == 1

    # Check parameters-related accessors.
    @test Catalyst.has_rxs(coupled_rs)
    @test issetequal(Catalyst.get_rxs(coupled_rs), eqs[1:3])
    @test issetequal(reactions(coupled_rs), eqs[1:3])
    @test issetequal(equations(coupled_rs), eqs)
    @test issetequal(nonreactions(coupled_rs), eqs[4:5])
    @test issetequal(reactionrates(coupled_rs), [p, d, d])
    @test numreactions(coupled_rs) == 3

    # Check parameters-related accessors.
    @test issetequal(parameters(coupled_rs), [p, d, v])
    @test numparams(coupled_rs) == 3

    # Check other accessors.
    @test !isspatial(coupled_rs)
end


### Species, Variables, and Parameter Handling ###

# Checks that coupled systems contain the correct species, variables, and parameters.
# Checks that species, variables, and parameters are inferred correctly from equations.
# Checks that non-default iv is inferred correctly from reactions/equations.
let
    # Create coupled model.
    @parameters τ
    @variables A(τ) B(τ)
    @species X(τ) X2(τ)
    @parameters k1 k2 k b1 b2
    D = Differential(τ)
    eqs = [
        Reaction(k1, [X], [X2], [2], [1]),
        Reaction(k2, [X2], [X], [1], [2]),
        D(A) ~ k*X2 - A,
        B + A ~ b1*X + b2*X2
    ]
    @named coupled_rs = ReactionSystem(eqs, τ)
    coupled_rs = complete(coupled_rs)

    # Checks that systems created from coupled reaction systems contain the correct content
    osys = ode_model(coupled_rs)
    ssys = sde_model(coupled_rs)
    nlsys = ss_ode_model(coupled_rs)
    initps = Initial.((X, X2, A, B))
    fullps = union(initps, [k1, k2, k, b1, b2])
    for sys in [coupled_rs, osys, ssys, nlsys]
        @test issetequal(parameters(sys), [k1, k2, k, b1, b2])
        @test issetequal(unknowns(sys), [A, B, X, X2])
    end
end

# Checks that parameters, species, and variables can be correctly accessed in coupled systems.
# Checks for both differential and algebraic equations.
# Checks for problems, integrators, and solutions yielded by coupled systems.
# Checks that metadata, types, and default values are carried through correctly.
let
    # Creates the model
    @parameters a1 [description="Parameter a1"] a2::Rational{Int64} a3=0.3 a4::Rational{Int64}=4//10 [description="Parameter a4"]
    @parameters b1 [description="Parameter b1"] b2::Int64 b3 = 3 b4::Int64=4 [description="Parameter b4"]
    @parameters c1 [description="Parameter c1"] c2::Float32 c3=30.0 c4::Float32=40.0 [description="Parameter c4"]
    @species A1(t) [description="Species A1"] A2(t)=0.2 A3(t)=0.3 [description="Species A3"] A4(t)
    @variables B1(t) [description="Variable B1"] B2(t)=2.0 B3(t)=3.0 [description="Variable B3"] B4(t)
    @variables C1(t) [description="Variable C1"] C2(t) C3(t) [description="Variable C3"] C4(t)
    eqs = [
        Reaction(a1, [A1], nothing),
        Reaction(a2, [A2], nothing),
        Reaction(a3, [A3], nothing),
        Reaction(a4, [A4], nothing),
        D(B1) ~ b1*B1,
        D(B2) ~ b2*B2,
        D(B3) ~ b3*B3,
        D(B4) ~ b4*B4,
        C1 ~ sqrt(c1 + B1^5),
        C2 ~ sqrt(c2 + B2^5),
        C3 ~ sqrt(c3 + B3^5),
        C4 ~ sqrt(c4 + B4^5)
    ]
    @named coupled_rs = ReactionSystem(eqs, t)
    coupled_rs = complete(coupled_rs)

    # Checks that the model has the correct content.
    @test issetequal(parameters(coupled_rs), [a1, a2, a3, a4, b1, b2, b3, b4, c1, c2, c3, c4])
    @test issetequal(species(coupled_rs), unknowns(coupled_rs)[1:4])
    @test issetequal(unknowns(coupled_rs)[1:4], [A1, A2, A3, A4])
    @test issetequal(unknowns(coupled_rs)[5:12], [B1, B2, B3, B4, C1, C2, C3, C4])
    @test issetequal(reactions(coupled_rs)[1:4], equations(coupled_rs)[1:4])
    @test issetequal(equations(coupled_rs)[1:4], eqs[1:4])
    @test issetequal(equations(coupled_rs)[5:12], eqs[5:12])

    # Checks that parameters, species, and variables carried the correct information.
    @test SymbolicUtils.symtype(coupled_rs.a1) == Real
    @test SymbolicUtils.symtype(coupled_rs.a2) == Rational{Int64}
    @test SymbolicUtils.symtype(coupled_rs.a3) == Real
    @test SymbolicUtils.symtype(coupled_rs.a4) == Rational{Int64}
    @test SymbolicUtils.symtype(coupled_rs.b1) == Real
    @test SymbolicUtils.symtype(coupled_rs.b2) == Int64
    @test SymbolicUtils.symtype(coupled_rs.b3) == Real
    @test SymbolicUtils.symtype(coupled_rs.b4) == Int64
    @test SymbolicUtils.symtype(coupled_rs.c1) == Real
    @test SymbolicUtils.symtype(coupled_rs.c2) == Float32
    @test SymbolicUtils.symtype(coupled_rs.c3) == Real
    @test SymbolicUtils.symtype(coupled_rs.c4) == Float32
    @test getdescription(coupled_rs.a1) == "Parameter a1"
    @test getdescription(coupled_rs.a4) == "Parameter a4"
    @test getdescription(coupled_rs.b1) == "Parameter b1"
    @test getdescription(coupled_rs.b4) == "Parameter b4"
    @test getdescription(coupled_rs.c1) == "Parameter c1"
    @test getdescription(coupled_rs.c4) == "Parameter c4"
    @test getdefault(coupled_rs.a3) == 0.3
    @test getdefault(coupled_rs.a4) == 4//10
    @test getdefault(coupled_rs.b3) == 3
    @test getdefault(coupled_rs.b4) == 4
    @test getdefault(coupled_rs.c3) == 30
    @test getdefault(coupled_rs.c4) == 40
    @test getdescription(coupled_rs.A1) == "Species A1"
    @test getdescription(coupled_rs.A3) == "Species A3"
    @test getdescription(coupled_rs.B1) == "Variable B1"
    @test getdescription(coupled_rs.B3) == "Variable B3"
    @test getdescription(coupled_rs.C1) == "Variable C1"
    @test getdescription(coupled_rs.C3) == "Variable C3"
    @test getdefault(coupled_rs.A2) == 0.2
    @test getdefault(coupled_rs.A3) == 0.3
    @test getdefault(coupled_rs.B2) == 2.0
    @test getdefault(coupled_rs.B3) == 3.0

    # Creates problem inputs.
    u0 = [A1 => 0.1, A4 => 0.4, B1 => 1.0, B4 => 4.0]
    u0_nlp = [A1 => 0.1, A4 => 0.4, B1 => 1.0, B4 => 4.0, C1 => sqrt(10.0 + 1.0^5), C2 => sqrt(20.0 + 2.0^5), C3 => sqrt(30.0 + 3.0^5), C4 => sqrt(40.0 + 4.0^5)]
    tspan = (0.0, 1.0)
    ps = [a1 => 0.1, a2 => 2//10, b1 => 1.0, b2 => 2, c1 => 10.0, c2 => 20.0]

    # Create ODE structures.
    oprob = ODEProblem(coupled_rs, u0, tspan, ps; structural_simplify = true, warn_initialize_determined = false)
    oint = init(oprob, Rosenbrock23())
    osol = solve(oprob, Rosenbrock23())

    # Create SDE structures.
    sprob = SDEProblem(coupled_rs, u0, tspan, ps; structural_simplify = true, warn_initialize_determined = false)
    sint = init(sprob, ImplicitEM())
    ssol = solve(sprob, ImplicitEM())

    # Creates Nonlinear structures.
    nlprob = NonlinearProblem(coupled_rs, u0_nlp, ps; structural_simplify = true, warn_initialize_determined = false)
    nlint = init(nlprob, NewtonRaphson())
    nlsol = solve(nlprob, NewtonRaphson())

    # Checks indexing.
    for mtk_struct in [oprob, oint, osol, sprob, sint, ssol, nlprob, nlint, nlsol]
        # Parameters.
        @test mtk_struct.ps[a1] == 0.1
        @test_broken mtk_struct.ps[a2] == 2//10 # An equivalent value (but different nominator and denominator) is generated. https://github.com/SciML/ModelingToolkit.jl/issues/4163.
        @test mtk_struct.ps[a3] == 0.3
        @test mtk_struct.ps[a4] == 4//10
        @test mtk_struct.ps[b1] == 1.0
        @test mtk_struct.ps[b2] == 2
        @test mtk_struct.ps[b3] == 3.0
        @test mtk_struct.ps[b4] == 4
        @test mtk_struct.ps[c1] == 10.0
        @test mtk_struct.ps[c2] == 20.0
        @test mtk_struct.ps[c3] == 30.0
        @test mtk_struct.ps[c4] == 40.0
    end
    for mtk_struct in [oprob, oint, sprob, sint, nlprob, nlint]

        # Species.
        @test mtk_struct[A1] == 0.1
        @test mtk_struct[A2] == 0.2
        @test mtk_struct[A3] == 0.3
        @test mtk_struct[A4] == 0.4

        # Variables.
        @test mtk_struct[B1] == 1.0
        @test mtk_struct[B2] == 2
        @test mtk_struct[B3] == 3.0
        @test mtk_struct[B4] == 4
        @test mtk_struct[C1] == sqrt(mtk_struct.ps[c1] + mtk_struct[B1]^5)
        @test mtk_struct[C2] == sqrt(mtk_struct.ps[c2] + mtk_struct[B2]^5)
        @test mtk_struct[C3] == sqrt(mtk_struct.ps[c3] + mtk_struct[B3]^5)
        @test mtk_struct[C4] == sqrt(mtk_struct.ps[c4] + mtk_struct[B4]^5)
    end
    for mtk_struct in [osol, ssol]

        # Species.
        @test mtk_struct[A1][1] == 0.1
        @test mtk_struct[A2][1] == 0.2
        @test mtk_struct[A3][1] == 0.3
        @test mtk_struct[A4][1] == 0.4

        # Variables.
        @test mtk_struct[B1][1] == 1.0
        @test mtk_struct[B2][1] == 2
        @test mtk_struct[B3][1] == 3.0
        @test mtk_struct[B4][1] == 4
        @test mtk_struct[C1][1] == sqrt(mtk_struct.ps[c1] + mtk_struct[B1][1]^5)
        @test mtk_struct[C2][1] == sqrt(mtk_struct.ps[c2] + mtk_struct[B2][1]^5)
        @test mtk_struct[C3][1] == sqrt(mtk_struct.ps[c3] + mtk_struct[B3][1]^5)
        @test mtk_struct[C4][1] == sqrt(mtk_struct.ps[c4] + mtk_struct[B4][1]^5)
    end
end


### Coupled SDESystem Tests ###

# Checks that a coupled SDE + differential equations works.
# Checks that CLE noise does not affect ODE part that should be deterministic.
# Only considers added differential equations without noise.
let
    # Creates coupled reactions system.
    @parameters p d k1 k2
    @species X(t)
    @variables A(t) B(t)
    eqs = [
        Reaction(p, nothing, [X]; metadata = [:noise_scaling => 0.1]),
        Reaction(d, [X], nothing; metadata = [:noise_scaling => 0.1]),
        D(A) ~ X - k1*A,
        D(B) ~ k2 - B
    ]
    @named coupled_rs = ReactionSystem(eqs, t)
    coupled_rs = complete(coupled_rs)

    # Set simulation inputs.
    u0 = [X => 100.0, A => 50.0, B => 2.0]
    tspan = (0.0, 10000.0)
    ps = [p => 10.0, d => 0.1, k1 => 2.0, k2 => 20.0]

    # Checks that the simulations have the expected means (or endpoint, for B).
    sprob = SDEProblem(coupled_rs, u0, tspan, ps)
    ssol = solve(sprob, ImplicitEM(); maxiters = 1e9, seed)
    @test mean(ssol[:X]) ≈ 100.0 atol = 1e-2 rtol = 1e-2
    @test mean(ssol[:A]) ≈ 50.0 atol = 1e-2 rtol = 1e-2
    @test ssol[:B][end] ≈ 20.0
end

# Checks that a coupled SDE + algebraic equations works.
# Checks that structural_simplify is required to simulate coupled SDE + algebraic equations.
let
    # Creates coupled reactions system.
    @parameters p d k1 k2
    @species X(t)
    @variables A(t)
    eqs = [
        Reaction(p, nothing, [X]),
        Reaction(d, [X], nothing),
        2 + k1 * A ~ 3 + k2 * X
    ]
    @named coupled_rs = ReactionSystem(eqs, t)
    coupled_rs = complete(coupled_rs)

    # Set simulation inputs.
    u0 = [X => 100.0]
    tspan = (0.0, 1000.0)
    ps = Dict([p => 1.0, d => 0.01, k1 => 3.0, k2 => 4.0])

    # Check that the structural_simplify argument is required.
    @test_throws Exception SDEProblem(coupled_rs, u0, tspan, ps)

    # Checks the algebraic equation holds.
    sprob = SDEProblem(coupled_rs, u0, tspan, ps; guesses = [A => 1.0], structural_simplify = true, warn_initialize_determined = false)
    ssol = solve(sprob, ImplicitEM(); abstol = 1e-5, reltol = 1e-5)
    @test (2 .+ ps[k1] * ssol[:A]) ≈ (3 .+ ps[k2] * ssol[:X]) atol = 1e-1 rtol = 1e-1
end


### Coupled NonlinearSystems Tests ###

# Checks that systems with weird differential equations yield errors.
let
    # This one is normal, and should not yield an error.
    begin
        rs = @reaction_network begin
            @equations D(V) ~ 1.0 - V
        end
        @test_nowarn ss_ode_model(rs)
    end

    # Higher-order differential on the lhs, should yield an error.
    begin
        rs = @reaction_network begin
            @differentials D = Differential(t)
            @variables V(t)
            @equations D(D(V)) ~ 1.0 - V
            (p,d), 0 <--> X
        end
        @test_throws Exception ss_ode_model(rs)
    end

    # Differential on the rhs, should yield an error.
    begin
        rs = @reaction_network begin
            @variables U(t)
            @equations D(V) ~ 1.0 - V + D(U)
            (p,d), 0 <--> X
        end
        @test_throws Exception ss_ode_model(rs)
    end

    # Non-differential term on the lhs, should yield an error.
    begin
        rs = @reaction_network begin
            @differentials D = Differential(t)
            @variables V(t)
            @equations D(V) + V ~ 1.0 - V
            (p,d), 0 <--> X
        end
        @test_throws Exception ss_ode_model(rs)
    end
end


### Unusual Differentials Tests ###

# Tests that coupled CRN/DAEs with higher order differentials can be created.
# Tests that these can be solved using ODEs, nonlinear solving, and steady state simulations.
@test_broken let # Cannot generate System, issue in  # https://github.com/SciML/ModelingToolkit.jl/issues/4186.
    # Create coupled model.
    @species X(t)
    @variables A(t) B(t)
    @parameters p d ω k
    eqs = [
        Reaction(p, nothing, [X]),
        Reaction(d, [X], nothing),
        D(D(A)) + 2ω*D(A) +(ω^2)*A ~ 0,
        A + k*(B + D(A)) ~ X
    ]
    @named coupled_rs = ReactionSystem(eqs, t)
    coupled_rs = complete(coupled_rs)
    u0 = [X => 1.0, A => 2.0, D(A) => 1.0]
    ps = [p => 2.0, d => 1.0, ω => 0.5, k => 2.0]

    # Checks that ODE an simulation of the system achieves the correct steady state.
    oprob = ODEProblem(coupled_rs, u0, (0.0, 1000.0), ps; structural_simplify = true)
    osol = solve(oprob, Vern7(); abstol = 1e-8, reltol = 1e-8)
    @test osol[X][end] ≈ 2.0
    @test osol[A][end] ≈ 0.0 atol = 1e-8
    @test osol[D(A)][end] ≈ 0.0 atol = 1e-8
    @test osol[B][end] ≈ 1.0

    # Checks that SteadyState simulation of the system achieves the correct steady state.
    ssprob = SteadyStateProblem(coupled_rs, u0, ps; structural_simplify = true)
    sssol = solve(ssprob, DynamicSS(Vern7()); abstol = 1e-8, reltol = 1e-8)
    @test sssol[X][end] ≈ 2.0
    @test sssol[A][end] ≈ 0.0 atol = 1e-8
    @test sssol[D(A)][end] ≈ 0.0 atol = 1e-8
    @test sssol[B][end] ≈ 1.0

    # Checks that the steady state can be found by solving a nonlinear problem.
    # Here `B => 0.1` has to be provided as well (and it shouldn't for the 2nd order ODE), hence the
    # separate `u0` declaration.
    u0 = [X => 1.0, A => 2.0, D(A) => 1.0, B => 0.1]
    nlprob = NonlinearProblem(coupled_rs, u0, ps; structural_simplify = true, all_differentials_permitted = true)
    nlsol = solve(nlprob)
    @test nlsol[X][end] ≈ 2.0
    @test nlsol[A][end] ≈ 0.0
    @test nlsol[B][end] ≈ 1.0
end


### DSL Tests ###

# Check that a coupled CRN/DAE created programmatically and via the DSL are identical.
# Checks where variables are implied from differential equations, and with variables/parameter
# default values, types, and metadata.
# Checks that generated system contents are correct, and ODE simulations are identical.
let
    # Creates the model programmatically.
    @species X1(t) X2(t) X3(t)
    @variables V(t)=5.0 [description="Volume"] N(t) X_conc(t) X_tot(t)
    @parameters p k1 k2 d v n x_scale::Float32
    eqs = [
        Reaction(p, nothing, [X1])
        Reaction(k1, [X1], [X2])
        Reaction(k2, [X2], [X3])
        Reaction(d, [X3], nothing)
        D(V) ~ X3/(1+X3) - v*V
        D(N) ~ - n*N*X3
        V*X_conc ~ x_scale*(X1 + X2 + X3)
        X_tot + X1 + X2 ~ -X3
    ]
    rs_prog = complete(ReactionSystem(eqs, t; name = :coupled_rs))

    # Creates the model via the DSL.
    rs_dsl = @reaction_network coupled_rs begin
        @variables X_conc(t) V(t)=5.0 [description="Volume"] X_tot(t)
        @parameters v n x_scale::Float32
        @equations begin
            D(V) ~ X3/(1+X3) - v*V
            D(N) ~ - n*N*X3
            V*X_conc ~ x_scale*(X1 + X2 + X3)
            X_tot + X1 + X2 ~ -X3
        end
        p, 0 --> X1
        k1, X1 --> X2
        k2, X2 --> X3
        d, X3 --> 0
    end

    # Checks that models are identical. Also checks that they have the correct content.
    @test Catalyst.isequivalent(rs_prog, rs_dsl)
    @test getdescription(rs_dsl.V) == "Volume"
    @test getdefault(rs_dsl.V) == 5.0
    @test SymbolicUtils.symtype(rs_dsl.x_scale) == Float32

    @test issetequal(parameters(rs_dsl), [p, k1, k2, d, v, n, x_scale])
    @test issetequal(species(rs_dsl), unknowns(rs_dsl)[1:3])
    @test issetequal(unknowns(rs_dsl)[1:3], [X1, X2, X3])
    @test issetequal(unknowns(rs_dsl)[4:7], [V, N, X_conc, X_tot])
    @test issetequal(reactions(rs_dsl), equations(rs_dsl)[1:4])
    @test issetequal(equations(rs_dsl)[1:4], eqs[1:4])
    @test issetequal(equations(rs_dsl)[5:7], eqs[5:7])


    # Checks that the models can be simulated and yield identical results.
    # Test most likely redundant, but seem useful to have one test like this to be sure.
    u0 = [X1 => 0.1, X2 => 0.2, X3 => 0.2, X_tot => 0.6, N => 10.0, X_conc => 10.0]
    ps = [p => 1.0, k1 => 1.2, k2 => 1.5, d => 2.0, v => 0.2, n => 0.5, x_scale => 2.0]
    oprob_prog = ODEProblem(rs_prog, u0, (0.0, 10.0), ps; structural_simplify = true, warn_initialize_determined = false)
    oprob_dsl = ODEProblem(rs_dsl, u0, (0.0, 10.0), ps; structural_simplify = true, warn_initialize_determined = false)
    @test solve(oprob_prog, Rosenbrock23()) == solve(oprob_dsl, Rosenbrock23())
end

# Checks that equations can both be declared in a single line, or within a `begin ... end` block.
let
    # Checks for system with a single differential equation.
    rs_1_line = @reaction_network rs_1 begin
        @equations D(M) ~ -M*I
        i, S + I --> 2I
        r, I --> R
    end
    rs_1_block = @reaction_network rs_1 begin
        @equations begin
            D(M) ~ -M*I
        end
        i, S + I --> 2I
        r, I --> R
    end
    @test Catalyst.isequivalent(rs_1_line, rs_1_block)

    # Checks for system with a single algebraic equation.
    rs_2_line = @reaction_network rs_2 begin
        @variables H(t)
        @equations H ~ 100 - I
        i, S + I --> 2I
        r, I --> R
    end
    rs_2_block = @reaction_network rs_2 begin
        @variables H(t)
        @equations begin
            H ~ 100 - I
        end
        i, S + I --> 2I
        r, I --> R
    end
    @test Catalyst.isequivalent(rs_2_line, rs_2_block)
end

# Checks that lhs variable is correctly inferred from differential equations.
let
    # Checks for system with a differential equation and an algebraic equation.
    # Here, `H` is defined using `@variables`, but M should be inferred.
    rs_1 = @reaction_network begin
        @variables H(t)
        @equations begin
            D(M) ~ -M*I
            H ~ 100 - I
        end
        i, S + I --> 2I
        r, I --> R
    end
    issetequal(species(rs_1), [rs_1.S, rs_1.I, rs_1.R])
    issetequal(unknowns(rs_1)[4:5], [rs_1.H, rs_1.M])

    # Checks for system with two differential equations, and which do not use `@variables`,
    rs_2 = @reaction_network coupled_rs begin
        @equations begin
            D(V) ~ X/(1+X) - V
            D(N) ~ - V
        end
        (p,d), 0 <--> X
    end
    issetequal(species(rs_2), [rs_2.X])
    issetequal(unknowns(rs_2)[2:3], [rs_2.V, rs_2.N])

    # Checks for system with two differential equations, where one is defined using `@variables`.
    rs_2 = @reaction_network coupled_rs begin
        @variables N(t)
        @equations begin
            D(V) ~ X/(1+X) - V
            D(N) ~ - V
        end
        (p,d), 0 <--> X
    end
    issetequal(species(rs_2), [rs_2.X])
    issetequal(unknowns(rs_2)[2:3], [rs_2.V, rs_2.N])
end

# Checks that variables that can be inferred from differential equations, but are also declared
# manually, have their additional inputs properly registered.
let
    rs = @reaction_network begin
        @variables V(t)=2.0 [description = "A variable"]
        @equations D(V) ~ -1
    end
    @test getdefault(rs.V) == 2.0
    @test getdescription(rs.V) == "A variable"
end

# Checks that equations can be formatted in various ways. Tries e.g. isolating a single number on
# either side of the equality.
# Checks that various weird function can be used within equations.
# Checks that special symbols, like π and t can be used within equations.
let
    # Declares models with a single equation, formatted in various ways.
    rs_1 = @reaction_network rs begin
        @parameters p q
        @species X(t)
        @variables A(t) B(t)
        @equations X^2 + log(A+X) ~ 1 - sqrt(B) + sin(p + X + π)/exp(A/(1+t)) + q
    end
    rs_2 = @reaction_network rs begin
        @parameters p q
        @species X(t)
        @variables A(t) B(t)
        @equations X^2 + log(A+X) + sqrt(B) - sin(p + X + π)/exp(A/(1+t)) - q ~ 1
    end
    rs_3 = @reaction_network rs begin
        @parameters p q
        @species X(t)
        @variables A(t) B(t)
        @equations X^2 + log(A+X) + sqrt(B) - sin(p + X + π)/exp(A/(1+t)) - 1 - q ~ 0
    end
    rs_4 = @reaction_network rs begin
        @parameters p q
        @species X(t)
        @variables A(t) B(t)
        @equations 0 ~ X^2 + log(A+X) + sqrt(B) - sin(p + X + π)/exp(A/(1+t)) - 1 - q
    end
    rs_5 = @reaction_network rs begin
        @parameters p q
        @species X(t)
        @variables A(t) B(t)
        @equations q ~ X^2 + log(A+X) + sqrt(B) - sin(p + X + π)/exp(A/(1+t)) - 1
    end
    rs_6 = @reaction_network rs begin
        @parameters p q
        @species X(t)
        @variables A(t) B(t)
        @equations X^2 + log(A+X) + (A + B)^p ~ 1 - sqrt(B) + sin(p + X + π)/exp(A/(1+t)) + q + (A + B)^p
    end

    # Uses a special function to check that all equations indeed are identical.
    function is_eqs_equal(rs1, rs2; eq_idx = 1)
        eq1 = equations(rs1)[eq_idx]
        eq2 = equations(rs2)[eq_idx]
        ModelingToolkitBase._iszero(eq1.lhs - eq1.rhs - eq2.lhs + eq2.rhs) && return true
        ModelingToolkitBase._iszero(eq1.lhs - eq1.rhs + eq2.lhs - eq2.rhs) && return true
        return false
    end
    @test_broken is_eqs_equal(rs_1, rs_2) # Test broken due to https://github.com/JuliaSymbolics/Symbolics.jl/issues/1739, not a Catalyst problem.
    @test_broken is_eqs_equal(rs_1, rs_3) # Test broken due to https://github.com/JuliaSymbolics/Symbolics.jl/issues/1739, not a Catalyst problem.
    @test_broken is_eqs_equal(rs_1, rs_4) # Test broken due to https://github.com/JuliaSymbolics/Symbolics.jl/issues/1739, not a Catalyst problem.
    @test_broken is_eqs_equal(rs_1, rs_5) # Test broken due to https://github.com/JuliaSymbolics/Symbolics.jl/issues/1739, not a Catalyst problem.
    @test is_eqs_equal(rs_1, rs_6)
end

# Checks that the default differential (`D`) uses a declared, non-default, independent variable.
# Check that inferred variables depends on declared time independent variables.
let
    # Declares model.
    rs = @reaction_network begin
        @ivs τ
        @equations D(V) ~ -1.0
    end

    # Checks that the default differential uses τ iv.
    Ds = Differential(ModelingToolkitBase.get_iv(rs))
    @test isequal(operation(equations(rs)[1].lhs), Ds)

    # Checks that the inferred variable depends on τ iv.
    @variables V($(ModelingToolkitBase.get_iv(rs)))
    @test isequal(V, rs.V)
end

# Checks that custom differentials can be declared.
# Checks that non-default ivs work, and that a new differential using this can overwrite the default one.
let
    # Declares the reaction system using the default differential and iv.
    rs_1 = @reaction_network begin
        @equations D(N) ~ -N
        (p,d), 0 <--> X
    end

    # Declares the reaction system using a new iv, and overwriting the default differential.
    rs_2 = @reaction_network begin
        @ivs τ
        @species X(τ)
        @variables N(τ)
        @differentials D = Differential(τ)
        @equations D(N) ~ -N
        (p,d), 0 <--> X
    end

    # Declares the reaction system using a new differential and iv.
    rs_3 = @reaction_network begin
        @ivs τ
        @species X(τ)
        @variables N(τ)
        @differentials Δ = Differential(τ)
        @equations Δ(N) ~ -N
        (p,d), 0 <--> X
    end

    # Simulates all three models, checking that the results are identical.
    u0 = [:X => 5.0, :N => 10.0]
    tspan = (0.0, 10.)
    ps = [:p => 1.0, :d => 0.2]
    oprob_1 = ODEProblem(rs_1, u0, tspan, ps)
    oprob_2 = ODEProblem(rs_2, u0, tspan, ps)
    oprob_3 = ODEProblem(rs_3, u0, tspan, ps)
    @test solve(oprob_1, Tsit5()) == solve(oprob_2, Tsit5()) == solve(oprob_3, Tsit5())
end

# Checks that various misformatted declarations yield errors.
let
    # Attempting to create a new differential from an unknown iv.
    @test_throws Exception @eval @reaction_network begin
        @differentials D = Differential(τ)
    end

    # Misformatted expression for a differential.
    @test_throws Exception @eval @reaction_network begin
        @variables D
        @differentials d ~ D
    end

    # Several equations without `begin ... end` block.
    @test_throws Exception @eval @reaction_network begin
        @equations D(V) + 1 ~ - 1.0
        @equations D(W) + 1 ~ - 1.0
    end

    # System using multiple ivs.
    @test_throws Exception @eval @reaction_network begin
        @ivs τ Τ
        @variables n(τ) N(Τ)
        @differentials begin
            δ = Differential(τ)
            Δ = Differential(Τ)
        end
        @equations begin
            δ(n) ~ -n
            Δ(N) ~ -N
        end
    end
end


### Error Tests ###

# Checks that various erroneous coupled system declarations yield errors.
let
    @parameters p1 p2
    @parameters τ
    @variables U1(τ) V1(t)
    @species R1(τ) R2(τ) S1(t) S2(t)
    E = Differential(τ)

    # Variables as reaction reactants.
    @test_throws Exception ReactionSystem([
        Reaction(p1, [S1], [V1])
    ], t; name = :rs)

    # Equation with variable using non-declared independent variable.
    @test_throws Exception ReactionSystem([
        Reaction(p1, [S1], [S2]),
        E(U1) ~ S1 + p2
    ], t; name = :rs)

    # Differential with respect to non-declared independent variable.
    @test_throws Exception ReactionSystem([
        Reaction(p1, [S1], [S2]),
        E(V1) ~ S1 + p2
    ], [t, τ]; name = :rs)
end

# Checks that various attempts to create `ODEProblem`s from faulty systems generate errors.
let
    @parameters p1 p2
    @variables V1(t)
    @species S1(t) S2(t)

    # Coupled system overconstrained due to additional algebraic equations (with variables).
    eqs = [
        Reaction(p1, [S1], [S2]),
        V1 ~ p2 - S1,
        S2 ~ V1^2 + sqrt(S2)
    ]
    @named rs = ReactionSystem(eqs, t)
    rs = complete(rs)
    u0 = [S1 => 1.0, S2 => 2.0, V1 => 0.1]
    ps = [p1 => 2.0, p2 => 3.0]
    @test_throws Exception ODEProblem(rs, u0, (0.0, 1.0), ps; structural_simplify = true)
end

# Checks that equations cannot contain differentials with respect to species.
let
    # Basic case.
    @test_throws Exception @eval @reaction_network begin
        @equations D(X) ~ 1.0
        d, X --> 0
    end

    # Complicated differential.
    @test_throws Exception @eval @reaction_network begin
        @equations V + X^2 ~ 1.0 - D(V + log(1 + X))
        d, X --> 0
    end

    # Case where the species is declared, but not part of a reaction.
    @test_throws Exception @eval @reaction_network begin
        @species Y(t)
        @equations D(Y) ~ 1.0
        d, X --> 0
    end

    # Case where the equation also declares a new non-species variable.
    # At some point something like this could be supported, however, not right now.
    @test_throws Exception @eval @reaction_network begin
        @equations D(V) ~ 1.0 + D(X)
        d, X --> 0
    end
end
