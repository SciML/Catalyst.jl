# Fetch packages.
using Catalyst, NonlinearSolve, OrdinaryDiffEq, Statistics, SteadyStateDiffEq, StochasticDiffEq, Test
using ModelingToolkit: getdefault
using Symbolics: BasicSymbolic, unwrap

# Sets stable rng number.
using StableRNGs
rng = StableRNG(12345)
seed = rand(rng, 1:100)

# Sets the default `t` to use.
t = default_t()
D = default_time_deriv()


### Basic Differential Hybrid Tests ###

# Tests hybrid CRN/ODE. Checks that known steady state is reached.
# Check that steady state can be found using NonlinearSolve and SteadyStateDiffEq.
let
    # Creates hybrid reactions system.
    @parameters p d k
    @species X(t)
    @variables A(t)
    eqs = [
        Reaction(p, nothing, [X]),
        Reaction(d, [X], nothing),
        D(A) ~ p*X - k*A
    ]
    @named hybrid_rs = ReactionSystem(eqs, t)
    hybrid_rs = complete(hybrid_rs)

    # Basic model checks.
    @test issetequal(parameters(hybrid_rs), [p, d, k])
    @test issetequal(species(hybrid_rs), unknowns(hybrid_rs)[1:1])
    @test issetequal(unknowns(hybrid_rs)[1:1], [X])
    @test issetequal(unknowns(hybrid_rs)[2:2], [A])
    @test issetequal(reactions(hybrid_rs), equations(hybrid_rs)[1:2])
    @test issetequal(equations(hybrid_rs)[1:2], eqs[1:2])
    @test issetequal(equations(hybrid_rs)[3:3], eqs[3:3])

    # Set simulation inputs.
    u0 = [X => 0.1, A => 10.0]
    tspan = (0.0, 1000.0)
    ps = [p => 1.0, d => 0.5, k => 2.0]

    # Checks that the correct steady state is found through ODEProblem.
    oprob = ODEProblem(hybrid_rs, u0, tspan, ps)
    osol = solve(oprob, Vern7(); abstol = 1e-8, reltol = 1e-8)
    @test osol[end] ≈ [2.0, 1.0]

    # Checks that the correct steady state is found through NonlinearProblem.
    nlprob = NonlinearProblem(hybrid_rs, u0, ps)
    nlsol = solve(nlprob; abstol = 1e-8, reltol = 1e-8)
    @test nlsol ≈ [2.0, 1.0]

    # Checks that the correct steady state is found through SteadyStateProblem.
    ssprob = SteadyStateProblem(hybrid_rs, u0, ps)
    sssol = solve(ssprob, DynamicSS(Rosenbrock23()); abstol = 1e-8, reltol = 1e-8)
    @test sssol ≈ [2.0, 1.0]
end

# Checks that hybrid systems created via the DSL, extension, and programmatically are identical.
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
    hybrid_rs_prog = ReactionSystem(eqs_prog, t, [A, B, X1, X2], [k1, k2, a, b]; name = :hybrid_rs)
    hybrid_rs_prog = complete(hybrid_rs_prog)

    # Creates model by extending a `ReactionSystem` with a ODESystem.
    rn_extended = @network_component begin
        ($k1*$A, $k2*$B), X1 <--> X2 
    end
    eqs_extended = [
        D(A) ~ X1 + a - A
        D(B) ~ X2 + b - B
    ]
    @named osys_extended = ODESystem(eqs_extended, t)
    hybrid_rs_extended = complete(extend(osys_extended, rn_extended; name = :hybrid_rs))

    # Creates the model through the DSL.
    hybrid_rs_dsl = @reaction_network hybrid_rs begin
        @parameters k1 k2 a b
        @equations begin
            D(A) ~ X1 + a - A
            D(B) ~ X2 + b - B
        end
        (k1*A, k2*B), X1 <--> X2 
    end

    # Checks that models are equivalent and contain the correct stuff.
    @test hybrid_rs_prog == hybrid_rs_extended == hybrid_rs_dsl
    @test issetequal(parameters(hybrid_rs_extended), [a, b, k1, k2])
    @test issetequal(species(hybrid_rs_extended), [X1, X2])
    @test issetequal(unknowns(hybrid_rs_extended)[1:2], [X1, X2])
    @test issetequal(unknowns(hybrid_rs_extended)[3:4], [A, B])
    @test issetequal(equations(hybrid_rs_extended)[3:4], eqs_extended)

    # Simulates the three models, checking that they all yield the correct end point.
    u0 = [A => 1.0, B => 1.0, X1 => 10.0, X2 => 10.0]
    tspan = (0.0, 100.)
    ps = [a => 1.0, b => 1.0, k1 => 1.0, k2 => 1.0]
    for hybrid_rs in [hybrid_rs_prog, hybrid_rs_extended, hybrid_rs_dsl]
        oprob = ODEProblem(hybrid_rs, u0, tspan, ps)
        osol = solve(oprob, Vern7(); abstol = 1e-8, reltol = 1e-8)
        osol[end] ≈ [10.0, 10.0, 11.0, 11.0]
    end
end


### Basic Algebraic Hybrid Tests ###

# Tests hybrid CRN/algebraic equation. Checks that known steady state is reached using ODE solve.
# Check that steady state can be found using NonlinearSolve and SteadyStateDiffEq.
# Checks that errors are given if `structural_simplify = true` argument is not given.
let 
    # Creates a simple hybrid model with an algebraic equation.
    @parameters p d a b
    @species X(t)
    @variables A(t)
    eqs = [
        Reaction(p, nothing, [X]),
        Reaction(d, [X], nothing),
        a*A^2 ~ X + b
    ]
    @named hybrid_rs = ReactionSystem(eqs, t)
    hybrid_rs = complete(hybrid_rs)

    # Check model content.
    @test issetequal(parameters(hybrid_rs), [p, d, a, b])
    @test issetequal(species(hybrid_rs), unknowns(hybrid_rs)[1:1])
    @test issetequal(unknowns(hybrid_rs)[1:1], [X])
    @test issetequal(unknowns(hybrid_rs)[2:2], [A])
    @test issetequal(reactions(hybrid_rs), equations(hybrid_rs)[1:2])
    @test issetequal(equations(hybrid_rs)[1:2], eqs[1:2])
    @test issetequal(equations(hybrid_rs)[3:3], eqs[3:3])

    # Set simulation inputs.
    u0 = [X => 0.1, A => 10.0]
    tspan = (0.0, 1000.0)
    ps = [p => 1.0, d => 0.5, a => 2.0, b => 16.0]

    # Checks not using `structural_simplify` argument yields an error.
    @test_throws Exception ODEProblem(hybrid_rs, u0, tspan, ps)
    @test_throws Exception SteadyStateProblem(hybrid_rs, u0, ps)

    # Checks that the correct steady state is found through ODEProblem.
    oprob = ODEProblem(hybrid_rs, u0, tspan, ps; structural_simplify = true)
    osol = solve(oprob, Rosenbrock23(); abstol = 1e-8, reltol = 1e-8)
    @test osol[end] ≈ [2.0, 3.0]

    # Checks that the correct steady state is found through NonlinearProblem.
    nlprob = NonlinearProblem(hybrid_rs, u0, ps)
    nlsol = solve(nlprob)
    @test nlsol ≈ [2.0, 3.0]

    # Checks that the correct steady state is found through SteadyStateProblem.
    ssprob = SteadyStateProblem(hybrid_rs, u0, ps; structural_simplify = true)
    sssol = solve(ssprob, DynamicSS(Rosenbrock23()); abstol = 1e-8, reltol = 1e-8)
    @test sssol ≈ [2.0, 3.0]
end


### Basic Combined Algebraic/Hybrid Hybrid Tests ###

# Checks that a combined reaction/differential/algebraic hybrid system can be created.
# Checks that it can its ODE, SteadyState, and Nonlinear problems all can be solved.
# Checks that Tuple u0/ps input, and non-default independent variables works.
# The system is mostly made up to be non-trivial, but reliably solvable.
let
    @parameters p d a b c
    @variables τ A(τ) B(τ) C(τ)
    @species X(τ)
    Δ = Differential(τ)
    eqs = [
        Δ(A) ~ b + X - A,
        Δ(B) ~ sqrt(A + X + b) - B,
        Reaction(p, nothing, [X], nothing, [2]),
        Reaction(d, [X], nothing),
        (X + C)*B ~ A
    ]
    @named hybrid_rs = ReactionSystem(eqs, τ)
    hybrid_rs = complete(hybrid_rs)

    # Set simulation inputs.
    u0 = (X => 1.0, A => 2.0, B => 3.0, C => 4.0)
    ps = (p => 1.0, d => 2.0, a => 3.0, b => 4.0, c => 5.0)

    # Creates and solves a ODE, SteadyState, and Nonlinear problems.
    # Success is tested by checking that the same steady state solution is found.
    oprob = ODEProblem(hybrid_rs, u0, (0.0, 1000.0), ps; structural_simplify = true)
    ssprob = SteadyStateProblem(hybrid_rs, u0, ps; structural_simplify = true)
    nlprob = NonlinearProblem(hybrid_rs, u0, ps)
    osol = solve(oprob, Rosenbrock23(); abstol = 1e-8, reltol = 1e-8)
    sssol = solve(ssprob, DynamicSS(Rosenbrock23()); abstol = 1e-8, reltol = 1e-8)
    nlsol = solve(nlprob; abstol = 1e-8, reltol = 1e-8)
    @test osol[end] ≈ sssol ≈ nlsol
end


### Species, Variables, and Parameter Handling ###

# Checks that hybrid systems contain the correct species, variables, and parameters.
# Checks that species, variables, and parameters are inferred correctly from equations.
# Checks that non-default iv is inferred correctly from reactions/equations.
let 
    # Create hybrid model.
    @variables τ A(τ) B(τ)
    @species X(τ) X2(τ)
    @parameters k1 k2 k b1 b2
    D = Differential(τ)
    eqs = [
        Reaction(k1, [X], [X2], [2], [1]),
        Reaction(k2, [X2], [X], [1], [2]),
        D(A) ~ k*X2 - A,
        B + A ~ b1*X + b2*X2
    ]
    @named hybrid_rs = ReactionSystem(eqs, τ)
    hybrid_rs = complete(hybrid_rs)

    # Checks that systems created from hybrid reaction systems contain the correct content (in the correct order).
    osys = convert(ODESystem, hybrid_rs)
    ssys = convert(SDESystem, hybrid_rs)
    nlsys = convert(NonlinearSystem, hybrid_rs)
    for sys in [hybrid_rs, osys, ssys, nlsys]
        @test issetequal(parameters(sys), [k1, k2, k, b1, b2])
        @test issetequal(unknowns(sys)[1:2], [X, X2])
        @test issetequal(unknowns(sys)[3:4], [A, B])
    end
end

# Checks that parameters, species, and variables can be correctly accessed in hybrid systems.
# Checks for both differential and algebraic equations.
# Checks for problems, integrators, and solutions yielded by hybrid systems.
# Checks that metadata, types, and default values are carried through correctly.
let
    # Creates the model
    @parameters a1 [description="Parameter a1"] a2::Rational{Int64} a3=0.3 a4::Rational{Int64}=4//10 [description="Parameter a4"]
    @parameters b1 [description="Parameter b1"] b2::Int64 b3 = 3 b4::Int64=4 [description="Parameter b4"]
    @parameters c1 [description="Parameter c1"] c2::Float32 c3=30.0 c4::Float32=40.0 [description="Parameter c4"]
    @species A1(t) [description="Species A1"] A2(t)=0.2 A3(t)=0.3 [description="Species A3"] A4(t)
    @variables B1(t) [description="Variable B1"] B2(t)=2.0 B3(t)=3.0 [description="Variable B3"] B4(t)
    @variables C1(t) [description="Variable C1"] C2(t)=20.0 C3(t)=30.0 [description="Variable C3"] C4(t)
    eqs = [
        Reaction(a1, nothing, [A1]),
        Reaction(a2, nothing, [A2]),
        Reaction(a3, nothing, [A3]),
        Reaction(a4, nothing, [A4]),
        D(B1) ~ b1*B1,
        D(B2) ~ b2*B2,
        D(B3) ~ b3*B3,
        D(B4) ~ b4*B4,
        C1^2 ~ c1 + B1^5,
        C2^2 ~ c2 + B2^5,
        C3^2 ~ c3 + B3^5,
        C4^2 ~ c4 + B4^5
    ]
    @named hybrid_rs = ReactionSystem(eqs, t)
    hybrid_rs = complete(hybrid_rs)

    # Checks that the model has the correct content.
    @test issetequal(parameters(hybrid_rs), [a1, a2, a3, a4, b1, b2, b3, b4, c1, c2, c3, c4])
    @test issetequal(species(hybrid_rs), unknowns(hybrid_rs)[1:4])
    @test issetequal(unknowns(hybrid_rs)[1:4], [A1, A2, A3, A4])
    @test issetequal(unknowns(hybrid_rs)[5:12], [B1, B2, B3, B4, C1, C2, C3, C4])
    @test issetequal(reactions(hybrid_rs)[1:4], equations(hybrid_rs)[1:4])
    @test issetequal(equations(hybrid_rs)[1:4], eqs[1:4])
    @test issetequal(equations(hybrid_rs)[5:12], eqs[5:12])

    # Checks that parameters, species, and variables carried the correct information.
    @test unwrap(hybrid_rs.a1) isa BasicSymbolic{Real}
    @test unwrap(hybrid_rs.a2) isa BasicSymbolic{Rational{Int64}}
    @test unwrap(hybrid_rs.a3) isa BasicSymbolic{Real}
    @test unwrap(hybrid_rs.a4) isa BasicSymbolic{Rational{Int64}}
    @test unwrap(hybrid_rs.b1) isa BasicSymbolic{Real}
    @test unwrap(hybrid_rs.b2) isa BasicSymbolic{Int64}
    @test unwrap(hybrid_rs.b3) isa BasicSymbolic{Real}
    @test unwrap(hybrid_rs.b4) isa BasicSymbolic{Int64}
    @test unwrap(hybrid_rs.c1) isa BasicSymbolic{Real}
    @test unwrap(hybrid_rs.c2) isa BasicSymbolic{Float32}
    @test unwrap(hybrid_rs.c3) isa BasicSymbolic{Real}
    @test unwrap(hybrid_rs.c4) isa BasicSymbolic{Float32}
    @test getdescription(hybrid_rs.a1) == "Parameter a1"
    @test getdescription(hybrid_rs.a4) == "Parameter a4"
    @test getdescription(hybrid_rs.b1) == "Parameter b1"
    @test getdescription(hybrid_rs.b4) == "Parameter b4"
    @test getdescription(hybrid_rs.c1) == "Parameter c1"
    @test getdescription(hybrid_rs.c4) == "Parameter c4"
    @test getdefault(hybrid_rs.a3) == 0.3
    @test getdefault(hybrid_rs.a4) == 4//10
    @test getdefault(hybrid_rs.b3) == 3
    @test getdefault(hybrid_rs.b4) == 4
    @test getdefault(hybrid_rs.c3) == 30
    @test getdefault(hybrid_rs.c4) == 40
    @test getdescription(hybrid_rs.A1) == "Species A1"
    @test getdescription(hybrid_rs.A3) == "Species A3"
    @test getdescription(hybrid_rs.B1) == "Variable B1"
    @test getdescription(hybrid_rs.B3) == "Variable B3"
    @test getdescription(hybrid_rs.C1) == "Variable C1"
    @test getdescription(hybrid_rs.C3) == "Variable C3"
    @test getdefault(hybrid_rs.A2) == 0.2
    @test getdefault(hybrid_rs.A3) == 0.3
    @test getdefault(hybrid_rs.B2) == 2.0
    @test getdefault(hybrid_rs.B3) == 3.0
    @test getdefault(hybrid_rs.C2) == 20.0
    @test getdefault(hybrid_rs.C3) == 30.0

    # Creates problem inputs.
    u0 = [a1 => 0.1, a2 => 2//10, b1 => 1.0, b2 => 2, c1 => 10.0, c2 => 20.0]
    tspan = (0.0, 1.0)
    ps = [A1 => 0.1, B1 => 1.0, C1 => 10.0]

    # Create ODE structures.
    oprob = ODEProblem(hybrid_rs, u0, tspan, ps; structural_simplify = true)
    oint = init(oprob, Tsit5())
    osol = solve(oprob, Tsit5())

    # Create SDE structures.
    sprob = SDEProblem(hybrid_rs, u0, tspan, ps)
    sint = init(oprob, ImplicitEM())
    ssol = solve(oprob, ImplicitEM())

    # Creates Nonlinear structures.
    nlprob = NonlinearProblem(hybrid_rs, u0, ps)
    nlint = init(nlprob, NewtonRaphson())
    nlsol = solve(nlprob, NewtonRaphson())

    # Checks indexing.
    for mtk_struct in [oprob, oint, osol, sprob, sint, ssol, nlprob, nlint, nlsol]
        # Parameters.
        @test mtk_struct[a1] == 0.1
        @test mtk_struct[a2] == 2//10
        @test mtk_struct[a3] == 0.3
        @test mtk_struct[a4] == 4//10
        @test mtk_struct[b1] == 1.0
        @test mtk_struct[b2] == 2
        @test mtk_struct[b3] == 3.0
        @test mtk_struct[b4] == 4
        @test mtk_struct[c1] == 10.0
        @test mtk_struct[c2] == 20.0
        @test mtk_struct[c3] == 30.0
        @test mtk_struct[c4] == 40.0

        # Species.
        @test mtk_struct[A1] == 0.1
        @test mtk_struct[A2] == 2//10
        @test mtk_struct[A3] == 0.3
        @test mtk_struct[A4] == 4//10

        # Variables.
        @test mtk_struct[B1] == 1.0
        @test mtk_struct[B2] == 2
        @test mtk_struct[B3] == 3.0
        @test mtk_struct[B4] == 4
        @test mtk_struct[C1] == 10.0
        @test mtk_struct[C2] == 20.0
        @test mtk_struct[C3] == 30.0
        @test mtk_struct[C4] == 40.0
    end
end


### Hybrid SDE Tests ###

# Checks that a hybrid SDE + differential equations works.
# Checks that CLE noise does not affect ODE part that should be deterministic.
# Only considers added differential equations without noise.
let
    # Creates hybrid reactions system.
    @parameters p d k1 k2
    @species X(t)
    @variables A(t) B(t)
    eqs = [
        Reaction(p, nothing, [X]; metadata = [:noise_scaling => 0.1]),
        Reaction(d, [X], nothing; metadata = [:noise_scaling => 0.1]),
        D(A) ~ X - k1*A,
        D(B) ~ k2 - B
    ]
    @named hybrid_rs = ReactionSystem(eqs, t)
    hybrid_rs = complete(hybrid_rs)

    # Set simulation inputs.
    u0 = [X => 100.0, A => 50.0, B => 2.0]
    tspan = (0.0, 10000.0)
    ps = [p => 10.0, d => 0.1, k1 => 2.0, k2 => 20.0]

    # Checks that the simulations have the expected means (or endpoint, for B).
    sprob = SDEProblem(hybrid_rs, u0, tspan, ps)
    ssol = solve(sprob, ImplicitEM(); maxiters = 1e9, seed)
    @test mean(ssol[:X]) ≈ 100.0 atol = 1e-2 rtol = 1e-2
    @test mean(ssol[:A]) ≈ 50.0 atol = 1e-2 rtol = 1e-2
    @test ssol[:B][end] ≈ 20.0
end

# Checks that a hybrid SDE + algebraic equations works.
# Checks that structural_simplify is required to simulate hybrid SDE + algebraic equations.
let 
    # Creates hybrid reactions system.
    @parameters p d k1 k2
    @species X(t)
    @variables A(t)
    eqs = [
        Reaction(p, nothing, [X]),
        Reaction(d, [X], nothing),
        2 + k1 * A ~ 3 + k2 * X
    ]
    @named hybrid_rs = ReactionSystem(eqs, t)
    hybrid_rs = complete(hybrid_rs)

    # Set simulation inputs.
    u0 = [X => 100.0, A => 10.0]
    tspan = (0.0, 1000.0)
    ps = Dict([p => 1.0, d => 0.01, k1 => 3.0, k2 => 4.0])

    # Check that the structural_simplify argument is required.
    @test_throws Exception SDEProblem(hybrid_rs, u0, tspan, ps)

    # Checks the algebraic equation holds.
    sprob = SDEProblem(hybrid_rs, u0, tspan, ps; structural_simplify = true)
    ssol = solve(sprob, ImplicitEM())
    @test 2 .+ ps[:k1] * ssol[:A] == 3 .+ ps[:k2] * ssol[:X]    
end

### Unusual Differentials Tests ###

# Tests that hybrid CRN/DAEs with higher order differentials can be created.
# Tests that these can be solved using ODEs, nonlinear solving, and steady state simulations.
let 
    # Create hybrid model.
    @species X(t)
    @variables A(t) B(t)
    @parameters p d ω k
    eqs = [
        Reaction(p, nothing, [X]),
        Reaction(d, [X], nothing),
        D(D(A)) + 2ω*D(A) +(ω^2)*A ~ 0,
        A + k*(B + D(A)) ~ X 
    ]
    @named hybrid_rs = ReactionSystem(eqs, t)
    hybrid_rs = complete(hybrid_rs)
    u0 = [X => 1.0, A => 2.0, D(A) => 1.0]
    ps = [p => 2.0, d => 1.0, ω => 0.5, k => 2.0]

    # Checks that ODE an simulation of the system achieves the correct steady state.
    oprob = ODEProblem(hybrid_rs, u0, (0.0, 1000.0), ps; structural_simplify = true)
    osol = solve(oprob, Vern7(); abstol = 1e-8, reltol = 1e-8)
    @test osol[X][end] ≈ 2.0
    @test osol[A][end] ≈ 0.0 atol = 1e-8
    @test osol[D(A)][end] ≈ 0.0 atol = 1e-8
    @test osol[B][end] ≈ 1.0

    # Checks that SteadyState simulation of the system achieves the correct steady state.
    # Currently broken due to MTK.
    @test_broken begin
        ssprob = SteadyStateProblem(hybrid_rs, u0, ps; structural_simplify = true)
        sssol = solve(oprob, DynamicSS(Vern7()); abstol = 1e-8, reltol = 1e-8)
        @test osol[X][end] ≈ 2.0
        @test osol[A][end] ≈ 0.0 atol = 1e-8
        @test osol[D(A)][end] ≈ 0.0 atol = 1e-8
        @test osol[B][end] ≈ 1.0
    end

    # Checks that the steady state can be found by solving a nonlinear problem.
    # Here `B => 0.1` has to be provided as well (and it shouldn't for the 2nd order ODE), hence the 
    # separate `u0` declaration.
    u0 = [X => 1.0, A => 2.0, D(A) => 1.0, B => 0.1]
    nlprob = NonlinearProblem(hybrid_rs, u0, ps; structural_simplify = true)
    nlsol = solve(nlprob)
    @test nlsol[X][end] ≈ 2.0
    @test nlsol[A][end] ≈ 0.0
    @test nlsol[B][end] ≈ 1.0
end

# Checks that DAEs are created properly when provided disorderly.
# Checks that differential equations can provided in a form no `D(...) ~ ...` (i.e. several
# differentials, not necessarily on the same side).
# Checks with non-default iv, and parameters/initial conditions given using Symbols.
# Checks with default value for algebraic variable.
let 
    # Prepares stuff common to both simulations.
    @parameters i r m1 m2 h_max
    u0 = [:S => 999.0, :I => 1.0, :R => 0.0, :M => 1000.0]
    tspan = 500.0
    ps = [:i => 0.0001, :r => 0.01, :m1 => 5000.0, :m2 => 9000//3, :h_max => 1500.0]

    # Declares the model in an ordered fashion, and simulates it.
    osol_ordered = let
        @variables M(t) H(t)=h_max
        @species S(t) I(t) R(t)
        eqs_ordered = [
            Reaction(i, [S, I], [I], [1, 1], [2]),
            Reaction(r, [I], [R]),
            D(M) ~ -I*M/(m1 + m2),
            H ~ h_max - I
        ]
        @named hybrid_sir_ordered = ReactionSystem(eqs_ordered, t)
        hybrid_sir_ordered = complete(hybrid_sir_ordered)

        # Checks that ODE an simulation of the system achieves the correct steady state.
        oprob_ordered = ODEProblem(hybrid_sir_ordered, u0, tspan, ps; structural_simplify = true)
        solve(oprob_ordered, Vern7(); abstol = 1e-8, reltol = 1e-8, saveat = 1.0)
    end

    # Declares the model in a messy fashion, and simulates it.
    osol_messy = let
        @variables τ M(τ) H(τ)=h_max
        @species S(τ) I(τ) R(τ)
        Δ = Differential(τ)
        eqs_messy = [
            Reaction(i, [S, I], [I], [1, 1], [2]),
            Reaction(r, [I], [R]),
            I*M + m1*Δ(M) ~ -m2*Δ(M),
            H ~ h_max - I
        ]
        @named hybrid_sir_messy = ReactionSystem(eqs_messy, τ)
        hybrid_sir_messy = complete(hybrid_sir_messy)

        # Checks that ODE an simulation of the system achieves the correct steady state.
        oprob_messy = ODEProblem(hybrid_sir_messy, u0, tspan, ps; structural_simplify = true)
        solve(oprob_messy, Vern7(); abstol = 1e-8, reltol = 1e-8, saveat = 1.0)
    end

    # Checks that the simulations are identical.
    # Some internal details will be different, however, the solutions should be identical.
    osol_messy[[:S, :I, :R, :M, :H]] ≈ osol_ordered[[:S, :I, :R, :M, :H]]
end

### DSL Tests ###

### Error Tests ###

# Checks that various erroneous hybrid system declarations yield errors.
let 
    @parameters p1 p2
    @variables τ  U1(τ) V1(t)
    @species R1(τ) R2(τ) S1(t) S2(t) 
    E = Differential(τ)

    # Variables as reaction reactants.
    @test_throws Exception ReactionSystem([
        Reaction(p1, [S1], [V1])
    ], t; name = :rs)

    # Species using non-declared independent variable.
    @test_throws Exception ReactionSystem([
        Reaction(p1, [R1], [R2])
    ], t; name = :rs)

    # Equation with variable using non-declared independent variable. 
    @test_throws Exception ReactionSystem([
        Reaction(p1, [S1], [S2]),
        U1 ~ S1 + p2
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

    # Hybrid system with additional differential equation for species.
    eqs = [
        Reaction(p1, [S1], [S2]),
        D(S1) ~ p2 - S1
    ]
    @named rs = ReactionSystem(eqs, t)
    rs = complete(rs)
    u0 = [S1 => 1.0, S2 => 2.0]
    ps = [p1 => 2.0, p2 => 3.0]
    @test_throws Exception ODEProblem(rs, u0, (0.0, 1.0), ps; structural_simplify = true)

    # Hybrid system overconstrained due to additional algebraic equations (without variables).
    eqs = [
        Reaction(p1, [S1], [S2]),
        S1 ~ p2 + S1,
    ]
    @named rs = ReactionSystem(eqs, t)
    rs = complete(rs)
    u0 = [S1 => 1.0, S2 => 2.0]
    ps = [p1 => 2.0, p2 => 3.0]
    @test_throws Exception ODEProblem(rs, u0, (0.0, 1.0), ps; structural_simplify = true)

    # Hybrid system overconstrained due to additional algebraic equations (with variables).
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