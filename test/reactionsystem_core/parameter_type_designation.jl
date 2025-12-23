### Fetch Packages and Set Global Variables ###

# Fetch packages.
using Catalyst, JumpProcesses, NonlinearSolve, OrdinaryDiffEqTsit5, StochasticDiffEq, Test
using Symbolics: BasicSymbolic, unwrap

# Sets stable rng number.
using StableRNGs
rng = StableRNG(12345)
seed = rand(rng, 1:100)

# Sets the default `t` to use.
t = default_t()

### Basic Tests ###

# Declares a simple model to run tests on.
begin
    t = default_t()
    @parameters p1 p2 p3::Int64 p4 p5::Rational{Int64}
    @parameters d1 d2 = 1.2 d3::Int64 = 2 [description = "A parameter"] d4::Rational{Int64} d5
    @species X1(t) X2(t) X3(t) X4(t) X5(t)

    rxs = [
        Reaction(p1, nothing, [X1]),
        Reaction(p2, nothing, [X2]),
        Reaction(p3, nothing, [X3]),
        Reaction(p4, nothing, [X4]),
        Reaction(p5, nothing, [X5]),
        Reaction(d1, [X1], nothing),
        Reaction(d2, [X2], nothing),
        Reaction(d3, [X3], nothing),
        Reaction(d4, [X4], nothing),
        Reaction(d5, [X5], nothing)
    ]
    @named rs = ReactionSystem(rxs)
    rs = complete(rs)

    # Declares initial condition and potential parameter sets.
    u0 = [X1 => 0.1, X2 => 0.2, X3 => 0.3, X4 => 0.4, X5 => 0.5]
    p_alts = [
        [p1 => 1.0, d1 => 1.0, p2 => 1.2, p3 => 2, p4 => 0.5, d4 => 1//2, p5 => 3//2, d5 => 1.5],
        (p1 => 1.0, d1 => 1.0, p2 => 1.2, p3 => 2, p4 => 0.5, d4 => 1//2, p5 => 3//2, d5 => 1.5),
        Dict([p1 => 1.0, d1 => 1.0, p2 => 1.2, p3 => 2, p4 => 0.5, d4 => 1//2, p5 => 3//2, d5 => 1.5])
    ]
end

# Tests that parameters stored in the system have the correct type.
let
    @test SymbolicUtils.symtype(rs.p1) == Real
    @test SymbolicUtils.symtype(rs.d1) == Real
    @test SymbolicUtils.symtype(rs.p2) == Real
    @test SymbolicUtils.symtype(rs.d2) == Real
    @test SymbolicUtils.symtype(rs.p3) == Int64
    @test SymbolicUtils.symtype(rs.d3) == Int64
    @test SymbolicUtils.symtype(rs.p4) == Real
    @test SymbolicUtils.symtype(rs.d4) == Rational{Int64}
    @test SymbolicUtils.symtype(rs.p5) == Rational{Int64}
    @test SymbolicUtils.symtype(rs.d5) == Real
end

# Tests that simulations with differentially typed variables yields correct results.
let
    for p in p_alts
        oprob = ODEProblem(rs, u0, (0.0, 1000.0), p; abstol = 1e-10, reltol = 1e-10)
        sol = solve(oprob, Tsit5())
        @test all(sol.u[end] .≈ 1.0)
    end
end

# Test that the various structures stores the parameters using the correct type.
let
    # Creates problems, integrators, and solutions.
    oprob = ODEProblem(rs, u0, (0.0, 1.0), p_alts[1])
    sprob = SDEProblem(rs, u0, (0.0, 1.0), p_alts[1])
    jprob = JumpProblem(rs, u0, (0.0, 1.0), p_alts[1]; rng)
    nprob = NonlinearProblem(rs, u0, p_alts[1])

    oinit = init(oprob, Tsit5())
    sinit = init(sprob, ImplicitEM())
    jinit = init(jprob, SSAStepper())
    ninit = init(nprob, NewtonRaphson())

    osol = solve(oprob, Tsit5())
    ssol = solve(sprob, ImplicitEM(); seed)
    jsol = solve(jprob, SSAStepper(); seed)
    nsol = solve(nprob, NewtonRaphson())

    # Checks all stored parameters.
    for mtk_struct in [oprob, sprob, jprob, nprob, oinit, sinit, jinit, ninit, osol, ssol, jsol, nsol]
        # Checks that all parameters have the correct type.
        @test unwrap(mtk_struct.ps[p1]) isa Float64
        @test unwrap(mtk_struct.ps[d1]) isa Float64
        @test unwrap(mtk_struct.ps[p2]) isa Float64
        @test unwrap(mtk_struct.ps[d2]) isa Float64
        @test unwrap(mtk_struct.ps[p3]) isa Int64
        @test unwrap(mtk_struct.ps[d3]) isa Int64
        @test unwrap(mtk_struct.ps[p4]) isa Float64
        @test unwrap(mtk_struct.ps[d4]) isa Rational{Int64}
        @test unwrap(mtk_struct.ps[p5]) isa Rational{Int64}
        @test unwrap(mtk_struct.ps[d5]) isa Float64

        # Checks that all parameters have the correct value.
        @test unwrap(mtk_struct.ps[p1]) == 1.0
        @test unwrap(mtk_struct.ps[d1]) == 1.0
        @test unwrap(mtk_struct.ps[p2]) == 1.2
        @test unwrap(mtk_struct.ps[d2]) == 1.2
        @test unwrap(mtk_struct.ps[p3]) == 2
        @test unwrap(mtk_struct.ps[d3]) == 2
        @test unwrap(mtk_struct.ps[p4]) == Float32(0.5)
        @test unwrap(mtk_struct.ps[d4]) == 1//2
        @test unwrap(mtk_struct.ps[p5]) == 3//2
        @test unwrap(mtk_struct.ps[d5]) == Float32(1.5)
    end

    # Checks all stored variables (these should always be `Float64`).
    for mtk_struct in [oprob, sprob, jprob, nprob, oinit, sinit, jinit, ninit]
        # Checks that all variables have the correct type.
        @test unwrap(mtk_struct[X1]) isa Float64
        @test unwrap(mtk_struct[X2]) isa Float64
        @test unwrap(mtk_struct[X3]) isa Float64
        @test unwrap(mtk_struct[X4]) isa Float64
        @test unwrap(mtk_struct[X5]) isa Float64

        # Checks that all variables have the correct value.
        @test unwrap(mtk_struct[X1]) == 0.1
        @test unwrap(mtk_struct[X2]) == 0.2
        @test unwrap(mtk_struct[X3]) == 0.3
        @test unwrap(mtk_struct[X4]) == 0.4
        @test unwrap(mtk_struct[X5]) == 0.5
    end
end
