### Prepares Tests ###

# Fetch packages.
using Catalyst, JumpProcesses, OrdinaryDiffEq, StableRNGs, Test
using ModelingToolkitBase: Pre
const MT = ModelingToolkitBase

# Sets the default `t` and `D` to use, and a stable rng.
t = default_t()
D = default_time_deriv()
rng = StableRNG(12345)

# Whether `x` is an (unscalarized) array symbolic.
is_array_symbolic(x) = Symbolics.scalarize(Symbolics.unwrap(x)) isa AbstractArray

# Array unknowns are scalarized whichever way a `ReactionSystem` is built (the DSL, the
# constructor inferring the unknowns, or the explicit constructor). ModelingToolkitBase doesn't
# handle unscalarized array unknowns correctly in jumps and events, and reactions refer to array
# species through their elements (SciML/ModelingToolkit.jl#5428, SciML/Catalyst.jl#1607). The
# array species and variables sit between scalar ones, so that a wrong position in the state
# would hit a different unknown.

### Reaction Networks With Equations and Events ###

# The same model, built in the three ways and with scalar species and variables instead (the
# reference). The discrete event changes `B`, after the array species `X`; the continuous event's
# condition reads `W`, after the array variable `V`.
@species A(t) X(t)[1:2] B(t)
@variables V(t)[1:2] W(t)
@parameters k1 k2 k3 d
rxs = [Reaction(k1, [A], [X[1]]), Reaction(k2, [X[1]], [X[2]]), Reaction(k3, [X[2]], [B]),
    Reaction(d, [B], nothing)]
eqs = [D(V[1]) ~ -V[1], D(V[2]) ~ -2V[2], D(W) ~ 1]
devs = [[1.0] => [B ~ Pre(B) + 10]]
cevs = [[W ~ 1.5] => [A ~ Pre(A) + 5]]
rs_explicit = complete(ReactionSystem([rxs; eqs], t, [A, X, B, V, W], [k1, k2, k3, d];
    discrete_events = devs, continuous_events = cevs, name = :rs))
rs_inferred = complete(ReactionSystem([rxs; eqs], t; discrete_events = devs,
    continuous_events = cevs, name = :rs))
rs_dsl = @reaction_network rs begin
    @species A(t) X(t)[1:2] B(t)
    @variables V(t)[1:2] W(t)
    @equations begin
        D(V[1]) ~ -V[1]
        D(V[2]) ~ -2V[2]
        D(W) ~ 1
    end
    @discrete_events [1.0] => [B => B + 10]
    @continuous_events [W ~ 1.5] => [A => A + 5]
    k1, A --> X[1]
    k2, X[1] --> X[2]
    k3, X[2] --> B
    d, B --> 0
end

@species X1(t) X2(t)
@variables V1(t) V2(t)
rs_ref = complete(ReactionSystem(
    [Reaction(k1, [A], [X1]), Reaction(k2, [X1], [X2]), Reaction(k3, [X2], [B]),
        Reaction(d, [B], nothing), D(V1) ~ -V1, D(V2) ~ -2V2, D(W) ~ 1], t,
    [A, X1, X2, B, V1, V2, W], [k1, k2, k3, d];
    discrete_events = [[1.0] => [B ~ Pre(B) + 10]],
    continuous_events = [[W ~ 1.5] => [A ~ Pre(A) + 5]], name = :rs_ref))
to_ref = Dict(A => A, X[1] => X1, X[2] => X2, B => B, V[1] => V1, V[2] => V2, W => W)

let
    ps = [k1 => 1.0, k2 => 2.0, k3 => 3.0, d => 0.5]
    u0_ref = [A => 4.0, X1 => 1.0, X2 => 2.0, B => 3.0, V1 => 1.0, V2 => 2.0, W => 0.0]
    sol_ref = solve(ODEProblem(rs_ref, u0_ref, (0.0, 3.0), ps), Tsit5(); abstol = 1e-10,
        reltol = 1e-10)
    # The explicit constructor keeps the given order, with the species first.
    @test isequal(unknowns(rs_explicit), [A, X[1], X[2], B, V[1], V[2], W])
    rs_interleaved = ReactionSystem([rxs; eqs], t, [V, A, W, X, B], [k1, k2, k3, d]; name = :rs)
    @test isequal(unknowns(rs_interleaved), [A, X[1], X[2], B, V[1], V[2], W])
    for rs in (rs_explicit, rs_inferred, rs_dsl)
        # The unknowns are scalarized, with the species first.
        @test !any(is_array_symbolic, unknowns(rs))
        @test issetequal(unknowns(rs), [A, X[1], X[2], B, V[1], V[2], W])
        @test length(unknowns(rs)) == 7
        @test issetequal(species(rs), [A, X[1], X[2], B])
        @test isequal(unknowns(rs)[1:numspecies(rs)], species(rs))
        @test numspecies(rs) == 4
        @test !any(is_array_symbolic, unknowns(ode_model(rs)))

        # The stoichiometry, ODEs (and events) are those of the reference.
        perm = [findfirst(isequal(to_ref[sp]), species(rs_ref)) for sp in species(rs)]
        @test netstoichmat(rs) == netstoichmat(rs_ref)[perm, :]
        u0 = [A => 4.0, X => [1.0, 2.0], B => 3.0, V => [1.0, 2.0], W => 0.0]
        sol = solve(ODEProblem(rs, u0, (0.0, 3.0), ps), Tsit5(); abstol = 1e-10, reltol = 1e-10)
        for u in unknowns(rs)
            @test sol[u][end] ≈ sol_ref[to_ref[u]][end]
        end
        @test sol[rs.X][end] ≈ [sol_ref[X1][end], sol_ref[X2][end]]
    end
end

# Equations for whole arrays are not supported, and give an error whichever way the system is
# built; one equation per element (here from `Symbolics.scalarize`) works.
let
    rx = Reaction(k1, nothing, [A])
    @test_throws "for a whole array" ReactionSystem([rx, D(V) ~ -V], t; name = :rs)
    @test_throws "for a whole array" ReactionSystem([rx, D(V) ~ -V], t, [A, V], [k1]; name = :rs)
    @test_throws "for a whole array" @reaction_network begin
        @variables V(t)[1:2]
        @equations D(V) ~ -V
        k1, 0 --> A
    end

    rs = complete(ReactionSystem([rx; Symbolics.scalarize(D(V) ~ -V); D(W) ~ 1], t;
        discrete_events = [[1.0] => [W ~ Pre(W) + 10]], name = :rs))
    @test issetequal(unknowns(rs), [A, V[1], V[2], W])
    sol = solve(ODEProblem(rs, [A => 0.0, V => [1.0, 2.0], W => 0.0], (0.0, 2.0), [k1 => 1.0]),
        Tsit5(); abstol = 1e-10, reltol = 1e-10)
    @test sol[rs.V][end] ≈ [exp(-2.0), 2exp(-2.0)]
    @test sol[rs.W][end] ≈ 12.0
end

### Jump Problems ###

# Both reactions consume `X[2]`, which starts at 1, so only one of them can occur: every
# aggregator must update the other reaction's rate when `X[2]` changes (the dependency maps
# indexed `unknowns` instead of the state, SciML/ModelingToolkit.jl#5428). The non-reaction jump
# changes `B`, after the array species `X`. `X[1]` takes part in a reaction that never occurs, so
# that the constructor inferring the unknowns includes it too.
let
    @parameters k k0 kB
    jrxs = [Reaction(k, [X[2]], [A]), Reaction(k, [X[2]], nothing), Reaction(k0, nothing, [X[1]])]
    crj = ConstantRateJump(kB, [B ~ Pre(B) + 1])
    aggregators = (Direct(), NRM(), SortingDirect(), RDirect(), DirectCR(), RSSA(), RSSACR())

    # The states the system can be in (`X[2]` and `A`): the dependency maps must not let both
    # reactions occur, and a jump mustn't change the wrong species.
    valid_states = ((1, 0), (0, 0), (0, 1))

    # Reactions only (mass action jumps), with the explicit constructor: the reproducer of
    # SciML/ModelingToolkit.jl#5428.
    rs = complete(ReactionSystem(jrxs[1:2], t, [X, A], [k]; name = :rs))
    @test !any(is_array_symbolic, unknowns(jump_model(rs)))
    for aggregator in aggregators
        jprob = JumpProblem(rs, [X => [20, 1], A => 0], (0.0, 20.0), [k => 1.0]; aggregator)
        nconsumed = 0
        for _ in 1:10
            sol = solve(jprob, SSAStepper(); rng)
            @test all(in(valid_states), zip(sol[rs.X[2]], sol[rs.A]))
            @test all(==(20), sol[rs.X[1]])
            nconsumed += (sol[rs.X[2]][end] == 0)
        end
        @test nconsumed > 0
    end

    # With a non-reaction jump as well, with the explicit and the inferring constructors.
    jnets = [
        ReactionSystem(jrxs, t, [X, A, B], [k, k0, kB]; jumps = [crj], name = :rs),
        ReactionSystem(jrxs, t; jumps = [crj], name = :rs),
    ]
    for rs in complete.(jnets), aggregator in aggregators
        @test !any(is_array_symbolic, unknowns(rs))
        jprob = JumpProblem(rs, [X => [20, 1], A => 0, B => 0], (0.0, 20.0),
            [k => 1.0, k0 => 0.0, kB => 5.0]; aggregator)
        nconsumed = 0
        for _ in 1:10
            sol = solve(jprob, SSAStepper(); rng)
            @test all(in(valid_states), zip(sol[rs.X[2]], sol[rs.A]))
            @test all(==(20), sol[rs.X[1]])
            @test sol[rs.B][end] > 0 # (the jump fires at rate 5 over time 20)
            nconsumed += (sol[rs.X[2]][end] == 0)
        end
        @test nconsumed > 0
    end

    # The same reactions with the DSL (which has no option for non-reaction jumps).
    rs = @reaction_network begin
        @species X(t)[1:2] A(t)
        k, X[2] --> A
        k, X[2] --> 0
        k0, 0 --> X[1]
    end
    @test !any(is_array_symbolic, unknowns(rs))
    for aggregator in aggregators
        jprob = JumpProblem(rs, [:X => [20, 1], :A => 0], (0.0, 20.0), [:k => 1.0, :k0 => 0.0];
            aggregator)
        nconsumed = 0
        for _ in 1:10
            sol = solve(jprob, SSAStepper(); rng)
            @test all(in(valid_states), zip(sol[rs.X[2]], sol[rs.A]))
            @test all(==(20), sol[rs.X[1]])
            nconsumed += (sol[rs.X[2]][end] == 0)
        end
        @test nconsumed > 0
    end
end

### Other Inputs ###

# Array species can be given values as arrays or by element, and are read in both ways.
let
    @species Y(t)[1:3]
    @parameters p
    rs = complete(ReactionSystem([Reaction(p, [Y[1]], [Y[2]]), Reaction(p, [Y[2]], [Y[3]])], t,
        [Y], [p]; name = :rs))
    @test isequal(species(rs), [Y[1], Y[2], Y[3]])
    prob = ODEProblem(rs, [:Y => [1.0, 2.0, 3.0]], 1.0, [:p => 1.0])
    @test prob[rs.Y] == [1.0, 2.0, 3.0]
    @test prob[rs.Y[2]] == 2.0
    prob2 = ODEProblem(rs, [Y[1] => 1.0, Y[2] => 2.0, Y[3] => 3.0], 1.0, [p => 1.0])
    @test prob2[rs.Y] == [1.0, 2.0, 3.0]

    # Listing an element as well as its array gives the element once.
    rs2 = ReactionSystem([Reaction(p, [Y[1]], [Y[2]])], t, [Y, Y[1]], [p]; name = :rs2)
    @test isequal(species(rs2), [Y[1], Y[2], Y[3]])
end

# The inferring constructor includes only the elements the reactions use; listing the array in the
# explicit constructor includes all of them.
let
    @parameters p
    @test isequal(unknowns(ReactionSystem([Reaction(p, [X[2]], nothing)], t; name = :rs)), [X[2]])
    @test isequal(unknowns(ReactionSystem([Reaction(p, [X[2]], nothing)], t, [X], [p]; name = :rs)),
        [X[1], X[2]])
end

# Matrices of species are scalarized in column-major order, and can be given and read as matrices.
let
    @species C(t) M(t)[1:2, 1:2]
    @parameters p
    mrxs = [Reaction(p, [M[1, 2]], [C])]
    rs_dsl = @reaction_network rs begin
        @species C(t) M(t)[1:2, 1:2]
        p, M[1, 2] --> C
    end
    rs_explicit = complete(ReactionSystem(mrxs, t, [C, M], [p]; name = :rs))
    @test isequal(unknowns(rs_explicit), [C, M[1, 1], M[2, 1], M[1, 2], M[2, 2]])
    @test isequal(unknowns(rs_dsl), [C, M[1, 1], M[2, 1], M[1, 2], M[2, 2]])
    for rs in (rs_explicit, rs_dsl)
        prob = ODEProblem(rs, [C => 0.0, M => [1.0 2.0; 3.0 4.0]], 1.0, [p => 1.0])
        @test prob[rs.M] == [1.0 2.0; 3.0 4.0]
        @test prob[rs.M[2, 1]] == 3.0
        sol = solve(prob, Tsit5(); abstol = 1e-10, reltol = 1e-10)
        @test sol[rs.C][end] ≈ 2.0 * (1 - exp(-1.0))
        @test sol[rs.M][end] ≈ [1.0 2exp(-1.0); 3.0 4.0]
    end
end

# A default value declared for an array variable is used when no value is given.
let
    @variables U(t)[1:2] = [1.0, 2.0]
    @parameters p
    rs = complete(ReactionSystem([Reaction(p, nothing, [A]), D(U[1]) ~ -U[1], D(U[2]) ~ -U[2]], t,
        [A, U], [p]; name = :rs))
    @test isequal(unknowns(rs), [A, U[1], U[2]])
    sol = solve(ODEProblem(rs, [A => 0.0], 1.0, [p => 1.0]), Tsit5(); abstol = 1e-10, reltol = 1e-10)
    @test sol[rs.U][end] ≈ [exp(-1.0), 2exp(-1.0)]
end

# A subsystem using the elements of an array variable of its parent (through `ParentScope`): the
# completed (flattened) system, and the conversions of one completed without flattening, have them.
let
    @parameters p
    V1p, V2p = ParentScope(V[1]), ParentScope(V[2])
    @named child = ReactionSystem([Reaction(p, nothing, [A]), D(V1p) ~ -V1p, D(V2p) ~ -2V2p], t)
    @named parent = ReactionSystem(Reaction[], t)
    rs = complete(compose(parent, [child]))
    @test !any(is_array_symbolic, unknowns(rs))
    @test issetequal(unknowns(rs), [rs.child.A, V[1], V[2]])
    hier = complete(compose(parent, [child]); flatten = false)
    @test !any(is_array_symbolic, unknowns(ode_model(hier)))
    @test issetequal(unknowns(ode_model(hier)), [hier.child.A, V[1], V[2]])
    sol = solve(ODEProblem(rs, [V => [1.0, 2.0], rs.child.A => 0.0], 1.0, [rs.child.p => 1.0]),
        Tsit5(); abstol = 1e-10, reltol = 1e-10)
    @test sol[V][end] ≈ [exp(-1.0), 2exp(-2.0)]
end

# A repeated unknown is listed once.
let
    @parameters p
    rs = complete(ReactionSystem([Reaction(p, [A], nothing)], t, [A, A], [p]; name = :rs))
    @test isequal(unknowns(rs), [A])
    sol = solve(ODEProblem(rs, [A => 1.0], 1.0, [p => 1.0]), Tsit5(); abstol = 1e-10, reltol = 1e-10)
    @test sol[rs.A][end] ≈ exp(-1.0)
end

# Array elements given as unknowns are kept as they are (also with indices not starting at 1),
# while a whole array gives all its elements.
let
    @species Z(t)[0:2]
    @parameters p
    rx = Reaction(p, [Z[0]], [Z[2]])
    @test isequal(unknowns(ReactionSystem([rx], t, [Z[0], Z[2]], [p]; name = :rs)), [Z[0], Z[2]])
    @test isequal(unknowns(ReactionSystem([rx], t, [Z], [p]; name = :rs)), [Z[0], Z[1], Z[2]])
end

# The elements of an array of boundary condition species are boundary condition species.
let
    @species Q(t)[1:2] [isbcspecies = true]
    @parameters p
    @test Catalyst.isbc(Q[1]) && Catalyst.isbc(Q[2])
    rs_explicit = complete(ReactionSystem([Reaction(p, [Q[1]], [Q[1], A]), Q[1] ~ 3, Q[2] ~ 4], t,
        [Q, A], [p]; name = :rs))
    rs_dsl = @reaction_network rs begin
        @species Q(t)[1:2] [isbcspecies = true] A(t)
        @equations begin
            Q[1] ~ 3
            Q[2] ~ 4
        end
        p, Q[1] --> Q[1] + A
    end
    for rs in (rs_explicit, rs_dsl)
        @test issetequal(unknowns(rs), [Q[1], Q[2], A])
        sol = solve(ODEProblem(rs, [A => 0.0], 1.0, [p => 1.0]; mtkcompile = true), Tsit5())
        @test sol[rs.A][end] ≈ 3.0
        @test sol[rs.Q][end] ≈ [3.0, 4.0]
    end
end
