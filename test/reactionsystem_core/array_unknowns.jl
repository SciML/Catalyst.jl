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
    for rs in (rs_explicit, rs_inferred, rs_dsl)
        # The unknowns are scalarized, with the species first.
        @test !any(is_array_symbolic, unknowns(rs))
        @test issetequal(unknowns(rs), [A, X[1], X[2], B, V[1], V[2], W])
        @test issetequal(species(rs), [A, X[1], X[2], B])
        @test numspecies(rs) == 4

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

# An equation for a whole array variable also gives scalarized unknowns.
let
    rs = complete(ReactionSystem([Reaction(k1, nothing, [A]), D(V) ~ -V, D(W) ~ 1], t;
        discrete_events = [[1.0] => [W ~ Pre(W) + 10]], name = :rs))
    @test !any(is_array_symbolic, unknowns(rs))
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
    jnets = [
        ReactionSystem(jrxs, t, [X, A, B], [k, k0, kB]; jumps = [crj], name = :rs),
        ReactionSystem(jrxs, t; jumps = [crj], name = :rs),
    ]
    aggregators = (Direct(), NRM(), SortingDirect(), RDirect(), DirectCR(), RSSA(), RSSACR())
    for rs in complete.(jnets), aggregator in aggregators
        @test !any(is_array_symbolic, unknowns(rs))
        jprob = JumpProblem(rs, [X => [20, 1], A => 0, B => 0], (0.0, 20.0),
            [k => 1.0, k0 => 0.0, kB => 1.0]; aggregator)
        for _ in 1:10
            sol = solve(jprob, SSAStepper(); rng)
            @test minimum(sol[rs.X[2]]) == 0
            @test sol[rs.X[2]][end] == 0
            @test sol[rs.A][end] in (0, 1)
            @test all(==(20), sol[rs.X[1]])
            @test sol[rs.B][end] > 0
        end
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
        for _ in 1:10
            sol = solve(jprob, SSAStepper(); rng)
            @test minimum(sol[rs.X[2]]) == 0
            @test sol[rs.A][end] in (0, 1)
            @test all(==(20), sol[rs.X[1]])
        end
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
