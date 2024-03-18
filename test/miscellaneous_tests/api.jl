#! format: off

### Fetch Packages and Test Networks ###
using Catalyst, DiffEqBase, ModelingToolkit, Test, OrdinaryDiffEq, NonlinearSolve
using StochasticDiffEq
using LinearAlgebra: norm
using SparseArrays
using ModelingToolkit: value
t = default_t()

include("../test_networks.jl")

### Base Tests ###

let
    @parameters k1 k2
    @species S(t) I(t) R(t)
    rxs = [Reaction(k1, [S, I], [I], [1, 1], [2]),
        Reaction(k2, [I], [R])]
    @named rs = ReactionSystem(rxs, t, [S, I, R], [k1, k2])

    specset = Set([value(S) => 1, value(I) => 2, value(R) => 3])
    @test issetequal(specset, speciesmap(rs))

    pset = Set([value(k1) => 1, value(k2) => 2])
    @test issetequal(pset, paramsmap(rs))

    rxs = [Reaction(k1, [S, I], [I], [1, 1], [2]),
        Reaction(k2, [I], [R])]
    @named rs = ReactionSystem(rxs, t, [S, I, R], [k1, k2])
    rxs2 = [Reaction(k2, [I], [R], [1], [1]),
        Reaction(k1, [S, I], [I], [1, 1], [2])]
    rs2 = ReactionSystem(rxs2, t, [R, I, S], [k2, k1], name = :rs)
    @test rs == rs2

    rs3 = make_empty_network()
    @parameters k3 k4
    @species D(t)
    addspecies!(rs3, S)
    addspecies!(rs3, D)
    addparam!(rs3, k3)
    addparam!(rs3, k4)
    @test issetequal(species(rs3), [S, D])
    @test issetequal(parameters(rs3), [k3, k4])
    addreaction!(rs3, Reaction(k3, [S], [D]))
    addreaction!(rs3, Reaction(k4, [S, I], [D]))
    merge!(rs, rs3)
    addspecies!(rs2, S)
    addspecies!(rs2, D)
    addparam!(rs2, k3)
    addparam!(rs2, k4)
    addreaction!(rs2, Reaction(k3, [S], [D]))
    addreaction!(rs2, Reaction(k4, [S, I], [D]))
    @test rs2 == rs

    rxs = [Reaction(k1, [S, I], [I], [1, 1], [2]),
        Reaction(k2, [I], [R])]
    @named rs = ReactionSystem(rxs, t, [S, I, R], [k1, k2])
    rs3 = make_empty_network()
    addspecies!(rs3, S)
    addspecies!(rs3, D)
    addparam!(rs3, k3)
    addparam!(rs3, k4)
    addreaction!(rs3, Reaction(k3, [S], [D]))
    addreaction!(rs3, Reaction(k4, [S, I], [D]))
    rs4 = extend(rs, rs3)
    @test rs2 == rs4

    rxs = [Reaction(k1 * S, [S, I], [I], [2, 3], [2]),
        Reaction(k2 * R, [I], [R])]
    @named rs = ReactionSystem(rxs, t, [S, I, R], [k1, k2])
    deps = dependents(rxs[2], rs)
    @test isequal(deps, [R, I])
    @test isequal(dependents(rxs[1], rs), dependants(rxs[1], rs))
    addspecies!(rs, S)
    @test numspecies(rs) == 3
    addspecies!(rs, S, disablechecks = true)
    @test numspecies(rs) == 4
    addparam!(rs, k1)
    @test numparams(rs) == 2
    @test numreactionparams(rs) == 2
    addparam!(rs, k1, disablechecks = true)
    @test numparams(rs) == 3
    @test numreactionparams(rs) == 3
end

let
    rnmat = @reaction_network begin
        α, S + 2I --> 2I
        β, 3I --> 2R + S
    end

    smat = [1 0;
            2 3;
            0 0]
    pmat = [0 1;
            2 0;
            0 2]
    @test smat == substoichmat(rnmat) == Matrix(substoichmat(rnmat, sparse = true))
    @test pmat == prodstoichmat(rnmat) == Matrix(prodstoichmat(rnmat, sparse = true))
end


### Test Intermediate Complexes Reaction Networks ###

# Test function.
function testnetwork(rn, B, Z, Δ, lcs, d, subrn, lcd; skiprxtest = false)
    B2 = reactioncomplexes(rn)[2]
    @test B == B2 == Matrix(reactioncomplexes(rn, sparse = true)[2])
    @test B == incidencemat(rn)
    @test Z == complexstoichmat(rn) == Matrix(complexstoichmat(rn, sparse = true))
    @test Δ == complexoutgoingmat(rn) == Matrix(complexoutgoingmat(rn, sparse = true))
    ig = incidencematgraph(rn)
    lcs2 = linkageclasses(rn)
    @test lcs2 == linkageclasses(incidencematgraph(sparse(B))) == lcs
    @test deficiency(rn) == d
    if !skiprxtest
        @test all(issetequal.(subrn, reactions.(subnetworks(rn))))
    end
    @test linkagedeficiencies(rn) == lcd
    @test sum(linkagedeficiencies(rn)) <= deficiency(rn)
end


# Mass-action non-catalytic.
let
    rn = @reaction_network begin
        k₁, 2A --> B
        k₂, A --> C
        k₃, C --> D
        k₄, B + D --> E
    end
    Z = [2 0 1 0 0 0 0;
        0 1 0 0 0 1 0;
        0 0 0 1 0 0 0;
        0 0 0 0 1 1 0;
        0 0 0 0 0 0 1]
    B = [-1 0 0 0;
        1 0 0 0;
        0 -1 0 0;
        0 1 -1 0;
        0 0 1 0;
        0 0 0 -1;
        0 0 0 1]
    Δ = [-1 0 0 0; 0 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 0; 0 0 0 -1; 0 0 0 0]
    lcs = [[1, 2], [3, 4, 5], [6, 7]]
    r = reactions(rn)
    subrn = [[r[1]], [r[2], r[3]], [r[4]]]
    lcd = [0, 0, 0]
    testnetwork(rn, B, Z, Δ, lcs, 0, subrn, lcd)

    # Constant and BC species test.
    @parameters F [isconstantspecies = true]
    crn = @reaction_network begin
        k₁, 2A --> B
        k₂, A --> C + $F
        k₃, C --> D
        k₄, B + D + $F --> E
    end
    testnetwork(crn, B, Z, Δ, lcs, 0, subrn, lcd; skiprxtest = true)
end

# Mass-action rober.
let
    rn = @reaction_network begin
        k₁, A --> B
        k₂, B + B --> C + B
        k₃, B + C --> A + C
    end
    Z = [1 0 0 0 1;
        0 1 2 1 0;
        0 0 0 1 1]
    B = [-1 0 0;
        1 0 0;
        0 -1 0;
        0 1 -1;
        0 0 1]
    Δ = [-1 0 0; 0 0 0; 0 -1 0; 0 0 -1; 0 0 0]
    lcs = [[1, 2], [3, 4, 5]]
    r = reactions(rn)
    subrn = [[r[1]], [r[2], r[3]]]
    lcd = [0, 0]
    testnetwork(rn, B, Z, Δ, lcs, 1, subrn, lcd)
end

# Some rational functions as rates.
let
    rn = @reaction_network begin
        k₁, ∅ --> X₁
        (k₂ / (1 + X₁ * X₂ + X₃ * X₄), k₃ / (1 + X₁ * X₂ + X₃ * X₄)), 2X₁ + X₂ ↔ 3X₃ + X₄
        k₄, X₄ --> ∅
    end
    Z = [0 1 2 0 0;
        0 0 1 0 0;
        0 0 0 3 0;
        0 0 0 1 1]
    B = [-1 0 0 1;
        1 0 0 0;
        0 -1 1 0;
        0 1 -1 0;
        0 0 0 -1]
    Δ = [-1 0 0 0; 0 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 -1]
    lcs = [[1, 2, 5], [3, 4]]
    r = reactions(rn)
    subrn = [[r[1], r[4]], [r[2], r[3]]]
    lcd = [0, 0]
    testnetwork(rn, B, Z, Δ, lcs, 0, subrn, lcd)
end

# Repressilator.
let
    rn = @reaction_network begin
        hillr(P₃, α, K, n), ∅ --> m₁
        hillr(P₁, α, K, n), ∅ --> m₂
        hillr(P₂, α, K, n), ∅ --> m₃
        (δ, γ), m₁ ↔ ∅
        (δ, γ), m₂ ↔ ∅
        (δ, γ), m₃ ↔ ∅
        β, m₁ --> m₁ + P₁
        β, m₂ --> m₂ + P₂
        β, m₃ --> m₃ + P₃
        μ, P₁ --> ∅
        μ, P₂ --> ∅
        μ, P₃ --> ∅
    end
    Z = [0 1 0 0 1 0 0 0 0 0;
        0 0 1 0 0 1 0 0 0 0;
        0 0 0 1 0 0 1 0 0 0;
        0 0 0 0 1 0 0 1 0 0;
        0 0 0 0 0 1 0 0 1 0;
        0 0 0 0 0 0 1 0 0 1]
    B = [-1 -1 -1 1 -1 1 -1 1 -1 0 0 0 1 1 1;
        1 0 0 -1 1 0 0 0 0 -1 0 0 0 0 0;
        0 1 0 0 0 -1 1 0 0 0 -1 0 0 0 0;
        0 0 1 0 0 0 0 -1 1 0 0 -1 0 0 0;
        0 0 0 0 0 0 0 0 0 1 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 1 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1]
    Δ = [-1 -1 -1 0 -1 0 -1 0 -1 0 0 0 0 0 0; 0 0 0 -1 0 0 0 0 0 -1 0 0 0 0 0;
        0 0 0 0 0 -1 0 0 0 0 -1 0 0 0 0;
        0 0 0 0 0 0 0 -1 0 0 0 -1 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1]
    lcs = [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]]
    r = reactions(rn)
    subrn = [[r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8], r[9], r[13], r[14], r[15],
        r[10], r[11], r[12]]]
    lcd = [3]
    testnetwork(rn, B, Z, Δ, lcs, 3, subrn, lcd)
end

# Brusselator.
let
    rn = @reaction_network begin
        A, ∅ → X
        1, 2X + Y → 3X
        B, X → Y
        1, X → ∅
    end
    Z = [0 1 2 3 0;
        0 0 1 0 1]
    B = [-1 0 0 1;
        1 0 -1 -1;
        0 -1 0 0;
        0 1 0 0;
        0 0 1 0]
    Δ = [-1 0 0 0; 0 0 -1 -1; 0 -1 0 0; 0 0 0 0; 0 0 0 0]
    lcs = [[1, 2, 5], [3, 4]]
    r = reactions(rn)
    subrn = [[r[1], r[4], r[3]], [r[2]]]
    lcd = [0, 0]
    testnetwork(rn, B, Z, Δ, lcs, 1, subrn, lcd)
end

# Some rational functions as rates.
let
    rn = @reaction_network begin
        (k₁, k₋₁), X₁ + X₂ <--> X₃ + 2X₄
        (k₂ / (1 + X₄ * X₅ + X₆ * X₇), k₋₂ / (1 + X₄ * X₅ + X₆ * X₇)), 3X₄ + X₅ <--> X₆ + X₇
        (k₃ / (1 + X₇ + X₈ + X₉ + X₁₀), k₋₃ / (1 + X₇ + X₈ + X₉ + X₁₀)), 5X₇ + X₈ <--> X₉ + X₁₀
    end
    Z = [1 0 0 0 0 0;
        1 0 0 0 0 0;
        0 1 0 0 0 0;
        0 2 3 0 0 0;
        0 0 1 0 0 0;
        0 0 0 1 0 0;
        0 0 0 1 5 0;
        0 0 0 0 1 0;
        0 0 0 0 0 1;
        0 0 0 0 0 1]
    B = [-1 1 0 0 0 0;
        1 -1 0 0 0 0;
        0 0 -1 1 0 0;
        0 0 1 -1 0 0;
        0 0 0 0 -1 1;
        0 0 0 0 1 -1]
    Δ = [-1 0 0 0 0 0; 0 -1 0 0 0 0; 0 0 -1 0 0 0; 0 0 0 -1 0 0; 0 0 0 0 -1 0; 0 0 0 0 0 -1]
    lcs = [[1, 2], [3, 4], [5, 6]]
    r = reactions(rn)
    subrn = [[r[1], r[2]], [r[3], r[4]], [r[5], r[6]]]
    lcd = [0, 0, 0]
    testnetwork(rn, B, Z, Δ, lcs, 0, subrn, lcd)
end

### Testing Reversibility ###

# Test function.
function testreversibility(rn, B, rev, weak_rev)
    @test isreversible(rn) == rev
    subrn = subnetworks(rn)
    @test isweaklyreversible(rn, subrn) == weak_rev
end

let
    rn = @reaction_network begin
        (k2, k1), A1 <--> A2 + A3
        k3, A2 + A3 --> A4
        k4, A4 --> A5
        (k6, k5), A5 <--> 2A6
        k7, 2A6 --> A4
        k8, A4 + A5 --> A7
    end
    rev = false
    weak_rev = false
    testreversibility(rn, reactioncomplexes(rn)[2], rev, weak_rev)
end

let
    rn = @reaction_network begin
        (k2, k1), A1 <--> A2 + A3
        k3, A2 + A3 --> A4
        k4, A4 --> A5
        (k6, k5), A5 <--> 2A6
        k7, A4 --> 2A6
        (k9, k8), A4 + A5 <--> A7
    end
    rev = false
    weak_rev = false
    testreversibility(rn, reactioncomplexes(rn)[2], rev, weak_rev)
end

let
    rn = @reaction_network begin
        k1, A --> B
        k2, A --> C
    end
    rev = false
    weak_rev = false
    testreversibility(rn, reactioncomplexes(rn)[2], rev, weak_rev)
end

let
    rn = @reaction_network begin
        k1, A --> B
        k2, A --> C
        k3, B + C --> 2A
    end
    rev = false
    weak_rev = false
    testreversibility(rn, reactioncomplexes(rn)[2], rev, weak_rev)
end

let
    rn = @reaction_network begin
        (k2, k1), A <--> 2B
        (k4, k3), A + C --> D
        k5, D --> B + E
        k6, B + E --> A + C
    end
    rev = false
    weak_rev = true
    testreversibility(rn, reactioncomplexes(rn)[2], rev, weak_rev)
end

let
    rn = @reaction_network begin
        (k2, k1), A + E <--> AE
        k3, AE --> B + E
    end
    rev = false
    weak_rev = false
    testreversibility(rn, reactioncomplexes(rn)[2], rev, weak_rev)
end

let
    rn = @reaction_network begin
        (k2, k1), A + E <--> AE
        (k4, k3), AE <--> B + E
    end
    rev = true
    weak_rev = true
    testreversibility(rn, reactioncomplexes(rn)[2], rev, weak_rev)
end

let
    rn = @reaction_network begin (k2, k1), A + B <--> 2A end
    rev = true
    weak_rev = true
    testreversibility(rn, reactioncomplexes(rn)[2], rev, weak_rev)
end

let
    rn = @reaction_network begin
        k1, A + B --> 3A
        k2, 3A --> 2A + C
        k3, 2A + C --> 2B
        k4, 2B --> A + B
    end
    rev = false
    weak_rev = true
    testreversibility(rn, reactioncomplexes(rn)[2], rev, weak_rev)
end

let
    rn = @reaction_network begin
        (k2, k1), A + E <--> AE
        (k4, k3), AE <--> B + E
        k5, B --> 0
        k6, 0 --> A
    end
    rev = false
    weak_rev = false
    testreversibility(rn, reactioncomplexes(rn)[2], rev, weak_rev)
end

# ------------------------------------------------------------------------- #

### More Tests ###

let
    myrn = [reaction_networks_standard; reaction_networks_hill; reaction_networks_real]
    for i in 1:length(myrn)
        local rcs, B = reactioncomplexes(myrn[i])
        @test B == Matrix(reactioncomplexes(myrn[i], sparse = true)[2])
        local Z = complexstoichmat(myrn[i])
        @test Z == Matrix(complexstoichmat(myrn[i], sparse = true))
        @test Z * B == netstoichmat(myrn[i]) == Matrix(netstoichmat(myrn[i], sparse = true))
    end
end

# Test defaults.
let
    rn = @reaction_network begin
        α, S + I --> 2I
        β, I --> R
    end
    p = [0.1 / 1000, 0.01]
    tspan = (0.0, 250.0)
    u0 = [999.0, 1.0, 0.0]
    op = ODEProblem(rn, species(rn) .=> u0, tspan, parameters(rn) .=> p)
    sol = solve(op, Tsit5())  # old style
    setdefaults!(rn, [:S => 999.0, :I => 1.0, :R => 0.0, :α => 1e-4, :β => 0.01])
    op = ODEProblem(rn, [], tspan, [])
    sol2 = solve(op, Tsit5())
    @test norm(sol.u - sol2.u) ≈ 0
    @test all(p -> p[1] isa Symbolics.Symbolic, collect(ModelingToolkit.defaults(rn)))

    rn = @reaction_network begin
        α, S + I --> 2I
        β, I --> R
    end
    @parameters α β
    @species S(t) I(t) R(t)
    setdefaults!(rn, [S => 999.0, I => 1.0, R => 0.0, α => 1e-4, β => 0.01])
    op = ODEProblem(rn, [], tspan, [])
    sol2 = solve(op, Tsit5())
    @test norm(sol.u - sol2.u) ≈ 0

    # Rest unpacking variables.
    function unpacktest(rn)
        Catalyst.@unpacksys rn
        u₀ = [S1 => 999.0, I1 => 1.0, R1 => 0.0]
        p = [α1 => 1e-4, β1 => 0.01]
        op = ODEProblem(rn, u₀, (0.0, 250.0), p)
        solve(op, Tsit5())
    end
    rn = @reaction_network begin
        α1, S1 + I1 --> 2I1
        β1, I1 --> R1
    end
    sol3 = unpacktest(rn)
    @test norm(sol.u - sol3.u) ≈ 0

    # Test symmap_to_varmap.
    sir = @reaction_network sir begin
        β, S + I --> 2I
        ν, I --> R
    end
    subsys = @reaction_network subsys begin k, A --> B end
    @named sys = compose(sir, [subsys])
    symmap = [:S => 1.0, :I => 1.0, :R => 1.0, :subsys₊A => 1.0, :subsys₊B => 1.0]
    u0map = symmap_to_varmap(sys, symmap)
    pmap = symmap_to_varmap(sys, [:β => 1.0, :ν => 1.0, :subsys₊k => 1.0])
    @test isequal(u0map[4][1], subsys.A)
    @test isequal(u0map[1][1], @nonamespace sir.S)

    u0map = symmap_to_varmap(sir, [:S => 999.0, :I => 1.0, :R => 0.0])
    pmap = symmap_to_varmap(sir, [:β => 1e-4, :ν => 0.01])
    op = ODEProblem(sir, u0map, tspan, pmap)
    sol4 = solve(op, Tsit5())
    @test norm(sol.u - sol4.u) ≈ 0

    u0map = [:S => 999.0, :I => 1.0, :R => 0.0]
    pmap = (:β => 1e-4, :ν => 0.01)
    op = ODEProblem(sir, u0map, tspan, pmap)
    sol5 = solve(op, Tsit5())
    @test norm(sol.u - sol5.u) ≈ 0
end

# Test conservation law elimination.
let
    rn = @reaction_network begin
        (k1, k2), A + B <--> C
        (m1, m2), D <--> E
        b12, F1 --> F2
        b23, F2 --> F3
        b31, F3 --> F1
    end
    osys = convert(ODESystem, rn; remove_conserved = true)
    @unpack A, B, C, D, E, F1, F2, F3, k1, k2, m1, m2, b12, b23, b31 = osys
    u0 = [A => 10.0, B => 10.0, C => 0.0, D => 10.0, E => 0.0, F1 => 8.0, F2 => 0.0,
        F3 => 0.0]
    p = [k1 => 1.0, k2 => 0.1, m1 => 1.0, m2 => 2.0, b12 => 1.0, b23 => 2.0, b31 => 0.1]
    tspan = (0.0, 20.0)
    oprob = ODEProblem(osys, u0, tspan, p)
    sol = solve(oprob, Tsit5(); abstol = 1e-10, reltol = 1e-10)
    oprob2 = ODEProblem(rn, u0, tspan, p)
    sol2 = solve(oprob2, Tsit5(); abstol = 1e-10, reltol = 1e-10)
    oprob3 = ODEProblem(rn, u0, tspan, p; remove_conserved = true)
    sol3 = solve(oprob3, Tsit5(); abstol = 1e-10, reltol = 1e-10)

    tv = range(tspan[1], tspan[2], length = 101)
    for s in species(rn)
        @test isapprox(sol(tv, idxs = s), sol2(tv, idxs = s))
        @test isapprox(sol2(tv, idxs = s), sol2(tv, idxs = s))
    end

    nsys = convert(NonlinearSystem, rn; remove_conserved = true)
    nprob = NonlinearProblem{true}(nsys, u0, p)
    nsol = solve(nprob, NewtonRaphson(); abstol = 1e-10)
    nprob2 = ODEProblem(rn, u0, (0.0, 100.0 * tspan[2]), p)
    nsol2 = solve(nprob2, Tsit5(); abstol = 1e-10, reltol = 1e-10)
    nprob3 = NonlinearProblem(rn, u0, p; remove_conserved = true)
    nsol3 = solve(nprob3, NewtonRaphson(); abstol = 1e-10)
    for s in species(rn)
        @test isapprox(nsol[s], nsol2(tspan[2], idxs = s))
        @test isapprox(nsol2(tspan[2], idxs = s), nsol3[s])
    end

    u0 = [A => 100.0, B => 20.0, C => 5.0, D => 10.0, E => 3.0, F1 => 8.0, F2 => 2.0,
        F3 => 20.0]
    ssys = convert(SDESystem, rn; remove_conserved = true)
    sprob = SDEProblem(ssys, u0, tspan, p)
    sprob2 = SDEProblem(rn, u0, tspan, p)
    sprob3 = SDEProblem(rn, u0, tspan, p; remove_conserved = true)
    ists = ModelingToolkit.get_unknowns(ssys)
    sts = ModelingToolkit.get_unknowns(rn)
    istsidxs = findall(in(ists), sts)
    u1 = copy(sprob.u0)
    u2 = sprob2.u0
    u3 = copy(sprob3.u0)
    du1 = similar(u1)
    du2 = similar(u2)
    du3 = similar(u3)
    g1 = zeros(length(u1), numreactions(rn))
    g2 = zeros(length(u2), numreactions(rn))
    g3 = zeros(length(u3), numreactions(rn))
    sprob.f(du1, u1, sprob.p, 1.0)
    sprob2.f(du2, u2, sprob2.p, 1.0)
    sprob3.f(du3, u3, sprob3.p, 1.0)
    @test isapprox(du1, du2[istsidxs])
    @test isapprox(du2[istsidxs], du3)
    sprob.g(g1, u1, sprob.p, 1.0)
    sprob2.g(g2, u2, sprob2.p, 1.0)
    sprob3.g(g3, u3, sprob3.p, 1.0)
    @test isapprox(g1, g2[istsidxs, :])
    @test isapprox(g2[istsidxs, :], g3)
end


# Non-integer stoichiometry.
let
    function test_stoich(T, rn)
        @test eltype(substoichmat(rn)) == T
        @test eltype(prodstoichmat(rn)) == T
        @test eltype(netstoichmat(rn)) == T
        @test eltype(substoichmat(rn; sparse = true)) == T
        @test eltype(prodstoichmat(rn; sparse = true)) == T
        @test eltype(netstoichmat(rn; sparse = true)) == T
        nothing
    end

    rn = @reaction_network ABtoC begin
        (k₊,k₋), 3.4*A + 2B <--> 2.5*C
    end
    test_stoich(Float64, rn)

    rn2 = @reaction_network ABtoC begin
        (k₊,k₋), 3*A + 2B <--> 2*C
    end
    test_stoich(Int, rn2)
end

### Miscelenesous Tests ###

# Tests various additional API functions.
let
    rn = @reaction_network begin
        (p,d), 0 <--> X
        (kB,kD), 2X <--> X
    end
    @unpack p, d, kB, kD = rn
    isequal(reactionparamsmap(rn), Dict([p => 1, d => 2, kB => 3, kD => 4]))
    issetequal(reactionrates(rn), [p, d, kB, kD])
    isequal(symmap_to_varmap(rn, [:p => 1.0, :kB => 3.0]), [p => 1.0, kB => 3.0])
end

### Test Polynomial Transformation Functionality ###

# Tests normal network.
let
    rn = @reaction_network begin
        (p,d), 0 <--> X
        (kB,kD), 2X <--> X2
    end
    ns = convert(NonlinearSystem, rn)
    neweqs = getfield.(equations(ns),:rhs)
    poly = Catalyst.to_multivariate_poly(neweqs)
    @test length(poly) == 2
end

# Tests network with a fraction.
let
    rn = @reaction_network begin
        (p/X,d), 0 <--> X
    end
    ns = convert(NonlinearSystem, rn)
    neweqs = getfield.(equations(ns),:rhs)
    poly = Catalyst.to_multivariate_poly(neweqs)
    @test length(poly) == 1
end

# Test empty network.
let
    rn = @reaction_network
    ns = convert(NonlinearSystem, rn)
    neweqs = getfield.(equations(ns),:rhs)
    @test_throws MethodError Catalyst.to_multivariate_poly(neweqs)
end
