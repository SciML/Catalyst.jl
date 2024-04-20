### Prepares Tests ###

# Fetch packages.
using Catalyst, LinearAlgebra, Test

# Fetch test networks.
include("../test_networks.jl")

### Basic Tests ###

# Tests basic functionality on system with known conservation laws.
let
    rn = @reaction_network begin
        (1, 2), A + B <--> C
        (3, 2), D <--> E
        (0.1, 0.2), E <--> F
        6, F --> G
        7, H --> G
        5, 0 --> K
        3, A + B --> Z + C
    end

    S = netstoichmat(rn)
    C = conservationlaws(S)
    @test size(C, 1) == 3
    b = [0, 0, 0, 1, 1, 1, 1, 1, 0, 0]
    @test any(b == C[i, :] for i in 1:size(C, 1))

    # For the A + B <--> C subsystem one of these must occur
    # as a conservation law.
    D = [1 -1 0 0 0 0 0 0 0 0;
         -1 1 0 0 0 0 0 0 0 0
         0 1 1 0 0 0 0 0 0 0]
    @test any(D[j, :] == C[i, :] for i in 1:size(C, 1), j in 1:size(D, 1))

    C = conservationlaws(rn)
    @test size(C, 1) == 3
    @test Catalyst.get_networkproperties(rn).nullity == 3
    @test any(b == C[i, :] for i in 1:size(C, 1))
    @test any(D[j, :] == C[i, :] for i in 1:size(C, 1), j in 1:size(D, 1))
end

# Tests conservation law computation on large number of networks where we know which have conservation laws.
let
    Cs_standard = map(conservationlaws, reaction_networks_standard)
    @test all(size(C, 1) == 0 for C in Cs_standard)

    Cs_hill = map(conservationlaws, reaction_networks_hill)
    @test all(size(C, 1) == 0 for C in Cs_hill)

    function consequiv(A, B)
        rank([A; B]) == rank(A) == rank(B)
    end
    Cs_constraint = map(conservationlaws, reaction_networks_constraint)
    @test all(consequiv.(Matrix{Int}.(Cs_constraint), reaction_network_constraints))
end

# Tests additional conservation law-related functions.
let
    rn = @reaction_network begin
        (k1, k2), X1 <--> X2
        (k3, k4), X3 <--> X4
    end
    cons_laws = conservationlaws(rn)
    cons_eqs = conservedequations(rn)
    cons_laws_constants = conservationlaw_constants(rn)
    conserved_quantity = conservedquantities(cons_laws[1, :], unknowns(rn)[1])

    @test sum(cons_laws) == 4
    @test size(cons_laws) == (2, 4)
    @test length(cons_eqs) == 2
    @test length(conserved_quantity) == 4
    @test length(cons_laws_constants) == 2
    @test count(isequal.(conserved_quantity, Num(0))) == 2
end

# Conservation law simulations for vectorised species.
let 
    # Prepares the model.
    t = default_t()
    @species X(t)[1:2]
    @parameters k[1:2]
    rxs = [
        Reaction(k[1], [X[1]], [X[2]]),
        Reaction(k[2], [X[2]], [X[1]])
    ]
    @named rs = ReactionSystem(rxs, t)
    rs = complete(rs)

    # Checks that simulation reaches known equilibrium
    u0 = [:X => [3.0, 9.0]]
    ps = [:k => [1.0, 2.0]]
    oprob = ODEProblem(rs, u0, (0.0, 1000.0), ps; remove_conserved = true)
    sol = solve(oprob, Vern7())
    @test sol[X[1]] ≈ 8.0
    @test sol[X[2]] ≈ 4.0
end