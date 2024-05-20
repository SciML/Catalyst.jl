### Prepares Tests ###

# Fetch packages.
using Catalyst, LinearAlgebra, NonlinearSolve, OrdinaryDiffEq

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

# Test conservation law elimination.
let
    rn = @reaction_network begin
        (k1, k2), A + B <--> C
        (m1, m2), D <--> E
        b12, F1 --> F2
        b23, F2 --> F3
        b31, F3 --> F1
    end
    osys = complete(convert(ODESystem, rn; remove_conserved = true))
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

    nsys = complete(convert(NonlinearSystem, rn; remove_conserved = true))
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
    ssys = complete(convert(SDESystem, rn; remove_conserved = true))
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

### ConservedQuantity Metadata Tests ###

# Checks that `conservationquantity` metadata is added correctly to parameters.
# Checks that the `isconservationquantity` getter function works correctly.
let
    # Creates ODESystem with conserved quantities.
    rs = @reaction_network begin
        (k1,k2), X1 <--> X2
        (k1,k2), Y1 <--> Y2
    end
    osys = convert(ODESystem, rs)

    # Checks that the correct parameters have the `conservationquantity` metadata.
    @test Catalyst.isconservationquantity(osys.Γ[1])
    @test Catalyst.isconservationquantity(osys.Γ[2])
    @test !Catalyst.isconservationquantity(osys.k1)
    @test !Catalyst.isconservationquantity(osys.k2)
end
