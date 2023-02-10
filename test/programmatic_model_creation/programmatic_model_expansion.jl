#! format: off

### Fetch required packages and reaction networks ###
using Catalyst, Test
using ModelingToolkit: get_ps, get_states, get_eqs, get_systems, get_iv
include("../test_networks.jl")

using StableRNGs
rng = StableRNG(12345)

function unpacksys(sys)
    get_eqs(sys), get_iv(sys), get_ps(sys), nameof(sys), get_systems(sys)
end

### Tests construction of empty reaction networks ###
empty_network_1 = @reaction_network
eqs, iv, ps, name, systems = unpacksys(empty_network_1)
@test length(eqs) == 0
@test nameof(iv) == :t
@test length(get_states(empty_network_1)) == 0
@test length(ps) == 0

empty_network_2 = @reaction_network
@variables p1 p2 p3 p4 p5
addparam!(empty_network_2, p1)
addparam!(empty_network_2, p2)
addparam!(empty_network_2, p3)
addparam!(empty_network_2, p4)
addparam!(empty_network_2, p5)
eqs, iv, ps, name, systems = unpacksys(empty_network_2)
@test length(eqs) == 0
@test nameof(iv) == :t
@test length(get_states(empty_network_2)) == 0
@test length(ps) == 5
@test all(getproperty.(ps, :name) .== [:p1, :p2, :p3, :p4, :p5])

### Tests accessing parameters and species added with network API ###
empty_network_3 = @reaction_network
@variables x p
addspecies!(empty_network_3, x)
addparam!(empty_network_3, p)
@test isequal(empty_network_3.x, states(empty_network_3, x))
@test isequal(empty_network_3.p, parameters(empty_network_3, p))

### Tests creating a network and adding reactions ###
unfinished_network = @reaction_network begin
    @parameters k0 k1 k2 k3 k4
    (k1, k2), X1 ↔ X2
    (k3, k4), X3 ↔ X4
end
@add_reactions unfinished_network begin
    @parameters k0 k3 k4 k5 k6
    (k3, k4), X3 ↔ X4
    (k5, k6), X5 ↔ X6
end
@add_reactions unfinished_network begin
    @parameters k0 k5 k6 k7 k8
    (k5, k6), X5 ↔ X6
    (k7, k8), X7 ↔ X8
end
@test length(get_states(unfinished_network)) == 8
@test length(get_ps(unfinished_network)) == 9

### Compares test network to identical network constructed via @add_reactions ###
identical_networks = Vector{Pair}()

step_by_step_network_1 = @reaction_network rns5 begin
    @parameters p k1 k2 k3 k4
    p, ∅ → X1
    (k1, k2), X1 ⟷ X2
    (k3, k4), X2 ⟷ X3
end
@add_reactions step_by_step_network_1 begin
    @parameters k5 k6
    (k5, k6), X3 ⟷ X4
end
@add_reactions step_by_step_network_1 begin
    @parameters d
    d, X4 → ∅
end
push!(identical_networks, reaction_networks_standard[5] => step_by_step_network_1)

step_by_step_network_2 = @reaction_network rns7 begin (p1, p2, p3), ∅ → (X1, X2, X3) end
@add_reactions step_by_step_network_2 begin
    @parameters k1 k2 k3 v1 K1 d1 d2 d3 d4 d5
    (k1, k2), X1 + 2X2 ⟷ X4
    (mm(X3, v1, K1), k3), X4 ⟷ X5
    (d1, d2, d3, d4, d5), (X1, X2, X3, X4, X5) → ∅
end
push!(identical_networks, reaction_networks_standard[7] => step_by_step_network_2)

step_by_step_network_3 = @reaction_network rns10 begin
    @parameters p k1 k2 k3 k4
    p, ∅ ⟶ X1
    (k1, k2), X1 → X2
end
@add_reactions step_by_step_network_3 begin
    @parameters p k1 k2 k3 k4 k5 k6
    (k3, k4), X2 → X3
    (k5, k6), X3 → X4
end
@add_reactions step_by_step_network_3 begin
    @parameters p k1 k2 k3 k4 k5 k6 k7 k8 d
    (k7, k8), X4 → X5
    d, X5 ⟶ ∅
end
push!(identical_networks, reaction_networks_standard[10] => step_by_step_network_3)

step_by_step_network_4 = @reaction_network rnh7 begin
    @parameters v K n k1 k2 k3 d
    v / 10 + hill(X1, v, K, n),
    ∅ → X1 + X2
end
@add_reactions step_by_step_network_4 begin
    @parameters k1 k2
    k1, X1 + X2 → X3
    k2, X3 → X1 + X2
end
@add_reactions step_by_step_network_4 begin
    @parameters k2 k3 d
    k3, X3 → X1
    d, (X1, X2, X3) → ∅
end
push!(identical_networks, reaction_networks_hill[7] => step_by_step_network_4)

step_by_step_network_5 = @reaction_network rnc1 begin
    @parameters k1 k2 k3 k4 k5 k6
    (k1, k2), X1 ↔ X2
    (k3, k4), X2 ↔ X3
end
@add_reactions step_by_step_network_5 begin
    @parameters k5
    k5, X3 → X1
end
@add_reactions step_by_step_network_5 begin
    @parameters k6
    k6, X1 → X3
end
push!(identical_networks, reaction_networks_constraint[1] => step_by_step_network_5)

step_by_step_network_6 = @reaction_network rnc5 begin
    @parameters k1 k2 k3 k4 k5 k6
    (k1, k2), X1 ↔ 2X2
end
@add_reactions step_by_step_network_6 begin
    @parameters k3 k4
    (k3, k4), 2X2 ↔ 3X3
end
@add_reactions step_by_step_network_6 begin
    @parameters k5 k6
    (k5, k6), 3X3 ↔ 4X4
end
push!(identical_networks, reaction_networks_constraint[5] => step_by_step_network_6)

step_by_step_network_7 = @reaction_network rnr3 begin
    @parameters k1 k2p k2pp k3p k3pp A J3 k4 m J4
    k2p, Y → 0
end
@add_reactions step_by_step_network_7 begin
@parameters k3p k3pp A J3 k2pp
    (k3p + k3pp * A) / (J3 + Po), Po → P
    k2pp * P, Y → 0
end
@add_reactions step_by_step_network_7 begin
    @parameters k1
    k1, 0 → Y
end
@add_reactions step_by_step_network_7 begin
    @parameters k4 m J4
    (k4 * m) / (J4 + P), Y + P → Y + Po
end
push!(identical_networks, reaction_networks_real[3] => step_by_step_network_7)

@parameters k1
step_by_step_network_8 = @reaction_network rnw7
addparam!(step_by_step_network_8, k1)
@add_reactions step_by_step_network_8 begin
    @parameters k1
    k1, X1 → X2
    0, X2 → X3
end
@add_reactions step_by_step_network_8 begin
    @parameters k2 k3
    k2, X3 → X4
    k3, X4 → X5
end
push!(identical_networks, reaction_networks_weird[7] => step_by_step_network_8)

step_by_step_network_9 = @reaction_network rnw10
@add_reactions step_by_step_network_9 begin
    @parameters d
    d, 5X1 → 4X1
end
push!(identical_networks, reaction_networks_weird[10] => step_by_step_network_9)

function permute_ps(pvals, rn1, rn2)
    ps1 = parameters(rn1)
    ps2 = parameters(rn2)
    pvals2 = similar(pvals)
    for (i, p) in enumerate(ps2)
        pidx = findfirst(isequal(p), ps1)
        pvals2[i] = pvals[pidx]
    end
    pvals2
end

for networks in identical_networks
    f1 = ODEFunction(convert(ODESystem, networks[1]), jac = true)
    f2 = ODEFunction(convert(ODESystem, networks[2]), jac = true)
    g1 = SDEFunction(convert(SDESystem, networks[1]))
    g2 = SDEFunction(convert(SDESystem, networks[2]))
    @test networks[1] == networks[2]
    for factor in [1e-2, 1e-1, 1e0, 1e1]
        u0 = factor * rand(rng, length(get_states(networks[1])))
        pp = factor * rand(rng, length(get_ps(networks[1])))

        # needed because this code assumes an ordering of the parameters and species...
        pp2 = permute_ps(pp, networks[1], networks[2])
        t = rand(rng)
        @test all(abs.(f1(u0, pp, t) .- f2(u0, pp2, t)) .< 1000 * eps())
        @test all(abs.(f1.jac(u0, pp, t) .- f2.jac(u0, pp2, t)) .< 1000 * eps())
        @test all(abs.(g1(u0, pp, t) .- g2(u0, pp2, t)) .< 1000 * eps())
    end
end
