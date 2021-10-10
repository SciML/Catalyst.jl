### Fetch required packages and reaction networks ###
using DiffEqBase, Catalyst, Random, Test
using ModelingToolkit: operation, Sym, istree, get_states, get_ps, get_eqs, get_systems, get_iv

using StableRNGs
rng = StableRNG(12345)

include("test_networks.jl")

function unpacksys(sys)
    get_eqs(sys),get_iv(sys),get_states(sys),get_ps(sys),nameof(sys),get_systems(sys)
end

### Debug functions ###
opname(x) = istree(x) ? nameof(operation(x)) : nameof(x)
alleq(xs,ys) = all(isequal(x,y) for (x, y) in zip(xs, ys))

# Gets all the reactants in a set of equations.
function all_reactants(eqs)
    all_reactants = []
    for eq in eqs
        append!(all_reactants,opname.(eq.substrates))
        append!(all_reactants,opname.(eq.products))
    end
    return Set{Symbol}(unique(all_reactants))
end

# Gets all parameters (where every reaction rate is constant)
function all_parameters(eqs)
    return Set(unique(map(eq -> opname(eq.rate),eqs)))
end

### Test basic properties of networks ###
function basic_test(rn, N, states_syms, p_syms)
    eqs,iv,states,ps,name,systems = unpacksys(rn)
    @test length(eqs) == N
    @test opname(iv) == :t
    @test length(states) == length(states_syms)
    @test issetequal(map(opname, states), states_syms)
    @test all_reactants(eqs) == Set(states_syms)
    @test length(ps) == length(p_syms)
    @test issetequal(map(opname, ps), p_syms)
end
basic_test(reaction_networks_standard[1], 10, [:X1,:X2,:X3], [:p1,:p2,:p3,:k1,:k2,:k3,:k4,:d1,:d2,:d3])
@test all_parameters(get_eqs(reaction_networks_standard[1])) == Set([:p1,:p2,:p3,:k1,:k2,:k3,:k4,:d1,:d2,:d3])
basic_test(reaction_networks_standard[2], 3, [:X1,:X2], [:v1,:K1,:v2,:K2,:d])   
basic_test(reaction_networks_standard[3], 10, [:X1,:X2,:X3,:X4], [:v1,:K1,:v2,:K2,:k1,:k2,:k3,:k4,:d])
basic_test(reaction_networks_standard[4], 8, [:X1,:X2,:X3,:X4], [:v1,:K1,:v2,:K2,:v3,:K3,:v4,:K4,:d1,:d2,:d3,:d4])
basic_test(reaction_networks_standard[5], 8, [:X1,:X2,:X3,:X4], [:p,:k1,:k2,:k3,:k4,:k5,:k6,:d])
@test all_parameters(get_eqs(reaction_networks_standard[5])) == Set([:p,:k1,:k2,:k3,:k4,:k5,:k6,:d])
basic_test(reaction_networks_hill[1], 4, [:X1,:X2], [:v1,:v2,:K1,:K2,:n1,:n2,:d1,:d2])
basic_test(reaction_networks_constraint[1], 6, [:X1,:X2,:X3], [:k1,:k2,:k3,:k4,:k5,:k6])
basic_test(reaction_networks_real[1], 4, [:X,:Y], [:A,:B])
basic_test(reaction_networks_weird[1], 2, [:X], [:p,:d])
basic_test(reaction_networks_weird[2], 4, [:X,:Y,:Z], [:k1,:k2,:k3,:k4])


### Tries making various systems ###
identical_networks_1 = Vector{Pair}()

different_arrow_1 = @reaction_network begin
    (p1,p2,p3), ∅ > (X1,X2,X3)
    (k1,k2), X2 ↔ X1 + 2X3
    (k3,k4), X1 ⟷ X3
    (d1,d2,d3), (X1,X2,X3) → ∅
end p1 p2 p3 k1 k2 k3 k4 d1 d2 d3
push!(identical_networks_1, reaction_networks_standard[1] => different_arrow_1)

different_arrow_2 = @reaction_network begin
    mmr(X2,v1,K1), ∅ → X1
    mm(X1,v2,K2), ∅ ↣ X2
    d, X1+X2 ↦ ∅
end v1 K1 v2 K2 d
push!(identical_networks_1, reaction_networks_standard[2] => different_arrow_2)

different_arrow_3 = @reaction_network begin
    mm(X2,v1,K1), ∅ ⇾ X1
    mm(X3,v2,K2), ∅ ⟶ X2
    (k1,k2), X1 ⇄ X3
    (k3,k4), X3 + X2 ⇆ X4 +X1
    d, (X1,X2,X3,X4) ⟼ ∅
end v1 K1 v2 K2 k1 k2 k3 k4 d
push!(identical_networks_1, reaction_networks_standard[3] => different_arrow_3)

different_arrow_4 = @reaction_network begin
    mmr(X4,v1,K1), ∅ ⥟ X1
    mmr(X1,v2,K2), ∅ ⥟ X2
    mmr(X2,v3,K3), ∅ ⇀ X3
    mmr(X3,v4,K4), ∅ ⇁ X4
    (d1,d2,d3,d4), (X1,X2,X3,X4) --> ∅
end v1 K1 v2 K2 v3 K3 v4 K4 d1 d2 d3 d4
push!(identical_networks_1, reaction_networks_standard[4] => different_arrow_4)

# Yes the name is different, I wanted one with several single direction arrows.
different_arrow_8 = @reaction_network begin
    p, 2X1 < ∅
    k1, X2 ← X1
    (k2, k3), X3 ⟻ X2
    d, ∅ ↼ X3
end p k1 k2 k3 d
push!(identical_networks_1, reaction_networks_standard[8] => different_arrow_8)

for networks in identical_networks_1
    f1 = ODEFunction(convert(ODESystem,networks[1]),jac=true)
    f2 = ODEFunction(convert(ODESystem,networks[2]),jac=true)
    g1 = SDEFunction(convert(SDESystem,networks[1]))
    g2 = SDEFunction(convert(SDESystem,networks[2]))
    for factor in [1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3]
        u0 = factor*rand(rng,length(get_states(networks[1])))
        p = factor*rand(rng,length(get_ps(networks[1])))
        t = rand(rng)
        @test all(abs.(f1(u0,p,t) .≈ f2(u0,p,t)))
        @test all(abs.(f1.jac(u0,p,t) .≈ f2.jac(u0,p,t)))
        @test all(abs.(g1(u0,p,t) .≈ g2(u0,p,t)))
    end
end


### Tests that networks expressed in different ways are identical ###
identical_networks_2 = Vector{Pair}()

# Different parameter and variable names.
differently_written_5 = @reaction_network begin
    q, ∅ → Y1
    (l1,l2), Y1 ⟷ Y2
    (l3,l4), Y2 ⟷ Y3
    (l5,l6), Y3 ⟷ Y4
    c, Y4 → ∅
end q l1 l2 l3 l4 l5 l6 c
push!(identical_networks_2, reaction_networks_standard[5] => differently_written_5)

# Unfold reactions.
differently_written_6 = @reaction_network begin
    p1, ∅ → X1
    p2, ∅ → X2
    k1, 2X1 → X3
    k2, X3 → 2X1
    k3, X2 + X3 → 4X4
    k4, 4X4 → X2 + X3
    k5, X4 + X1 → 2X3
    k6, 2X3 → X4 + X1
    d, X1 → ∅
    d, X2 → ∅
    d, X3 → ∅
    d, X4 → ∅
end p1 p2 k1 k2 k3 k4 k5 k6 d
push!(identical_networks_2, reaction_networks_standard[6] => differently_written_6)

# Ignore mass action.
differently_written_7 = @reaction_network begin
    (p1,p2,p3), ∅ ⇒ (X1,X2,X3)
    (k1*X1*X2^2/2,k2*X4), X1 + 2X2 ⟺ X4
    (mm(X3,v1,K1)*X4,k3*X5), X4 ⇔ X5
    (d1*X1,d2*X2,d3*X3,d4*X4,d5*X5), ∅ ⟽ (X1,X2,X3,X4,X5)
end p1 p2 p3 k1 k2 k3 v1 K1 d1 d2 d3 d4 d5
push!(identical_networks_2, reaction_networks_standard[7] => differently_written_7)

for networks in identical_networks_2
    f1 = ODEFunction(convert(ODESystem,networks[1]),jac=true)
    f2 = ODEFunction(convert(ODESystem,networks[2]),jac=true)
    g1 = SDEFunction(convert(SDESystem,networks[1]))
    g2 = SDEFunction(convert(SDESystem,networks[2]))
    for factor in [1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3]
        u0 = factor*rand(rng,length(get_states(networks[1])))
        p = factor*rand(rng,length(get_ps(networks[1])))
        t = rand(rng)
        @test all(f1(u0,p,t) .≈ f2(u0,p,t))
        @test all(f1.jac(u0,p,t) .≈ f2.jac(u0,p,t))
        @test all(g1(u0,p,t) .≈ g2(u0,p,t))
    end
end


### Test networks without parameters ###
identical_networks_3 = Vector{Pair}()
parameter_sets = []

# Different parameter and variable names.
no_parameters_9 = @reaction_network begin
    (1.5,1,2), ∅ ⟶ (X1,X2,X3)
    (0.01,2.3,1001), (X1,X2,X3) ⟶ ∅
    (π,42), X1 + X2 ⟷ X3
    (19.9,999.99), X3 ⟷ X4
    (sqrt(3.7),exp(1.9)), X4 ⟷ X1 + X2
end
push!(identical_networks_3, reaction_networks_standard[9] => no_parameters_9)
push!(parameter_sets, [1.5,1,2,0.01,2.3,1001,π,42,19.9,999.99,sqrt(3.7),exp(1.9)])

no_parameters_10 = @reaction_network begin
    0.01, ∅ ⟶ X1
    (3.1,3.2), X1 → X2
    (0.,2.1), X2 → X3
    (901.,63.5), X3 → X4
    (7,8), X4 → X5
    1., X5  ⟶ ∅
end
push!(identical_networks_3, reaction_networks_standard[10] => no_parameters_10)
push!(parameter_sets, [0.01,3.1,3.2,0.,2.1,901.,63.5,7,8,1.])

for (i,networks) in enumerate(identical_networks_3)
    f1 = ODEFunction(convert(ODESystem,networks[1]),jac=true)
    f2 = ODEFunction(convert(ODESystem,networks[2]),jac=true)
    g1 = SDEFunction(convert(SDESystem,networks[1]))
    g2 = SDEFunction(convert(SDESystem,networks[2]))
    for factor in [1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3]
        u0 = factor*rand(rng,length(get_states(networks[1])))
        t = rand(rng)
        @test f1(u0,parameter_sets[i],t) ≈ f2(u0,[],t)
        @test f1.jac(u0,parameter_sets[i],t) ≈ f2.jac(u0,[],t)
        @test g1(u0,parameter_sets[i],t) ≈ g2(u0,[],t)
    end
end


### Tests that Reaction System created manually and through macro are identical ###
identical_networks_4 = Vector{Pair}()
@parameters v1 K1 v2 K2 k1 k2 k3 k4 k5 p d t
@variables X1(t) X2(t) X3(t) X4(t) X5(t)

rxs_1 = [Reaction(p, nothing, [X1], nothing, [2]),
       Reaction(k1, [X1], [X2], [1], [1]),
       Reaction(k2, [X2], [X3], [1], [1]),
       Reaction(k3, [X2], [X3], [1], [1]),
       Reaction(d, [X3], nothing, [1], nothing)]
@named rs_1 = ReactionSystem(rxs_1 , t, [X1,X2,X3], [p,k1,k2,k3,d])
push!(identical_networks_4, reaction_networks_standard[8] => rs_1)

rxs_2 = [Reaction(k1, [X1], [X2], [1], [1]),
         Reaction(k2*X5, [X2], [X1], [1], [1]),
         Reaction(k3*X5, [X3], [X4], [1], [1]),
         Reaction(k4, [X4], [X3], [1], [1]),
         Reaction(p+k5*X2*X3, nothing, [X5], nothing, [1]),
         Reaction(d, [X5], nothing, [1], nothing)]
@named rs_2 = ReactionSystem(rxs_2, t, [X1,X2,X3,X4,X5], [k1,k2,k3,k4,p,k5,d])
push!(identical_networks_4, reaction_networks_constraint[3] => rs_2)

rxs_3 = [Reaction(k1, [X1], [X2], [1], [1]),
         Reaction(0, [X2], [X3], [1], [1]),
         Reaction(k2, [X3], [X4], [1], [1]),
         Reaction(k3, [X4], [X5], [1], [1])]
@named rs_3 = ReactionSystem(rxs_3, t, [X1,X2,X3,X4,X5], [k1,k2,k3])
push!(identical_networks_4, reaction_networks_weird[7] => rs_3)

for networks in identical_networks_4
    @test isequal(get_iv(networks[1]), get_iv(networks[2]))
    @test alleq(get_states(networks[1]), get_states(networks[2]))
    @test alleq(get_ps(networks[1]), get_ps(networks[2]))
    @test ModelingToolkit.get_systems(networks[1]) == ModelingToolkit.get_systems(networks[2])
    @test length(get_eqs(networks[1])) == length(get_eqs(networks[2]))
    for (e1,e2) in zip(get_eqs(networks[1]),get_eqs(networks[2]))
        @test isequal(e1.rate,e2.rate)
        @test isequal(e1.substrates,e2.substrates)
        @test isequal(e1.products,e2.products)
        @test isequal(e1.substoich,e2.substoich)
        @test isequal(e1.prodstoich,e2.prodstoich)
        @test isequal(e1.netstoich,e2.netstoich)
        @test isequal(e1.only_use_rate,e2.only_use_rate)
    end
end


### Tests that time is handled properly ###
time_network = @reaction_network begin
    (t,k2), X1 ↔ X2
    (k3,t), X2 ↔ X3
    (t,k6), X3 ↔ X1
end k2 k3 k6

f1 = ODEFunction(convert(ODESystem,reaction_networks_constraint[1]),jac=true)
f2 = ODEFunction(convert(ODESystem,time_network),jac=true)
g1 = SDEFunction(convert(SDESystem,reaction_networks_constraint[1]))
g2 = SDEFunction(convert(SDESystem,time_network))
for factor in [1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3]
    u0 = factor*rand(rng,length(get_states(time_network)))
    κ2 = factor*rand(rng); κ3 = factor*rand(rng); κ6 = factor*rand(rng);
    τ = rand(rng)
    p1 = [τ, κ2, κ3, τ, τ, κ6]; p2 = [κ2, κ3, κ6];
    @test all(f1(u0,p1,τ) .≈ f2(u0,p2,τ))
    @test all(f1.jac(u0,p1,τ) .≈ f2.jac(u0,p2,τ))
    @test all(g1(u0,p1,τ) .≈ g2(u0,p2,τ))
end


### Test various names as varriables ###
test_network = @reaction_network begin
    (a,A),  n ⟷ N
    (b,B),  o ⟷ O
    (c,C),  p ⟷ P
    (d,D),  q ⟷ Q
    (e,E),  r ⟷ R
    (f,F),  s ⟷ S
    (g,G),  u ⟷ U
    (h,H),  v ⟷ V
    (j,J),  w ⟷ W
    (k,K),  x ⟷ X
    (l,L),  y ⟷ Y
    (m,M),  z ⟷ Z
end a A b B c C d D e E f F g G h H j J k K l L m M

test_network = @reaction_network begin
    (1.,1.),  i ⟷ T
end

test_network = @reaction_network begin
    (å,Å), ü ⟷ Ü
    (ä,Ä), ñ ⟷ Ñ
    (ö,Ö), æ ⟷ Æ
end å Å ä Ä ö Ö

test_network = @reaction_network begin
    (α,Α), ν ⟷ Ν
    (β,Β), ξ ⟷ Ξ
    (γ,Γ), ο ⟷ Ο
    (δ,Δ), Π ⟷ Π
    (ϵ,Ε), ρ ⟷ Ρ
    (ζ,Ζ), σ ⟷ Σ
    (η,Η), τ ⟷ Τ
    (θ,Θ), υ ⟷ Υ
    (ι,Ι), ϕ ⟷ Φ
    (κ,Κ), χ ⟷ Χ
    (λ,Λ), ψ ↔ Ψ
    (μ,Μ), ω ⟷ Ω
end α Α β Β γ Γ δ Δ ϵ Ε ζ Ζ η Η θ Θ ι Ι κ Κ λ Λ μ Μ

# Networks containing t, im, and π should generate errors.

# test I works
rn = @reaction_network begin
    k1, S + I --> 2I
    k2, I --> R
end k1 k2
@test isequal(opname(species(rn)[2]),:I)


# test names work
rn = @reaction_network SIR1 begin
    k1, S + I --> 2I
    k2, I --> R
end k1 k2
@test nameof(rn) == :SIR1


