### Fetch required packages and reaction networks ###
using DiffEqBase, DiffEqBiological, Random, Test, UnPack
include("test_networks.jl")


### Debugg functions ###

# Gets all the reactants in a set of equations.
function all_reactants(eqs)
    all_reactants = []
    for eq in eqs
        append!(all_reactants,getproperty.(getproperty.(eq.substrates,:op),:name))
        append!(all_reactants,getproperty.(getproperty.(eq.products,:op),:name))
    end
    return Set{Symbol}(unique(all_reactants))
end

# Gets all parameters (where every reaction rate is constant)
function all_parameters(eqs)
    return Set(unique(map(eq -> eq.rate.op.name,eqs)))
end
using UnPack


### Test basic properties of networks ###
@unpack eqs,iv,states,ps,name,systems = reaction_networks_standard[1]
@test length(eqs) == 10
@test iv.name == :t
@test length(states) == 3
@test all(getproperty.(states,:name) .== [:X1,:X2,:X3])
@test all_reactants(eqs) == Set([:X1,:X2,:X3])
@test length(ps) == 10
@test all(getproperty.(ps,:name) .== [:p1,:p2,:p3,:k1,:k2,:k3,:k4,:d1,:d2,:d3])
@test all_parameters(eqs) == Set([:p1,:p2,:p3,:k1,:k2,:k3,:k4,:d1,:d2,:d3])

@unpack eqs,iv,states,ps,name,systems = reaction_networks_standard[2]
@test length(eqs) == 3
@test iv.name == :t
@test length(states) == 2
@test all(getproperty.(states,:name) .== [:X1,:X2])
@test all_reactants(eqs) == Set([:X1,:X2])
@test length(ps) == 5
@test all(getproperty.(ps,:name) .== [:v1,:K1,:v2,:K2,:d])

@unpack eqs,iv,states,ps,name,systems = reaction_networks_standard[3]
@test length(eqs) == 10
@test iv.name == :t
@test length(states) == 4
@test all(getproperty.(states,:name) .== [:X1,:X2,:X3,:X4])
@test all_reactants(eqs) == Set([:X1,:X2,:X3,:X4])
@test length(ps) == 9
@test all(getproperty.(ps,:name) .== [:v1,:K1,:v2,:K2,:k1,:k2,:k3,:k4,:d])

@unpack eqs,iv,states,ps,name,systems = reaction_networks_standard[4]
@test length(eqs) == 8
@test iv.name == :t
@test length(states) == 4
@test all(getproperty.(states,:name) .== [:X1,:X2,:X3,:X4])
@test all_reactants(eqs) == Set([:X1,:X2,:X3,:X4])
@test length(ps) == 12
@test all(getproperty.(ps,:name) .== [:v1,:K1,:v2,:K2,:v3,:K3,:v4,:K4,:d1,:d2,:d3,:d4])

@unpack eqs,iv,states,ps,name,systems = reaction_networks_standard[5]
@test length(eqs) == 8
@test iv.name == :t
@test length(states) == 4
@test all(getproperty.(states,:name) .== [:X1,:X2,:X3,:X4])
@test all_reactants(eqs) == Set([:X1,:X2,:X3,:X4])
@test length(ps) == 8
@test all(getproperty.(ps,:name) .== [:p,:k1,:k2,:k3,:k4,:k5,:k6,:d])
@test all_parameters(eqs) == Set([:p,:k1,:k2,:k3,:k4,:k5,:k6,:d])

@unpack eqs,iv,states,ps,name,systems = reaction_networks_hill[1]
@test length(eqs) == 4
@test iv.name == :t
@test length(states) == 2
@test all(getproperty.(states,:name) .== [:X1,:X2])
@test all_reactants(eqs) == Set([:X1,:X2])
@test length(ps) == 8
@test all(getproperty.(ps,:name) .== [:v1,:v2,:K1,:K2,:n1,:n2,:d1,:d2])

@unpack eqs,iv,states,ps,name,systems = reaction_networks_constraint[1]
@test length(eqs) == 6
@test iv.name == :t
@test length(states) == 3
@test all(getproperty.(states,:name) .== [:X1,:X2,:X3])
@test all_reactants(eqs) == Set([:X1,:X2,:X3])
@test length(ps) == 6
@test all(getproperty.(ps,:name) .== [:k1,:k2,:k3,:k4,:k5,:k6])

@unpack eqs,iv,states,ps,name,systems = reaction_networks_real[1]
@test length(eqs) == 4
@test iv.name == :t
@test length(states) == 2
@test all(getproperty.(states,:name) .== [:X,:Y])
@test all_reactants(eqs) == Set([:X,:Y])
@test length(ps) == 2
@test all(getproperty.(ps,:name) .== [:A,:B])

@unpack eqs,iv,states,ps,name,systems = reaction_networks_weird[1]
@test length(eqs) == 2
@test iv.name == :t
@test length(states) == 1
@test all(getproperty.(states,:name) .== [:X])
@test all_reactants(eqs) == Set([:X])
@test length(ps) == 2
@test all(getproperty.(ps,:name) .== [:p,:d])

@unpack eqs,iv,states,ps,name,systems = reaction_networks_weird[2]
@test length(eqs) == 4
@test iv.name == :t
@test length(states) == 3
@test all(getproperty.(states,:name) .== [:X,:Y,:Z])
@test all_reactants(eqs) == Set([:X,:Y,:Z])
@test length(ps) == 4
@test all(getproperty.(ps,:name) .== [:k1,:k2,:k3,:k4])


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
    mmR(X2,v1,K1), ∅ → X1
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
    mmR(X4,v1,K1), ∅ ⥟ X1
    mmR(X1,v2,K2), ∅ ⥟ X2
    mmR(X2,v3,K3), ∅ ⇀ X3
    mmR(X3,v4,K4), ∅ ⇁ X4
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

@test_broken if false # Causes weird error, see ModelingToolkit issue #450.
    for networks in identical_networks_1
        f1 = ODEFunction(convert(ODESystem,networks[1]),jac=true)
        f2 = ODEFunction(convert(ODESystem,networks[2]),jac=true)
        g1 = SDEFunction(convert(SDESystem,networks[1]))
        g2 = SDEFunction(convert(SDESystem,networks[2]))
        for factor in [1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3]
            u0 = factor*rand(length(networks[1].states))
            p = factor*rand(length(networks[1].ps))
            t = rand()
            @test all(abs.(f1(u0,p,t) .- f2(u0,p,t)) .< 100*eps())
            @test all(abs.(f1.jac(u0,p,t) .- f2.jac(u0,p,t)) .< 100*eps())
            @test all(abs.(g1(u0,p,t) .- g2(u0,p,t)) .< 100*eps())
        end
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

@test_broken if false # Causes weird error, see ModelingToolkit issue #450.
    for networks in identical_networks_2
        f1 = ODEFunction(convert(ODESystem,networks[1]),jac=true)
        f2 = ODEFunction(convert(ODESystem,networks[2]),jac=true)
        g1 = SDEFunction(convert(SDESystem,networks[1]))
        g2 = SDEFunction(convert(SDESystem,networks[2]))
        for factor in [1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3]
            u0 = factor*rand(length(networks[1].states))
            p = factor*rand(length(networks[1].ps))
            t = rand()
            @test all(abs.(f1(u0,p,t) .- f2(u0,p,t)) .< 100*eps())
            @test all(abs.(f1.jac(u0,p,t) .- f2.jac(u0,p,t)) .< 100*eps())
            @test all(abs.(g1(u0,p,t) .- g2(u0,p,t)) .< 100*eps())
        end
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

@test_broken if false # Causes weird error, see ModelingToolkit issue #450.
    for (i,networks) in enumerate(identical_networks_3)
        f1 = ODEFunction(convert(ODESystem,networks[1]),jac=true)
        f2 = ODEFunction(convert(ODESystem,networks[2]),jac=true)
        g1 = SDEFunction(convert(SDESystem,networks[1]))
        g2 = SDEFunction(convert(SDESystem,networks[2]))
        for factor in [1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3]
            u0 = factor*rand(length(networks[1].states))
            t = rand()
            @test all(abs.(f1(u0,parameter_sets[i],t) .- f2(u0,[],t)) .< 100*eps())
            @test all(abs.(f1.jac(u0,parameter_sets[i],t) .- f2.jac(u0,[],t)) .< 100*eps())
            @test all(abs.(g1(u0,parameter_sets[i],t) .- g2(u0,[],t)) .< 100*eps())
        end
    end
end


### Tests that time is handled properly ###

time_network = @reaction_network begin
    (t,k2), X1 ↔ X2
    (k3,t), X2 ↔ X3
    (t,k6), X3 ↔ X1
end k2 k3 k6

@test_broken if false # Causes weird error, see ModelingToolkit issue #450.
    f1 = ODEFunction(convert(ODESystem,reaction_networks_constraint[1]),jac=true)
    f2 = ODEFunction(convert(ODESystem,time_network),jac=true)
    g1 = SDEFunction(convert(SDESystem,reaction_networks_constraint[1]))
    g2 = SDEFunction(convert(SDESystem,time_network))
    for factor in [1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3]
        u0 = factor*rand(length(time_network.states))
        k2 = factor*rand(); k3 = factor*rand(); k6 = factor*rand();
        t = rand()
        p1 = [t, k2, k3, t, t, k6]; p2 = [k2, k3, k6];
        @test all(abs.(f1(u0,p1,t) .- f2(u0,p2,t)) .< 100*eps())
        @test all(abs.(f1.jac(u0,p1,t) .- f2.jac(u0,p2,t)) .< 100*eps())
        @test all(abs.(g1(u0,p1,t) .- g2(u0,p2,t)) .< 100*eps())
    end
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

# Networks containing t, I, im, and π should generate errors.
