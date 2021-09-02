using Catalyst, DiffEqBase, ModelingToolkit, Test
using SparseArrays
using ModelingToolkit: value

@parameters t k1 k2
@variables S(t) I(t) R(t)
rxs = [Reaction(k1, [S,I], [I], [1,1], [2]),
       Reaction(k2, [I], [R]) ]
@named rs = ReactionSystem(rxs, t, [S,I,R], [k1,k2])

specset = Set([value(S) =>1, value(I) => 2, value(R) => 3])
@test issetequal(specset, speciesmap(rs))

pset = Set([value(k1) => 1, value(k2) => 2])
@test issetequal(pset, paramsmap(rs))

rxs2 = [Reaction(k2, [I], [R], [1], [1]),
        Reaction(k1, [S,I], [I], [1,1], [2])]
@named rs2 = ReactionSystem(rxs2, t, [R,I,S], [k2,k1])
@test rs == rs2

rs3 = make_empty_network()
@parameters k3 k4
@variables D(t)
addspecies!(rs3, S)
addspecies!(rs3, D)
addparam!(rs3, k3)
addparam!(rs3, k4)
@test issetequal(species(rs3), [S, D])
@test issetequal(params(rs3), [k3, k4])
addreaction!(rs3, Reaction(k3, [S], [D]))
addreaction!(rs3, Reaction(k4, [S,I], [D]))
merge!(rs, rs3)
addspecies!(rs2, S)
addspecies!(rs2, D)
addparam!(rs2, k3)
addparam!(rs2, k4)
addreaction!(rs2, Reaction(k3, [S], [D]))
addreaction!(rs2, Reaction(k4, [S,I], [D]))
@test rs2 == rs

rxs = [Reaction(k1, [S,I], [I], [1,1], [2]),
       Reaction(k2, [I], [R]) ]
@named rs = ReactionSystem(rxs, t, [S,I,R], [k1,k2])
rs3 = make_empty_network()
addspecies!(rs3, S)
addspecies!(rs3, D)
addparam!(rs3, k3)
addparam!(rs3, k4)
addreaction!(rs3, Reaction(k3, [S], [D]))
addreaction!(rs3, Reaction(k4, [S,I], [D]))
rs4 = merge(rs, rs3)
@test rs2 == rs4

rxs = [Reaction(k1*S, [S,I], [I], [2,3], [2]),
       Reaction(k2*R, [I], [R]) ]
@named rs = ReactionSystem(rxs, t, [S,I,R], [k1,k2])
deps = dependents(rxs[2], rs)
@test isequal(deps, [R,I])
@test isequal(dependents(rxs[1], rs), dependants(rxs[1], rs))
addspecies!(rs, S)
@test numspecies(rs) == 3
addspecies!(rs, S, disablechecks=true)
@test numspecies(rs) == 4
addparam!(rs, k1)
@test numparams(rs) == 2
addparam!(rs, k1, disablechecks=true)
@test numparams(rs) == 3


rnmat = @reaction_network begin
       α, S + 2I --> 2I
       β, 3I --> 2R + S
   end α β

smat = [1 0;
        2 3;
        0 0]
pmat = [0 1;
        2 0;
        0 2]
@test smat == substoichmat(rnmat) == Matrix(substoichmat(rnmat, sparse=true))
@test pmat == prodstoichmat(rnmat) == Matrix(prodstoichmat(rnmat, sparse=true))
              
############## testing newly added intermediate complexes reaction networks##############
rns  = Vector{ReactionSystem}(undef,6)
# mass-action non-catalytic
rns[1] = @reaction_network begin
    k₁, 2A --> B
    k₂, A --> C
    k₃, C --> D
    k₄, B + D --> E
end k₁ k₂ k₃ k₄
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
rcs,B2 = reactioncomplexes(rns[1])
@test B == B2 == Matrix(reactioncomplexes(rns[1], sparse=true)[2])
@test Z == complexstoichmat(rns[1]) == Matrix(complexstoichmat(rns[1], sparse=true))
@test Δ == complexoutgoingmat(rns[1]) == Matrix(complexoutgoingmat(rns[1], sparse=true))
@test linkageclasses(incidencematgraph(B)) == linkageclasses(incidencematgraph(sparse(B))) == [[1,2],[3,4,5],[6,7]]
@test deficiency(rns[1]) == 0

# mass-action rober
rns[2] = @reaction_network begin
    k₁, A --> B
    k₂, B + B --> C + B
    k₃, B + C --> A + C
end k₁ k₂ k₃
Z = [1 0 0 0 1;
      0 1 2 1 0;
      0 0 0 1 1]
B = [-1 0 0;
      1 0 0;
      0 -1 0;
      0 1 -1;
      0 0 1]
Δ = [-1 0 0; 0 0 0; 0 -1 0; 0 0 -1; 0 0 0]
rcs,B2 = reactioncomplexes(rns[2])
@test B == B2 == Matrix(reactioncomplexes(rns[2], sparse=true)[2])
@test Z == complexstoichmat(rns[2]) == Matrix(complexstoichmat(rns[2], sparse=true))
@test Δ == complexoutgoingmat(rns[2]) == Matrix(complexoutgoingmat(rns[2], sparse=true))
@test linkageclasses(incidencematgraph(B))==linkageclasses(incidencematgraph(sparse(B))) == [[1,2],[3,4,5]]
@test deficiency(rns[2]) == 1



#  some rational functions as rates
rns[3] = @reaction_network begin
    k₁ , ∅ --> X₁
    ( k₂/(1 + X₁*X₂ + X₃*X₄ ), k₃/(1 + X₁*X₂ + X₃*X₄ )), 2X₁ + X₂ ↔ 3X₃ + X₄
    k₄, X₄ --> ∅
end k₁ k₂ k₃ k₄
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
rcs,B2 = reactioncomplexes(rns[3])
@test B == B2 == Matrix(reactioncomplexes(rns[3], sparse=true)[2])
@test Z == complexstoichmat(rns[3]) == Matrix(complexstoichmat(rns[3], sparse=true))
@test Δ == complexoutgoingmat(rns[3]) ==  Matrix(complexoutgoingmat(rns[3], sparse=true))
@test linkageclasses(incidencematgraph(B))==linkageclasses(incidencematgraph(sparse(B))) == [[1,2,5],[3,4]]
@test deficiency(rns[3]) == 0


# repressilator
rns[4]  = @reaction_network begin
   hillr(P₃,α,K,n), ∅ --> m₁
   hillr(P₁,α,K,n), ∅ --> m₂
   hillr(P₂,α,K,n), ∅ --> m₃
   (δ,γ), m₁ ↔ ∅
   (δ,γ), m₂ ↔ ∅
   (δ,γ), m₃ ↔ ∅
   β, m₁ --> m₁ + P₁
   β, m₂ --> m₂ + P₂
   β, m₃ --> m₃ + P₃
   μ, P₁ --> ∅
   μ, P₂ --> ∅
   μ, P₃ --> ∅
end α K n δ γ β μ;
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
Δ = [-1 -1 -1 0 -1 0 -1 0 -1 0 0 0 0 0 0; 0 0 0 -1 0 0 0 0 0 -1 0 0 0 0 0; 0 0 0 0 0 -1 0 0 0 0 -1 0 0 0 0;
       0 0 0 0 0 0 0 -1 0 0 0 -1 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
       0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1]
rcs,B2 = reactioncomplexes(rns[4])
@test B == B2 == Matrix(reactioncomplexes(rns[4], sparse=true)[2])
@test Z == complexstoichmat(rns[4]) == Matrix(complexstoichmat(rns[4], sparse=true))
@test Δ == complexoutgoingmat(rns[4]) == Matrix(complexoutgoingmat(rns[4], sparse=true))
@test linkageclasses(incidencematgraph(B)) ==linkageclasses(incidencematgraph(sparse(B)))==  [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]]
@test deficiency(rns[4]) == 3


#brusselator
rns[5] = @reaction_network begin
   A, ∅ → X
   1, 2X + Y → 3X
   B, X → Y
   1, X → ∅
end A B
Z = [0 1 2 3 0;
      0 0 1 0 1]
B = [-1 0 0 1;
      1 0 -1 -1;
      0 -1 0 0;
      0 1 0 0;
      0 0 1 0]
Δ = [-1 0 0 0; 0 0 -1 -1; 0 -1 0 0; 0 0 0 0; 0 0 0 0]
rcs,B2 = reactioncomplexes(rns[5])
@test B == B2 == Matrix(reactioncomplexes(rns[5], sparse=true)[2])
@test Z == complexstoichmat(rns[5]) == Matrix(complexstoichmat(rns[5], sparse=true))
@test Δ == complexoutgoingmat(rns[5]) == Matrix(complexoutgoingmat(rns[5], sparse=true))
@test linkageclasses(incidencematgraph(B))==linkageclasses(incidencematgraph(sparse(B))) == [[1,2,5],[3,4]]
@test deficiency(rns[5]) == 1

# some rational functions as rates
rns[6] = @reaction_network begin
    (k₁, k₋₁), X₁ + X₂ <--> X₃ + 2X₄
    (k₂/(1 + X₄*X₅ + X₆*X₇), k₋₂/(1 + X₄*X₅ + X₆*X₇)), 3X₄ + X₅ <--> X₆ + X₇
    (k₃/(1 + X₇ + X₈ + X₉ + X₁₀) , k₋₃/(1 + X₇ + X₈ + X₉ + X₁₀ )), 5X₇ + X₈ <--> X₉ + X₁₀
end k₁ k₋₁ k₂ k₋₂ k₃ k₋₃
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
rcs,B2 = reactioncomplexes(rns[6])
@test B == B2 == Matrix(reactioncomplexes(rns[6], sparse=true)[2])
@test Z == complexstoichmat(rns[6]) == Matrix(complexstoichmat(rns[6], sparse=true))
@test Δ == complexoutgoingmat(rns[6]) == Matrix(complexoutgoingmat(rns[6], sparse=true))
@test linkageclasses(incidencematgraph(B))==linkageclasses(incidencematgraph(sparse(B))) == [[1,2],[3,4],[5,6]]
@test deficiency(rns[6]) == 0
                     

reaction_networks_standard = Vector{ReactionSystem}(undef,10)
reaction_networks_hill = Vector{ReactionSystem}(undef,10)
reaction_networks_real = Vector{ReactionSystem}(undef,3)
### some more networks to test whether Z*B == netstoichmat(rn) or not. ###
reaction_networks_standard[1] = @reaction_network begin
    (p1,p2,p3), ∅ → (X1,X2,X3)
    (k1,k2), X2 ⟷ X1 + 2X3
    (k3,k4), X1 ⟷ X3
    (d1,d2,d3), (X1,X2,X3) → ∅
end p1 p2 p3 k1 k2 k3 k4 d1 d2 d3

reaction_networks_standard[2] = @reaction_network begin
    mmr(X2,v1,K1), ∅ → X1
    mm(X1,v2,K2), ∅ → X2
    d, X1+X2 → ∅
end v1 K1 v2 K2 d

reaction_networks_standard[3] = @reaction_network begin
    mm(X2,v1,K1), ∅ → X1
    mm(X3,v2,K2), ∅ → X2
    (k1,k2), X1 ⟷ X3
    (k3,k4), X3 + X2 ⟷ X4 + X1
    d, (X1,X2,X3,X4) → ∅
end v1 K1 v2 K2 k1 k2 k3 k4 d

reaction_networks_standard[4] = @reaction_network begin
    mmr(X4,v1,K1), ∅ → X1
    mmr(X1,v2,K2), ∅ → X2
    mmr(X2,v3,K3), ∅ → X3
    mmr(X3,v4,K4), ∅ → X4
    (d1,d2,d3,d4), (X1,X2,X3,X4) → ∅
end v1 K1 v2 K2 v3 K3 v4 K4 d1 d2 d3 d4

reaction_networks_standard[5] = @reaction_network begin
    p, ∅ → X1
    (k1,k2), X1 ⟷ X2
    (k3,k4), X2 ⟷ X3
    (k5,k6), X3 ⟷ X4
    d, X4 → ∅
end p k1 k2 k3 k4 k5 k6 d

reaction_networks_standard[6] = @reaction_network begin
    (p1,p2), ∅ → (X1,X2)
    (k1,k2), 2X1 ⟷ X3
    (k3,k4), X2 + X3 ⟷ 4X4
    (k5,k6), X4 + X1 ⟷ 2X3
    d, (X1,X2,X3,X4) → ∅
end p1 p2 k1 k2 k3 k4 k5 k6 d

reaction_networks_standard[7] = @reaction_network begin
    (p1,p2,p3), ∅ → (X1,X2,X3)
    (k1,k2), X1 + 2X2 ⟷ X4
    (mm(X3,v1,K1),k3), X4 ⟷ X5
    (d1,d2,d3,d4,d5), (X1,X2,X3,X4,X5) → ∅
end p1 p2 p3 k1 k2 k3 v1 K1 d1 d2 d3 d4 d5

reaction_networks_standard[8] = @reaction_network begin
    p, ∅ → 2X1
    k1, X1 → X2
    (k2, k3), X2 → X3
    d, X3 → ∅
end p k1 k2 k3 d

reaction_networks_standard[9] = @reaction_network begin
    (p1,p2,p3), ∅ ⟶ (X1,X2,X3)
    (d1,d2,d3), (X1,X2,X3) ⟶ ∅
    (k1,k2), X1 + X2 ⟷ X3
    (k3,k4), X3 ⟷ X4
    (k5,k6), X4 ⟷ X1 + X2
end p1 p2 p3 d1 d2 d3 k1 k2 k3 k4 k5 k6

reaction_networks_standard[10] = @reaction_network begin
    p, ∅ ⟶ X1
    (k1,k2), X1 → X2
    (k3,k4), X2 → X3
    (k5,k6), X3 → X4
    (k7,k8), X4 → X5
    d, X5  ⟶ ∅
end p k1 k2 k3 k4 k5 k6 k7 k8 d


### Network with hill functions ###.

reaction_networks_hill[1] = @reaction_network begin
    hillr(X2,v1,K1,n1),   ∅ → X1
    hillr(X1,v2,K2,n2),   ∅ → X2
    (d1,d2),          (X1,X2) → ∅
end v1 v2 K1 K2 n1 n2 d1 d2

reaction_networks_hill[2] = @reaction_network begin
    hillr(X3,v1,K1,n1),   ∅ → X1
    hillr(X1,v2,K2,n2),   ∅ → X2
    hillr(X2,v3,K3,n3),   ∅ → X3
    (d1,d2,d3),          (X1,X2,X3) → ∅
end v1 v2 v3 K1 K2 K3 n1 n2 n3 d1 d2 d3

reaction_networks_hill[3] = @reaction_network begin
    hillr(X2,v1,K1,n1), ∅ → X1
    hill(X1,v2,K2,n2), ∅ → X2
    d, X1+X2 → ∅
end v1 K1 n1 v2 K2 n2 d

reaction_networks_hill[4] = @reaction_network begin
    hillr(X2,v1,K1,n1)*hillr(X3,v1,K1,n1), ∅ → X1
    hillr(X1,v2,K2,n2)*hillr(X3,v2,K2,n2), ∅ → X2
    hillr(X1,v3,K3,n3)*hillr(X2,v3,K3,n3), ∅ → X3
    (d1,d2,d3), (X1,X2,X3)  ⟶ ∅
end v1 K1 n1 v2 K2 n2 v3 K3 n3 d1 d2 d3

reaction_networks_hill[5] = @reaction_network begin
    hillr(X2,v1,K1,n1)*hill(X4,v1,K1,n1), ∅ → X1
    hill(X5,v2,K2,n2), ∅ → X2
    hill(X3,v3,K3,n3), ∅ → X3
    hillr(X1,v4,K4,n4), ∅ → X4
    hill(X2,v5,K5,n5), ∅ → X5
    (k1,k2), X2 ⟷ X1 + 2X4
    (k3,k4), X4 ⟷ X3
    (k5,k6), 3X5 + X1 ⟷ X2
    (d1,d2,d3,d4,d5), (X1,X2,X3,X4,X5)  ⟶ ∅
end v1 K1 n1 v2 K2 n2 v3 K3 n3 v4 K4 n4 v5 K5 n5 k1 k2 k3 k4 k5 k6 d1 d2 d3 d4 d5

reaction_networks_hill[6] = @reaction_network begin
    v/10+hill(X1,v,K,n), ∅ → X1
    d, X1 → ∅
end v K n d

reaction_networks_hill[7] = @reaction_network begin
    v/10 + hill(X1,v,K,n), ∅ → X1 + X2
    (k1,k2), X1 + X2 ↔ X3
    k3, X3 → X1
    d, (X1,X2,X3) → ∅
end v K n k1 k2 k3 d

reaction_networks_hill[8] = @reaction_network begin
    hill(X2,v1,K1,n1), ∅ → X1
    hillr(X1,v2,K2,n2)*hill(X3,v3,K3,n3), ∅ → X2
    hill(X2,v4,K4,n4), ∅ → X3
    (d1,d2,d3), (X1,X2,X3) → ∅
end v1 K1 n1 v2 K2 n2 v3 K3 n3 v4 K4 n4 d1 d2 d3

reaction_networks_hill[9] = @reaction_network begin
    hill(X1,v1,K1,n1)*hillr(X1,v2,K2,n2), ∅ → X1
    d, X1 → ∅
end v1 K1 n1 v2 K2 n2 d

reaction_networks_hill[10] = @reaction_network begin
    hill(X2,v1,K1,n1), ∅ → X1
    hillr(X4,v2,K2,n2), ∅ → X2
    (k1,k2), 2X1 + X2 ⟷ X3
    (k3,k4), 2X2 + X3 ⟷ X4
    (k5,k6), X1 + X2 + X3 + X4 ⟷ X5 + X6
    (d1,d2), (X5,X6) → ∅
end v1 K1 n1 v2 K2 n2 k1 k2 k3 k4 k5 k6 d1 d2


reaction_networks_real[1] = @reaction_network begin
    v0 + hill(σ,v,K,n), ∅ → (σ+A)
    deg, (σ,A,Aσ) → ∅
    (kB,kD), A + σ ↔ Aσ
    S*kC, Aσ → σ
end v0 v K n kD kB kC deg S;

# A cell cycle model
reaction_networks_real[2] = @reaction_network begin
  k1, 0 --> Y
  k2p, Y --> 0
  k2pp*P, Y --> 0
  (k3p+k3pp*A)/(J3+Po), Po-->P
  (k4*m)/(J4+P), Y + P --> Y + Po
end k1 k2p k2pp k3p k3pp A J3 k4 m J4
#@add_constraint cc_network P+Po=1
#reaction_networks_real[3] = cc_network

# A bistable switch
reaction_networks_real[3] = @reaction_network begin
    d,    (X,Y) → ∅
    hillr(Y,v1,K1,n1), ∅ → X
    hillr(X,v2,K2,n2), ∅ → Y
end d v1 K1 n1 v2 K2 n2

myrn = [reaction_networks_standard;reaction_networks_hill;reaction_networks_real]
for i in 1:length(myrn)
    local rcs,B = reactioncomplexes(myrn[i])
    @test B == Matrix(reactioncomplexes(myrn[i], sparse=true)[2])
    local Z = complexstoichmat(myrn[i])
    @test Z == Matrix(complexstoichmat(myrn[i], sparse=true))
    @test Z*B == netstoichmat(myrn[i]) == Matrix(netstoichmat(myrn[i], sparse=true))
end


