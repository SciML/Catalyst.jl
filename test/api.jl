using Catalyst, DiffEqBase, ModelingToolkit, Test

using ModelingToolkit: value

@parameters t k1 k2
@variables S(t) I(t) R(t)
rxs = [Reaction(k1, [S,I], [I], [1,1], [2]),
       Reaction(k2, [I], [R]) ]
rs = ReactionSystem(rxs, t, [S,I,R], [k1,k2])

specset = Set([value(S) =>1, value(I) => 2, value(R) => 3])
@test issetequal(specset, speciesmap(rs))

pset = Set([value(k1) => 1, value(k2) => 2])
@test issetequal(pset, paramsmap(rs))

rxs2 = [Reaction(k2, [I], [R], [1], [1]),
        Reaction(k1, [S,I], [I], [1,1], [2])]
rs2 = ReactionSystem(rxs2, t, [R,I,S], [k2,k1])
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
rs = ReactionSystem(rxs, t, [S,I,R], [k1,k2])
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
rs = ReactionSystem(rxs, t, [S,I,R], [k1,k2])
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

smat = [1 2 0;
        0 3 0]
pmat = [0 2 0;
        1 0 2]
@test all(smat .== substoichmat(rnmat))
@test all(pmat .== prodstoichmat(rnmat))
