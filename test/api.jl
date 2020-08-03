using Catalyst, DiffEqBase, ModelingToolkit, Test

@parameters t k1 k2
@variables S I R
rxs = [Reaction(k1, [S,I], [I], [1,1], [2]),
       Reaction(k2, [I], [R]) ]
rs = ReactionSystem(rxs, t, [S,I,R], [k1,k2])

specset = Set([S.op =>1, I.op => 2, R.op => 3])
@test issetequal(specset, speciesmap(rs))

pset = Set([k1.op => 1, k2.op => 2])
@test issetequal(pset, paramsmap(rs))

rxs2 = [Reaction(k2, [I], [R], [1], [1]),
        Reaction(k1, [S,I], [I], [1,1], [2])]
rs2 = ReactionSystem(rxs2, t, [R,I,S], [k2,k1])
@test rs == rs2

rs3 = make_empty_network()
@parameters k3 k4
@variables D
addspecies!(rs3, S)
addspecies!(rs3, D.op)
addparam!(rs3, k3)
addparam!(rs3, k4.op)
@test issetequal(species(rs3), [S.op, D.op])
@test issetequal(params(rs3), [k3.op, k4.op])
addreaction!(rs3, Reaction(k3, [S], [D]))
addreaction!(rs3, Reaction(k4, [S,I], [D]))
merge!(rs, rs3)
addspecies!(rs2, S)
addspecies!(rs2, D.op)
addparam!(rs2, k3)
addparam!(rs2, k4.op)
addreaction!(rs2, Reaction(k3, [S], [D]))
addreaction!(rs2, Reaction(k4, [S,I], [D]))
@test rs2 == rs

rxs = [Reaction(k1*S, [S,I], [I], [2,3], [2]),
       Reaction(k2*R, [I], [R]) ]
rs = ReactionSystem(rxs, t, [S,I,R], [k1,k2])
deps = [s.op for s in dependents(rxs[2], rs)]
@test isequal(deps, [R.op,I.op])
addspecies!(rs, Variable(:S))
@test numspecies(rs) == 3
addspecies!(rs, Variable(:S), disablechecks=true)
@test numspecies(rs) == 4
addparam!(rs, Variable(:k1))
@test numparams(rs) == 2
addparam!(rs, Variable(:k1), disablechecks=true)
@test numparams(rs) == 3