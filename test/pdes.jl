using Catalyst, OrdinaryDiffEq, Test
using ModelingToolkit, DomainSets
const MT = ModelingToolkit

@parameters k[1:7] D[1:3] n0[1:3] A
@variables t
@species x y U(x, y, t) V(x, y, t) W(x, y, t)
rxs = [Reaction(k[1], [U, W], [V, W]),
    Reaction(k[2], [V], [W], [2], [1]),
    Reaction(k[3], [W], [V], [1], [2]),
    Reaction(k[4], [U], nothing),
    Reaction(k[5], nothing, [U]),
    Reaction(k[6], [V], nothing),
    Reaction(k[7], nothing, [V])]
pars = vcat(MT.scalarize(k), MT.scalarize(D), MT.scalarize(n0), [A])
@named bpm = ReactionSystem(rxs, t, [U, V, W], pars; spatial_ivs = [x, y])

@test isequal(MT.get_iv(bpm), t)
@test issetequal(get_sivs(bpm), [x, y])
@test isspatial(bpm)

rxeqs = Catalyst.assemble_oderhs(bpm, states(bpm), combinatoric_ratelaws = false)
eqs = Dict((U => (k[5] - k[4] * U - k[1] * U * W),
            V => (2 * k[3] * W + k[1] * U * W + k[7] - k[6] * V - 2 * (V^2) * k[2]),
            W => ((V^2) * k[2] - k[3] * W)))
@test all(isequal.((MT.unwrap(eqs[st]) for st in states(bpm)), rxeqs))

@test issetequal(species(bpm), [MT.unwrap(U), MT.unwrap(V), MT.unwrap(W)])

# test a few API functions
ns = [-1 0 0 -1 1 0 0;
      1 -2 2 0 0 -1 1;
      0 1 -1 0 0 0 0]
@test ns == netstoichmat(bpm)
bpm2 = deepcopy(bpm)
@test bpm == bpm2

# check we can build a PDESystem
∂t = Differential(t)
∂x = Differential(x)
∂y = Differential(y)
Δ(u) = (∂x^2)(u) + (∂y^2)(u)
eqs = Vector{Equation}(undef, 3)
bcs = Vector{Equation}()
smap = speciesmap(bpm)
evalat(u, a, b, t) = (operation(ModelingToolkit.unwrap(u)))(a, b, t)
function icfun(n, x, y, A)
    float(rand(Poisson(round(n * A * 10))) / A / 10)
end
@register icfun(n, x, y, A)
L = 32.0
tstop = 5e4
for (i, st) in enumerate(states(bpm))
    idx = smap[st]
    eqs[i] = ∂t(st) ~ D[idx] * Δ(st) + rxeqs[idx]
    newbcs = [evalat(st, x, y, 0.0) ~ icfun(n0[idx], x, y, A),
        evalat(st, 0.0, y, t) ~ evalat(st, L, y, t),
        evalat(st, x, 0.0, t) ~ evalat(st, x, L, t)]
    append!(bcs, newbcs)
end
domains = [x ∈ Interval(0.0, L),
    y ∈ Interval(0.0, L),
    t ∈ Interval(0.0, tstop)]
pmap = collect(MT.defaults(bpm))
@named bpmpdes = PDESystem(eqs, bcs, domains, [x, y, t], [U, V, W], pmap)

rxs = [Reaction(k[1] * x, [U, W], [V, W]),
    Reaction(k[2] * y, [V], [W], [2], [1]),
    Reaction(k[3] + t, [W], [V], [1], [2]),
    Reaction(k[4], [U], nothing),
    Reaction(k[5], nothing, [U]),
    Reaction(k[6], [V], nothing),
    Reaction(k[7], nothing, [V])]
@named bpm4 = ReactionSystem(rxs, t, [U, V, W], pars; spatial_ivs = [x, y])
@test !ismassaction(rxs[1], bpm4)
@test ismassaction(rxs[1], bpm4; ivset = Set([t, y]))
@test !ismassaction(rxs[2], bpm4)
@test !ismassaction(rxs[3], bpm4)
@test ismassaction(rxs[4], bpm4)
