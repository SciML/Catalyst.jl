# Partial Differential Equation Models
**Note** this functionality is work in progress in both Catalyst, ModelingToolkit, and MethodOfLines (generally across SciML). As such, the recommended workflows, API features, and observed simulation performance may change.

We'll simulate the Baras–Pearson–Mansour (BPM) pattern formation model using the parameters in (Kim et al., J. Chem. Phys., 146, 2017).

First we load the packages we'll use
```julia
using Catalyst, MethodOfLines, DomainSets, OrdinaryDiffEq, Plots, Random
using ModelingToolkit: scalarize
```

Next let's specify our default parameter values
```julia
# the reaction rates, k₁,...,k₇
rates = [2e-4, 2e-4, 1.0, 3.33e-3, 16.7, 3.67e-2, 4.44]

# diffusivities [DU, DV, DW]
Dvals = [.1, .01, .01]

# physical domain size is LxL
L = 32.0

# time to end the simulation
tstop = 5e4

# for calculating the initial condition
icfun(n, x, y) = n * (1 + .01 * randn())
@register icfun(n, x, y)
```

We now define the reaction model
```julia
@parameters begin
    k[1:7] = rates
    D[1:3] = Dvals
    n0[1:3] = [1686.0, 534.0, 56.4]    # average initial densities
end
@variables t x y U(x,y,t) V(x,y,t) W(x,y,t)
rxs = [Reaction(k[1], [U, W], [V, W]),
       Reaction(k[2], [V], [W], [2], [1]),
       Reaction(k[3], [W], [V], [1], [2]),
       Reaction(k[4], [U], nothing),
       Reaction(k[5], nothing, [U]),
       Reaction(k[6], [V], nothing),
       Reaction(k[7], nothing, [V])]
pars = vcat(scalarize(k), scalarize(D), scalarize(n0))
@named bpm = ReactionSystem(rxs, t, [U, V, W], pars; spatial_ivs = [x,y])
```

We now put together the symbolic PDE model
```julia
# get the reaction terms
rxeqs = Catalyst.assemble_oderhs(bpm, states(bpm))

# get the ordering of the variables within rxeqs
smap = speciesmap(bpm)

# for defining boundary conditions
evalat(u, a, b, t) = (operation(ModelingToolkit.unwrap(u)))(a, b, t)

# construct the PDEs and bcs
∂t = Differential(t)
∂x = Differential(x)
∂y = Differential(y)
Δ(u) = (∂x^2)(u) + (∂y^2)(u)
eqs = Vector{Equation}(undef, 3)
bcs = Vector{Equation}()
for (i,st) in enumerate(states(bpm))
    idx = smap[st]
    eqs[i] = ∂t(st) ~ D[idx] * Δ(st) + rxeqs[idx]
    newbcs = [evalat(st, x, y, 0) ~ icfun(n0[idx], x, y),
              evalat(st, 0, y, t) ~ evalat(st, L, y, t),
              evalat(st, x, 0, t) ~ evalat(st, x, L, t)]
    append!(bcs, newbcs)
end

# define the domains
domains = [x ∈ Interval(0, L), y ∈ Interval(0, L), t ∈ Interval(0, tstop)]

@named bpmpdes = PDESystem(eqs, bcs, domains, [x,y,t], [U, V, W], pars)
```

Discretization
```julia
N = 3     # number of mesh points
h = L / N  # mesh width
order = 2  # order of the discretization

discretization = MOLFiniteDifference([x => h, y => h], t; approx_order = order,
                                     grid_align = center_align)
prob = discretize(bpmpdes, discretization)
sol = solve(prob, TRBDF2(), saveat = (tstop/10))
```

Plotting
```julia

```