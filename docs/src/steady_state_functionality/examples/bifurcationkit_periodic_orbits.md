# [Computing Periodic Orbits (Oscillations) Using BifurcationKit.jl](@id bifurcationkit_periodic_orbits)
Previously, we have shown how to [compute bifurcation diagrams](@ref bifurcation_diagrams) using [BifurcationKit.jl](https://github.com/bifurcationkit/BifurcationKit.jl). In this example we will consider a system which exhibits an oscillation and show how to use BifurcationKit to track not just the system's (potentially unstable) steady state, but also the periodic orbit itself. More information on how to track periodic orbits can be found in the [BifurcationKit documentation](https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/tutorials/tutorials/#Periodic-orbits).

## [Computing the bifurcation diagram for the Repressilator](@id bifurcationkit_periodic_orbits_bifdia)
We will first compute the bifurcation diagram, using the same approach as in the [corresponding tutorial](@ref bifurcation_diagrams). For this example, we will use the oscillating [Repressilator](@ref basic_CRN_library_repressilator) model.
```@example bifurcationkit_periodic_orbits
using Catalyst
repressilator = @reaction_network begin
    hillr(Z,v,K,n), ∅ --> X
    hillr(X,v,K,n), ∅ --> Y
    hillr(Y,v,K,n), ∅ --> Z
    d, (X, Y, Z) --> ∅
end
```
Next, we create a `BifurcationProblem` for our model. We will compute the bifurcation diagram with respect to the parameter $v$, and plot the species $X$ in the diagram.
```@example bifurcationkit_periodic_orbits
using BifurcationKit
bif_par = :v
u_guess = [:X => 20.0, :Y => 15.0, :Z => 15.0]
p_start = [:v => 10.0, :K => 15.0, :n => 3, :d => 0.2]
plot_var = :X
bprob = BifurcationProblem(repressilator, u_guess, p_start, bif_par; plot_var)
nothing # hide
```
We compute the bifurcation diagram using the `bifurcationdiagram` function. We will compute it across the interval $v \in (0,20)$.
```@example bifurcationkit_periodic_orbits
v_span = (0.0, 20.0)
opts_br = ContinuationPar(p_min = v_span[1], p_max = v_span[2])
bifdia = bifurcationdiagram(bprob, PALC(), 3, opts_br; bothside = true)
nothing # hide
```
Finally, we plot the bifurcation diagram.
```@example bifurcationkit_periodic_orbits
using Plots
plot(bifdia; xguide = bif_par, yguide = plot_var, xlimit = v_span, branchlabel = "Steady state concentration",
    linewidthstable = 5)
```
We note that for small values of $v$ the system's single steady state is stable (where the line is thicker). After a [Hopf](https://en.wikipedia.org/wiki/Hopf_bifurcation) bifurcation (the red point), the state turns unstable (where the line is thinner). For chemical reaction networks (which mostly are well-behaved) a single unstable steady state typically corresponds to an oscillation. We can confirm that the system oscillates in the unstable region (while it reaches a stable steady state in the stable region) using simulations:
```@example bifurcationkit_periodic_orbits
using OrdinaryDiffEqDefault
p_nosc = [:v => 5.0, :K => 15.0, :n => 3, :d => 0.2]
p_osc = [:v => 15.0, :K => 15.0, :n => 3, :d => 0.2]
prob_nosc = ODEProblem(repressilator, u_guess, 80.0, p_nosc)
prob_osc = ODEProblem(repressilator, u_guess, 80.0, p_osc)
sol_nosc = OrdinaryDiffEqDefault.solve(prob_nosc)
sol_osc = OrdinaryDiffEqDefault.solve(prob_osc)
plot(plot(sol_nosc; title = "v = 5"), plot(sol_osc; title = "v = 15"), size = (1000,400), lw = 4)
```

## [Tracking the periodic orbits](@id bifurcationkit_periodic_orbits_pos)
Next, we will use BifurcationKit.jl's [`continuation` function](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/library/#BifurcationKit.continuation) (the [`bifurcationdiagram` function](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/library/#BifurcationKit.bifurcationdiagram), which we previously have used, works by calling `continuation` recursively) to track the periodic orbit which appears with the Hopf bifurcation point.

First, we set the options for the continuation. Just like for bifurcation diagrams we must set our [continuation parameters](@ref bifurcation_diagrams_continuationpar). Here we will use the same one as for the initial diagram (however, additional ones can be supplied).
```@example bifurcationkit_periodic_orbits
opts_po = ContinuationPar(opts_br)
nothing # hide
```
During the continuation we will compute the periodic orbit using the [Collocation](https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/periodicOrbit/#Collocation-method) method (however, the [Trapezoid](https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/periodicOrbit/#Trapezoid-method) and [Shooting](https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/periodicOrbit/#Shooting-method) method also exist, with each having their advantages and disadvantages). For this, we create a [`PeriodicOrbitOCollProblem`](https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/library/#BifurcationKit.PeriodicOrbitOCollProblem) which we will supply to our continuation computation.
```@example bifurcationkit_periodic_orbits
poprob = PeriodicOrbitOCollProblem(50, 4)
nothing # hide
```
Finally, we will also create a `record_from_solution` function. This is a function which records information of the solution at each step of the continuation (which we can later plot or investigate using other means).
```@example bifurcationkit_periodic_orbits
X_idx = findfirst(isequal(repressilator.X), unknowns(complete(make_rre_algeqs(repressilator))))
function record_from_solution(x, p; kwargs...)
 xtt = get_periodic_orbit(p.prob, x, p.p)
 min, max = extrema(xtt[1,:])
 period = getperiod(p.prob, x, p.p)
    return (; min, max, period)
end
nothing # hide
```
Here, `get_periodic_orbit` computes the system's periodic orbit. From it we extract the $X$ species's minimum and maximum values (using `extrema`) and the period length (using `getperiod`). We return these quantities as a [`NamedTuple`](https://docs.julialang.org/en/v1/base/base/#Core.NamedTuple).

We can now compute the periodic orbit. First we must extract a branch from our bifurcation diagram (using [`get_branch`](https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/library/#BifurcationKit.get_branch)), as `continuation` do not work on bifurcation diagrams directly (our bifurcation diagram consists of a single branch, which we here extract). Next, we can call `continuation` on it, designating that we wish to do continuation from the second point on the branch (which corresponds to the Hopf bifurcation, the first point is one of the branch's two endpoints).
```@example bifurcationkit_periodic_orbits
branch = get_branch(bifdia, ()).γ
br_po = continuation(branch, 2, opts_po, poprob; record_from_solution)
nothing # hide
```
Finally, we can plot the periodic orbit. We use the `vars` argument to designate what we wish to plot. Here we first provide `(:param, :min)` (designating that we wish to plot the `min` value returned by `record_from_solution`, i.e. the minimum value throughout the periodic orbit) against the continuation parameter's ($v$) value. Next, we plot the maximum periodic orbit value using `(:param, :max)`.
```@example bifurcationkit_periodic_orbits
plot(bifdia; xlimit = v_span, branchlabel = "Steady state concentration", linewidthstable = 5)
plot!(br_po, vars = (:param, :min); color = 2, branchlabel = "Oscillation amplitude (min/max)")
plot!(br_po, vars = (:param, :max); color = 2, xguide = bif_par, yguide = plot_var, branchlabel = "")
```
Here we can see that, as $v$ increases, the oscillation amplitude increases with it.

Previously, we had `record_from_solution` record the periodic orbit's period. This means that we can plot it as well. Here, we plot it against $v$ using `vars = (:param, :period)`.
```@example bifurcationkit_periodic_orbits
plot(br_po, vars = (:param, :period); xguide = bif_par, yguide = "Period length", xlimit = v_span, ylimit = (0.0, Inf))
```
In the plot we see that the period starts at around $18$ time units, and slowly increase with $v$.


---
## [Citation](@id bifurcationkit_periodic_orbits_citation)
If you use BifurcationKit.jl for your work, we ask that you **cite** the following paper!! Open source development strongly depends on this. It is referenced on [HAL-Inria](https://hal.archives-ouvertes.fr/hal-02902346) with *bibtex* entry [CITATION.bib](https://github.com/bifurcationkit/BifurcationKit.jl/blob/master/CITATION.bib).