# # [Tracking Bifurcation Point w.r.t. Secondary Parameters using BifurcationKit](@id bifurcationkit_bifpoint_continuation)
Previously, we have shown how to [compute bifurcation diagrams](@ref bifurcation_diagrams) using [BifurcationKit.jl](https://github.com/bifurcationkit/BifurcationKit.jl). In this example, we will show how, after computing the initial diagram, we can track how the position of a bifurcation point moves as a secondary parameter is changed (so-called codimensional 2 bifurcation analysis). More information on how to track bifurcation points along secondary parameters can be found in the [BifurcationKit documentation](https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/tutorials/ode/tutorialCO/#CO-oxidation-(codim-2)).

## [Computing the bifurcation diagram for the Repressilator](@id bifurcationkit_bifpoint_continuation_bifdia)
We will first compute the bifurcation diagram, using the same approach as in the [corresponding tutorial](@ref bifurcation_diagrams). For this example, we will use the oscillating [Repressilator](@ref basic_CRN_library_repressilator) model.
```@example bifurcationkit_bifpoint_continuation
using Catalyst
repressilator = @reaction_network begin
    hillr(Z,v,K,n), ∅ --> X
    hillr(X,v,K,n), ∅ --> Y
    hillr(Y,v,K,n), ∅ --> Z
    d, (X, Y, Z) --> ∅
end
```
Next, we create a `BifurcationProblem` for our model. We will compute the bifurcation diagram with respect to the parameter $v$, and plot the species $X$ in the diagram.
```@example bifurcationkit_bifpoint_continuation
using BifurcationKit
bif_par = :v
u_guess = [:X => 15.0, :Y => 15.0, :Z => 15.0]
p_start = [:v => 10.0, :K => 15.0, :n => 3, :d => 0.2]
plot_var = :X
bprob = BifurcationProblem(repressilator, u_guess, p_start, bif_par; plot_var)
nothing # hide
```
We compute the bifurcation diagram using the `bifurcationdiagram` function. We will compute it across the interval $v \in (0,20)$.
```@example bifurcationkit_bifpoint_continuation
v_span = (0.0, 20.0)
opts_br = ContinuationPar(p_min = v_span[1], p_max = v_span[2])
bifdia = bifurcationdiagram(bprob, PALC(), 3, opts_br; bothside = true)
nothing # hide
```
Finally, we plot the bifurcation diagram.
```@example bifurcationkit_bifpoint_continuation
using Plots
plot(bifdia; xguide = bif_par, yguide = plot_var, branchlabel = "Continuation of steady state w.r.t. v",
    linewidthstable = 6, linewidthunstable = 3, markersize = 5)
```
We note that for small values of $v$ the system's single steady state is stable (where the line is thicker). After a [Hopf](https://en.wikipedia.org/wiki/Hopf_bifurcation) bifurcation (the red point), the state turns unstable (where the line is thinner). For chemical reaction networks (which mostly are well-behaved) a single unstable steady state typically corresponds to an oscillation. We can confirm that the system oscillates in the unstable region (while it reaches a stable steady state in the stable region) using simulations:
```@example bifurcationkit_bifpoint_continuation
using OrdinaryDiffEqDefault
p_nosc = [:v => 5.0, :K => 15.0, :n => 3, :d => 0.2]
p_osc = [:v => 15.0, :K => 15.0, :n => 3, :d => 0.2]
prob_nosc = ODEProblem(repressilator, u_guess, 80.0, p_nosc)
prob_osc = ODEProblem(repressilator, u_guess, 80.0, p_osc)
sol_nosc = OrdinaryDiffEqDefault.solve(prob_nosc)
sol_osc = OrdinaryDiffEqDefault.solve(prob_osc)
plot(plot(sol_nosc; title = "v = 5"), plot(sol_osc; title = "v = 15"), size = (1000,400), lw = 4)
```

## [Tracking the bifurcation point w.r.t. a second parameter](@id bifurcationkit_bifpoint_continuation_codim2)
Next, we will investigate how the Hopf bifurcation point moves (in $v$-$X$ space) as a second parameter ($K$) is changed. To do this we will use BifurcationKit.jl's [`continuation` function](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/library/#BifurcationKit.continuation) (the [`bifurcationdiagram` function](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/library/#BifurcationKit.bifurcationdiagram), which we previously have used, works by calling `continuation` recursively). We will call it on the Hopf bifurcation point. First, we will designate the parameter we wish to change, as well as which interval to change it. For this, we compute the parameter index (`K_idx`) that BifurcationKit will see (a smoother interface for this will hopefully be added in the future). We also create a new `ContinuationPar`, to which we add our $K$ interval.
```@example bifurcationkit_bifpoint_continuation
K_span = (0.01, 27.0)
opts_br_2 = ContinuationPar(p_min = K_span[1], p_max = K_span[2])
K_idx = findfirst(isequal(repressilator.K), parameters(complete(convert(NonlinearSystem, repressilator))))
nothing # hide
```
Now we can compute the continuation of the Hopf bifurcation. First we must extract a branch from our bifurcation diagram (using [`get_branch`](https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/library/#BifurcationKit.get_branch)), as `continuation` cannot work on bifurcation diagrams directly (our bifurcation diagram consists of a single branch, which we here extract). Next, we can call `continuation` on it, designating that we wish to perform continuation from the second point on the branch (which corresponds to the Hopf bifurcation, the first point is one of the branch's two endpoints).
```@example bifurcationkit_bifpoint_continuation
branch = get_branch(bifdia, ()).γ
cont_hopf = continuation(branch, 2, (@optic _[K_idx]), opts_br_2; start_with_eigen = true, bothside = true)
nothing # hide
```
We can now plot how the position of the bifurcation point changes with $K$. Here, we use `vars = (:p1, :x)` to designate that we wish to plot (across the continuation branch) the plotting variable ($X$, which we designated when we created our `BifurcationProblem`) against the first parameter ($v$).
```@example bifurcationkit_bifpoint_continuation
plot(bifdia; branchlabel = "Continuation of steady state w.r.t. v, (K = 15)", linewidthstable = 6, linewidthunstable = 3, markersize = 5)
plot!(cont_hopf; vars = (:p1, :x), xlimit = v_span, xguide = bif_par, yguide = plot_var, branchlabel = "Continuation of Hopf bifurcation w.r.t. K")
```
In this case we cannot see directly which part of the $K$ continuation branch corresponds to low values, however, for low $K$ the Hopf bifurcation occurs for much lower values of $v$ (and corresponds to lower steady state values of $X$). We can check this by e.g. re-computing the Hopf branch for `K_span = (0.01, 20.0)` and see that the rightmost part of the branch is shortened. 

We can confirm that the new line corresponds to the Hopf Bifurcation point by recomputing the initial bifurcation diagram, but for a lower $K$ value.
```@example bifurcationkit_bifpoint_continuation
p_start_2 = [:v => 10.0, :K => 10.0, :n => 3, :d => 0.2]
bprob_2 = BifurcationProblem(repressilator, u_guess, p_start_2, bif_par; plot_var)
bifdia_2 = bifurcationdiagram(bprob_2, PALC(), 3, opts_br; bothside = true)
plot!(bifdia_2; xguide = bif_par, yguide = plot_var, branchlabel = "Continuation of steady state w.r.t. v, (K = 10)",
    linewidthstable = 6, linewidthunstable = 3, markersize = 5)
```
Here we see that the Hopf bifurcation point of the new diagram also lies on the Hopf continuation line.

Finally, we have already noted that the Hopf bifurcation splits parameter space into one part where the system oscillates and one where it doesn't. Previously we plotted the Hopf continuation in $v$-$X$ space, however, it is also possible to plot it in $v$-$K$ space using the `vars = (:p1, :p2)` argument:
```@example bifurcationkit_bifpoint_continuation
xlimit = extrema(getfield.(cont_hopf.branch, :p1))
ylimit = extrema(getfield.(cont_hopf.branch, :p2))
plot(cont_hopf; vars = (:p1, :p2), xlimit, ylimit, branchlabel = "Hopf bifurcation", xguide = "v", yguide = "K", lw = 6)
```
Next, we colour parameter space according to whether the steady state is stable (blue) or unstable (red). We also mark two sample values (one in each region).
```@example bifurcationkit_bifpoint_continuation
sample1 = (15.0, 10.0)
sample2 = (5.0, 15.0)
plot(cont_hopf; vars = (:p1, :p2), fillrange = ylimit[2])
plot!(cont_hopf; vars = (:p1, :p2), fillrange = ylimit[1], xlimit, ylimit)
scatter!(sample1; label = "Oscillatory parameter set", markersize = 7)
scatter!(sample2; label = "Non-oscillatory parameter set", markersize = 7)
```
Finally, we can perform one simulation using each of the parameter samples, confirming that one corresponds to an oscillation, while the other one does not.
```@example bifurcationkit_bifpoint_continuation
ps_osc = [:v => sample1[1], :K => sample1[2], :n => 3, :d => 0.2]
ps_nosc = [:v => sample2[1], :K => sample2[2], :n => 3, :d => 0.2]
oprob_osc = ODEProblem(repressilator, u_guess, 100.0, ps_osc)
oprob_nosc = ODEProblem(repressilator, u_guess, 100.0, ps_nosc)
sol_osc = OrdinaryDiffEqDefault.solve(oprob_osc)
sol_nosc = OrdinaryDiffEqDefault.solve(oprob_nosc)
plot(plot(sol_osc; title = "Oscillation"), plot(sol_nosc; title = "No oscillation"); size = (1000, 400), lw = 4)
```