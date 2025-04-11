# [Bifurcation Diagrams](@id bifurcation_diagrams)

Bifurcation diagrams describe how, for a dynamical system, the quantity and type of its steady states change as a parameter is varied[^1]. When using Catalyst-generated models, these can be computed with the [BifurcationKit.jl](https://github.com/bifurcationkit/BifurcationKit.jl) package. Catalyst provides a simple interface for creating BifurcationKit compatible [`BifurcationProblem`](https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/library/#BifurcationKit.BifurcationProblem)s from [`ReactionSystem`](@ref)s.

This tutorial briefly introduces how to use Catalyst with BifurcationKit through basic examples, with BifurcationKit.jl providing [a more extensive documentation](https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/). Especially for more complicated systems, where careful tuning of algorithm options might be required, reading the BifurcationKit documentation is recommended. Finally, BifurcationKit provides many additional features not described here, including [computation of periodic orbits](@ref bifurcationkit_periodic_orbits), [tracking of bifurcation points along secondary parameters](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/branchswitching/), and [bifurcation computations for PDEs](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/tutorials/#PDEs:-bifurcations-of-equilibria).

## [Basic example](@id bifurcation_diagrams_basic_example)
For this example, we will use a modified version of the [Wilhelm model](@ref basic_CRN_library_wilhelm) (which demonstrates a bistable switch as the parameter $k1$ is varied). We declare the model using Catalyst:
```@example ex1
using Catalyst
wilhelm_model = @reaction_network begin
    k1, Y --> 2X
    k2, 2X --> X + Y
    k3, X + Y --> Y
    k4, X --> 0
    k5, 0 --> X
end
```
Next we will create a `BifurcationProblem`. In addition to the `ReactionSystem`, we need to provide:
- The bifurcation parameter (the parameter which is varied in the bifurcation diagram).
- A full model parameter set. This includes the values of all non-bifurcation parameters, but also a value for the bifurcation parameter (which corresponds to the point in parameter space from which the computation of the bifurcation diagram starts).
- An initial guess of the steady state values of the system at the provided parameter set. Using this point as a starting guess for root finding, BifurcationKit calculates an initial steady state from which to compute the bifurcation diagram. Hence, this guess does not need to be very exact (but may be important if the system exhibits multistability for the initial parameter set).
- The species or statistic we wish to plot on the y-axis of the bifurcation diagram (customised properties can be recorded and plotted using the [`record_from_solution` argument](https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/periodicOrbit/#.-record*from*solution)).

We combine all this information to form a `BifurcationProblem`:
```@example ex1
using BifurcationKit
bif_par = :k1
u_guess = [:X => 5.0, :Y => 2.0]
p_start = [:k1 => 8.0, :k2 => 1.0, :k3 => 1.0, :k4 => 1.5, :k5 => 1.25]
plot_var = :X
bprob = BifurcationProblem(wilhelm_model, u_guess, p_start, bif_par; plot_var)
nothing # hide
```

BifurcationKit computes bifurcation diagrams using the `bifurcationdiagram` function. From an initial point in the diagram, it tracks the solution (using a continuation algorithm) until the entire diagram is computed (BifurcationKit's continuation can be used for other purposes, however, this tutorial focuses on bifurcation diagram computation). The continuation settings are provided in a `ContinuationPar` structure. In this example, we will only specify two settings (`p_min` and `p_max`), which sets the minimum and maximum values over which the bifurcation parameter is varied, respectively. We wish to compute a bifurcation diagram over the interval $(2.0,20.0)$, and will use the following settings:
```@example ex1
p_span = (2.0, 20.0)
opts_br = ContinuationPar(p_min = p_span[1], p_max = p_span[2])
nothing # hide
```

Finally, we compute our bifurcation diagram using:
```@example ex1
bif_dia = bifurcationdiagram(bprob, PALC(), 2, opts_br; bothside = true)
nothing # hide
```
Where `PALC()` designates that we wish to use the pseudo-arclength continuation method to track our solution. The third argument (`2`) designates the maximum number of recursions when branches of branches are computed (branches appear as continuation encounters certain bifurcation points). For diagrams with highly branched structures (rare for many common small chemical reaction networks) this input is important. Finally, `bothside = true` designates that we wish to perform continuation on both sides of the initial point (which is typically the case). 

We can plot our bifurcation diagram using the Plots.jl package:
```@example ex1
using Plots
plot(bif_dia; xguide = bif_par, yguide = plot_var, branchlabel = "Steady state concentration")
```
Here, the steady state concentration of $X$ is shown as a function of $k1$'s value. Stable steady states are shown with thick lines, unstable ones with thin lines. The two [fold bifurcation points](https://en.wikipedia.org/wiki/Saddle-node_bifurcation) are marked with "bp".

## [Additional `ContinuationPar` options](@id bifurcation_diagrams_continuationpar)
Most of the options required by the `bifurcationdiagram` function are provided through the `ContinuationPar` structure. For full details, please read the [BifurcationKit documentation](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/library/#BifurcationKit.ContinuationPar). However, a few common options, and how they affect the continuation computation, are described here:
- `p_min` and `p_max`: Set the interval over which the bifurcation diagram is computed (with the continuation stopping if it reaches these bounds).
- `dsmin` and `dsmax`: The minimum and maximum length of the continuation steps (in the bifurcation parameter's value).
- `ds`: The initial length of the continuation steps. This is especially important when `bothside = true` *is not* used, as teh sign of `ds` determines the direction from the initial point in which the continuation will proceed.
- `max_steps`: The maximum number of continuation steps. If a bifurcation diagram looks incomplete, try increasing this value.
- `newton_options`: Options for the Newton's method that BifurcationKit uses to find steady states. This can be created using `NewtonPar(tol = 1e-9, max_iterations = 100)` which here sets the tolerance (to `1e-9`) and the maximum number of newton iterations (to `100`).

The previous bifurcation diagram can be computed, with these various options specified, in the following way:
```@example ex1
p_span = (2.0, 20.0)
newton_options = NewtonPar(tol = 1e-9, max_iterations = 20)
opts_br = ContinuationPar(; p_min = p_span[1], p_max = p_span[2], ds = 0.005, 
                          dsmin = 0.001, dsmax = 0.01, max_steps = 1000, newton_options)
bif_dia = bifurcationdiagram(bprob, PALC(), 2, opts_br; bothside = true)
nothing # hide
```
(In this case, however, these additional settings have no significant effect on the result)

## [Bifurcation diagrams with disjoint branches](@id bifurcation_diagrams_disjoint_branches)
Let's consider the previous case, but instead compute the bifurcation diagram over the interval $(8.0,15.0)$:
```@example ex1
p_span = (8.0, 15.0)
opts_br = ContinuationPar(p_min = p_span[1], p_max = p_span[2])
bif_dia = bifurcationdiagram(bprob, PALC(), 2, opts_br; bothside = true)
plot(bif_dia; ylimit = (0.0, 17.0), xguide = bif_par, yguide = plot_var, branchlabel = "Steady state concentration")
```
Here, in the bistable region, we only see a single branch. The reason is that the continuation algorithm starts at our initial guess (here made at $k1 = 8.0$ for $(X,Y) = (5.0,2.0)$) and tracks the diagram from there. However, with the upper bound set at $k1 = 15.0$ the bifurcation diagram has a disjoint branch structure, preventing the full diagram from being computed by continuation alone. In this case it could be solved by increasing the bound from $k1 = 15.0$, however, this is not possible in all cases. Here we will describe a more general solution, the use of *deflation* to find multiple steady states.

We will not describe the details of how deflation works here (instead, please read [BifurcationKit's documentation](https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/DeflatedContinuation/)). In practise, in involves first creating a deflation continuation algorithm (`defalg`). This one, in turn, requires us to create [a deflation operator](https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/library/#BifurcationKit.DeflationOperator) and a perturbation function.
```@example ex1
deflation_operator = DeflationOperator(2, 0.001, [zeros(numspecies(wilhelm_model))])
perturb_solution(x, _, _) = (x  .+ 0.1 .* rand(length(x)))
defalg = DefCont(; deflation_operator, perturb_solution)
nothing # hide
```
Next we compute our bifurcation diagram using `defalg` (instead of `PALC()` as previous) as our method. We also use a range of additional continuation parameter options to ensure a smooth tracking of the solution. Here we have 
```@example ex1
newton_options = NewtonPar(max_iterations = 10)
cont_par = ContinuationPar(; p_min = p_span[1], p_max = p_span[2], ds = 0.001, max_steps = 10000, newton_options)
bif_dia = continuation(bprob, defalg, cont_par; plot = false, verbosity = 0)
nothing # hide
```
Finally we can plot the diagram, noting that (unlike previously) we are able to compute the full diagram.
```@example ex1
plot(bif_dia; xguide = bif_par, yguide = plot_var, ylimit = (0.0, 17.0), color = 1)
```

## [Systems with conservation laws](@id bifurcation_diagrams_cons_laws)
Some systems are under-determined at steady state, so that for a given parameter set they have an infinite number of possible steady state solutions, preventing bifurcation diagrams from being computed. Similar to when we [compute steady states for fixed parameter values](@ref homotopy_continuation_conservation_laws), we can utilise Catalyst's ability to [detect and eliminate conservation laws](@ref conservation_laws) to resolve this issue. This requires us to provide information of the species concentrations at which we wish to compute the bifurcation diagram (to determine the values of conserved quantities). These are provided to the `BifurcationProblem` using the `u0` argument.

To illustrate this, we will create a simple model of a kinase that is produced and degraded (at rates $p$ and $d$). The kinase facilitates the phosphorylation of a protein ($X$), which is dephosphorylated at a constant rate. For this system, we will compute a bifurcation diagram, showing how the concentration of the phosphorylated protein ($Xp$) depends on the degradation rate of the kinase ($d$). We will set the total amount of protein ($X + Xp$) to $1.0$.
```@example ex2
using BifurcationKit, Catalyst, Plots
kinase_model = @reaction_network begin
    (p, d), 0 <--> K
    (K*kP,kD), X <--> Xp
end

u_guess = [:K => 1.0, :X => 1.0, :Xp => 1.0]
p_start = [:p => 1.0, :d => 0.5, :kP => 2.0, :kD => 5.0]
u0 = [:X => 1.0, :Xp => 0.0]
bprob = BifurcationProblem(kinase_model, u_guess, p_start, :d; plot_var = :Xp, u0)

p_span = (0.1, 10.0)
opts_br = ContinuationPar(p_min = p_span[1], p_max = p_span[2])
bif_dia = bifurcationdiagram(bprob, PALC(), 2, opts_br; bothside = true)
plot(bif_dia; xguide = "d", yguide = "Xp")
```
This bifurcation diagram does not contain any interesting features (such as bifurcation points), and only shows how the steady state concentration of $Xp$ is reduced as $d$ increases. 

Finally, for additional clarity, we reiterate the purpose of the two `u` arguments used:
- `u_guess`: A guess of the initial steady states (which BifurcationKit uses to find its starting point). Typically, most trivial guesses work (e.g. setting all species concentrations to `1.0`). `u_guess` *does not* have to fulfil the conserved concentrations provided in `u0`.
- `u0`: Used to compute the concentrations of any conserved quantities (e.g. in our example $X + Xp = 1.0$). Technically, values are only required for species that are involved in conservation laws (in our case we do not need to provide a value for $K$). However, sometimes determining which species are actually involved in conservation laws can be difficult, and it might be easier to simply provide concentrations for all species.


---
## [Citation](@id bifurcationkit_periodic_orbits_citation)
If you use BifurcationKit.jl for your work, we ask that you **cite** the following paper!! Open source development strongly depends on this. It is referenced on [HAL-Inria](https://hal.archives-ouvertes.fr/hal-02902346) with *bibtex* entry [CITATION.bib](https://github.com/bifurcationkit/BifurcationKit.jl/blob/master/CITATION.bib).

---
## References
[^1]: [Yuri A. Kuznetsov, *Elements of Applied Bifurcation Theory*, Springer (2023).](https://link.springer.com/book/10.1007/978-3-031-22007-4)