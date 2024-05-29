# [Analysing model steady state properties with DynamicalSystems.jl](@id dynamical_systems)
The [DynamicalSystems.jl package](https://github.com/JuliaDynamics/DynamicalSystems.jl) implements a wide range of methods for analysing dynamical systems. This includes both continuous-time systems (i.e. ODEs) and discrete-times ones (difference equations, however, these are not relevant to chemical reaction network modelling). Here we give two examples of how DynamicalSystems.jl can be used, with the package's [documentation describing many more features](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/dynamicalsystems/dev/tutorial/). Finally, it should also be noted that DynamicalSystems.jl contain several tools for [analysing data measured from dynamical systems](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/dynamicalsystems/dev/contents/#Exported-submodules).

## [Finding basins of attraction](@id dynamical_systems_basins_of_attraction)
Given enough time, an ODE will eventually reach a so-called [*attractor*](https://en.wikipedia.org/wiki/Attractor). For chemical reaction networks (CRNs), this will typically be either a *steady state* or a *limit cycle*. Since ODEs are deterministic, which attractor a simulation will reach is uniquely determined by the initial condition (assuming parameter values are fixed). Conversely, each attractor is associated with a set of initial conditions such that model simulations originating in these will tend to that attractor. These sets are called *basins of attraction*. Here, phase space (the space of all possible states of the system) can be divided into a number of basins of attraction equal to the number of attractors. 

DynamicalSystems.jl provides a simple interface for finding an ODE's basins of attraction across any given subspace of phase space. In this example we will use the bistable [Wilhelm model](https://bmcsystbiol.biomedcentral.com/articles/10.1186/1752-0509-3-90) (which steady states we have previous [computed using homotopy continuation](@ref homotopy_continuation_basic_example)). As a first step, we create an `ODEProblem` corresponding to the model which basins of attraction we wish to compute. For this application, `u0` and `tspan` is unused, and their values are of little importance (the only exception is than `tspan`, for implementation reason, must provide a not too small interval, we recommend minimum `(0.0, 1.0)`). 
```@example dynamical_systems_basins
using Catalyst
wilhelm_model = @reaction_network begin
    k1, Y --> 2X
    k2, 2X --> X + Y
    k3, X + Y --> Y
    k4, X --> 0
end

u0 = [:X => 0.0, :Y => .0]
tspan = (0.0, 10.0)
ps = [:k1 => 8.0, :k2 => 2.0, :k3 => 1.0, :k4 => 1.5]
oprob = ODEProblem(wilhelm_model, u0, tspan, ps)
nothing # hide
```
Next, for any application of DynamicalSystems.jl, our `ODEProblem` must be converted into a so-called `CoupledODEs` structure. This is done by combining the ODE with the solver (and potential solver options) with which we wish to simulate it (just like when it is simulated using `solve`). Here, we will simply designate the `Tsit5` numeric solver (but provide no other options).
```@example dynamical_systems_basins
using DynamicalSystems, OrdinaryDiffEq
ds = CoupledODEs(oprob, (alg = Tsit5(),))
```
We can now compute the basins of attraction. This is done by first creating a grid that designates which subspace of phase-space we wish to investigate (here, the corresponding basin of attraction is found for every point on the grid). Next, we create a `AttractorsViaRecurrences` struct, that maps initial conditions to attractors, and then use that as input to the `basins_of_attraction` function.
```@example dynamical_systems_basins
# We provide one grid of values for each species. These are then bundled into a tuple.
x_grid = 0.0:0.03:6.0
y_grid = 0.0:0.03:9.0
grid = (x_grid, y_grid)
avr = AttractorsViaRecurrences(ds, grid)
basins, attractors = basins_of_attraction(avr, grid; show_progress = false)
attractors
```
Here, `attractors` is a dictionary that maps attractor labels (integers) to attractors. In this case we have two fixed points, one at $(0.0,0.0)$ and one at $(4.5,6.0)$. Next, `basins` is a matrix of equal size to `grid`, where each value is an integer describing to which attractor's basin that state belongs.

DynamicalSystems.jl also provides a simple interface for plotting the resulting basins. This uses [Makie.jl](https://docs.makie.org/stable/), an alternative plotting package to [Plots.jl](https://github.com/JuliaPlots/Plots.jl) (which is typically the preferred plotting package within the context of Catalyst). Generally, Makie is good at creating animations or interactive graphics (however, it is also a [popular competitor to Plots.jl for general-purpose plotting](https://juliapackagecomparisons.github.io/pages/plotting/)). 
```@example dynamical_systems_basins
using CairoMakie
heatmap_basins_attractors(grid, basins, attractors)
```
Here, in addition to the basins of attraction, the system's three steady states are marked (the one at the intersection of the two basins is unstable).

!!! warning
    Both Makie and Plots.jl exports a function called `plot`. Hence, if both these packages are imported into the same session, calls to `plot` must be prepended with the package one wishes to use (e.g. `Plots.plot(sol)`).

More information on how to compute basins of attractions for ODEs using DynamicalSystems.jl can be found [here](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/attractors/stable/basins/).

## [Computing Lyapunov exponents](@id dynamical_systems_lyapunov_exponents)
[*Lyapunov exponents*](https://en.wikipedia.org/wiki/Lyapunov_exponent) are scalar values that can be computed for any attractor of an ODE. For an ODE with $n$ variables, for each state, a total of $n$ Lyapunov exponents can be computed (and these are collectively called the *Lyapunov spectrum*). Positive Lyapunov exponents indicate that trajectories initially infinitesimally close diverge from each other. Conversely, negative Lyapunov exponents suggests that trajectories converge to each others. 

While Lyapunov exponents can be used for other purposes, they are primarily used to characterise [*chaotic behaviours*](https://en.wikipedia.org/wiki/Chaos_theory) (where small changes in initial conditions has large effect on the resulting trajectories). Generally, an ODE exhibit chaotic behaviour if its attractor(s) have *at least one* positive Lyapunov exponent. Practically, Lyapunov exponents can be computed using DynamicalSystems.jl's `lyapunovspectrum` function. Here we will use it to investigate two models, one which exhibits chaos and one which do not.

First, let us consider the [Willamowski–Rössler model](@ref basic_CRN_library_wr), which is known to exhibit chaotic behaviour.
```@example dynamical_systems_lyapunov
using Catalyst
wr_model = @reaction_network begin
    k1, 2X --> 3X
    k2, X --> 2X
    k3, Z + 2X --> 2Z
    k4, Y + X --> 2Y
    k5, Y --> ∅
    k6, 2Z --> ∅
    k7, Z --> ∅
end
```
We can simulate the model, noting that its behaviour seems chaotic.
```@example dynamical_systems_lyapunov
using OrdinaryDiffEq, Plots

u0 = [:X => 1.5, :Y => 1.5, :Z => 1.5]
tspan = (0.0, 100.0)
p = [:k1 => 2.1, :k2 => 0.7, :k3 => 2.9, :k4 => 1.1, :k5 => 1.0, :k6 => 0.5, :k7 => 2.7]

oprob = ODEProblem(wr_model, u0, tspan, p)
sol = solve(oprob, Rodas5P())
plot(sol; idxs=(:X, :Y, :Z))
```
Next, like when we [computed basins of attraction](@ref dynamical_systems_basins_of_attraction), we create a `CoupledODEs` corresponding to the model and state for which we wish to compute our Lyapunov spectrum. Lke previously, `tspan` must provide some small interval (at least `(0.0, 1.0)` is recommended), but else have no impact on the computed Lyapunov spectrum.
```@example dynamical_systems_lyapunov
using DynamicalSystems
ds = CoupledODEs(oprob, (alg = Rodas5P(autodiff = false),))
nothing # hide
```
Here, the `autodiff = false` argument is required when Lyapunov spectrums are computed. We can now provide our `CoupledODEs` (`ds`) to `lyapunovspectrum` to compute the lyapunov spectrum. This function requires a second argument (here set to `100`). Generally setting this to a higher value will increase accuracy, but also increase runtime (since `lyapunovspectrum` is fast for most systems, setting this to a large value is recommended).
```@example dynamical_systems_lyapunov
lyapunovspectrum(ds, 100)
```
Here, the largest exponent is positive, suggesting that the model is chaotic (or, more accurately, it has at least one chaotic attractor, to which is approached from the initial condition $(1.5,1.5,1.5)$).

Next, we consider the [Brusselator] model. First we simulate the model for two similar initial conditions, confirming that they converge to the same limit cycle:
```@example dynamical_systems_lyapunov
brusselator = @reaction_network begin
    A, ∅ --> X
    1, 2X + Y --> 3X
    B, X --> Y
    1, X --> ∅
end

u0_1 = [:X => 1.0, :Y => 1.0]
u0_2 = [:X => 1.2, :Y => 1.0]
tspan = (0., 25.)
ps = [:A => 1.0, :B => 4.0]

oprob1 = ODEProblem(brusselator, u0_1, tspan, ps)
oprob2 = ODEProblem(brusselator, u0_2, tspan, ps)
osol1  = solve(oprob1, Tsit5())
osol2  = solve(oprob2, Tsit5())
plot(osol1; idxs = (:X, :Y))
plot!(osol2; idxs = (:X, :Y))
```
Next, we compute the Lyapunov spectrum at one of the initial conditions:
```@example dynamical_systems_lyapunov
ds = CoupledODEs(oprob1, (alg = Rodas5P(autodiff = false),))
lyapunovspectrum(ds, 100)
```
Here, all Lyapunov exponents are negative, confirming that the brusselator is non-chaotic.

More details on how to compute Lyapunov exponents using DynamicalSystems.jl can be found [here](https://juliadynamics.github.io/ChaosTools.jl/stable/lyapunovs/). A full overview of tools for analysing chaotic behaviours (using the "ChaosTools.jl subpackage) can be found [here](https://juliadynamics.github.io/ChaosTools.jl/stable/).


---
## [Citations](@id dynamical_systems_citations)
If you use this functionality in your research, [in addition to Catalyst](@ref catalyst_citation), please cite the following paper to support the author of the DynamicalSystems.jl package:
```
@article{DynamicalSystems.jl-2018,
  doi = {10.21105/joss.00598},
  url = {https://doi.org/10.21105/joss.00598},
  year  = {2018},
  month = {mar},
  volume = {3},
  number = {23},
  pages = {598},
  author = {George Datseris},
  title = {DynamicalSystems.jl: A Julia software library for chaos and nonlinear dynamics},
  journal = {Journal of Open Source Software}
}
```


---
## Learning more

If you want to learn more about analysing dynamical systems, including chaotic behaviour, you can have a look at the textbook [Nonlinear Dynamics](https://link.springer.com/book/10.1007/978-3-030-91032-7). It utilizes DynamicalSystems.jl and provides a concise, hands-on approach to learning nonlinear dynamics and analysing dynamical systems [^1].


---
## References
[^1]: [G. Datseris, U. Parlitz, *Nonlinear dynamics: A concise introduction interlaced with code*, Springer (2022).](https://link.springer.com/book/10.1007/978-3-030-91032-7)
[^2]: [S. H. Strogatz, *Nonlinear Dynamics and Chaos*, Westview Press (1994).](http://users.uoa.gr/~pjioannou/nonlin/Strogatz,%20S.%20H.%20-%20Nonlinear%20Dynamics%20And%20Chaos.pdf)
[^3]: [A. M. Lyapunov, *The general problem of the stability of motion*, International Journal of Control (1992).](https://www.tandfonline.com/doi/abs/10.1080/00207179208934253)