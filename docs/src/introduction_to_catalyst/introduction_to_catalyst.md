# [Introduction to Catalyst](@id introduction_to_catalyst)
This tutorial provides a basic introduction on how to create chemical reaction network (CRN) models using Catalyst, and how to perform ODE, SDE, and jump simulation using these. [An alternative introduction](@ref catalyst_for_new_julia_users) is available for users with little familiarity with Julia programming. 

Before starting this tutorial, we must (unless this has already been done) [install and activate](@ref ref) the Catalyst package:
```julia
using Pkg
Pkg.add("Catalyst")
using Catalyst
```

## [Creating basic Catalyst models](@id introduction_to_catalyst_model_creation)
Catalyst exports the [`@reaction_network`](@ref) [macro]([macro](https://docs.julialang.org/en/v1/manual/metaprogramming/#man-macros)), which provides the main way of creating CRN models (alternative approaches are described [later](@ref ref)). This is followed by a `begin ... end` block which encapsulates the reactions of the system. E.g here we create a simple [birth-death model](@ref ref) (where a single species $X$ is both produced and degraded):
```@example intro_1
using Catalyst # hide
bd_model = @reaction_network begin
 p, 0 --> X
 d, X --> 0
end
```
Next, we create a simple [Michaelis-Menten enzyme kinetics model](@ref ref) (where an enzyme, $E$ converts a substrate, $S$, to a product, $P$):
```@example intro_1
mm_model = @reaction_network begin
 kB, S + E --> SE
 kD, SE --> S + E
 kP, SE --> P + E
end
```
And finally, a [SIR model of an infectious disease](@ref ref) (where susceptible individuals, $I$, get infected by infected individuals, $I$, and later recovers, $R$):
```@example intro_1
sir_model = @reaction_network begin
    β, S + I --> 2I
    γ, I --> R
end
```

In all three cases, each reaction consists of:
- A rate (at which the reaction occurs).
- A set of substrates (species that are consumed by the reaction).
- A set of products (species that are produced by the reaction).

For reactions with multiple substrates/products (e.g. `S + E --> SE`), these are separated by `+`. Where substrates/products contain multiple units of the same species (e.g. `S + I --> 2I`), this is indicated by pre-appending that species with the corresponding [stoichiometric coefficient](https://en.wikipedia.org/wiki/Stoichiometry). The absences of substrates/products (that is, in pure production/degradation reactions like `0 --> X` and `X --> 0`) are denoted with `0`. A more throughout description of how to create models using `@reaction_network` (also called *the Catalyst DSL*) is provided [here](@ref ref).

## [Simulating Catalyst models](@id introduction_to_catalyst_model_simulation)
There exist three primary modes of simulation CRN models:
1. Using ordinary differential equations (ODEs).
2. Using stochastic differential equations (SDEs).
3. Using jump process simulations.

Catalyst models can be simulated using all three modes. Below we give a simple introduction to each, with [more detailed documentation provided elsewhere](@ref ref).

### [ODE simulations](@id introduction_to_catalyst_model_simulation_ode)
Let us consider the infectious disease model declared in the previous section. To simulate it we need (in addition to the model):
* The simulation's initial condition. That is, the concentration (or copy numbers) of each species at the start of the simulation.
* The simulation's timespan. That is, the time frame over which we wish to run the simulation.
* The simulation's parameter values. That is, the values of the model's parameters for this simulation.

The initial conditions and parameter values are declared as vectors of pairs. Here, each entry pairs the symbol corresponding to a specific quantity's name to its value (e.g. `:S => 99` denotes that the initial concentration of `S` should be `99`). The timespan is simply a tuple with the simulations *starting* and *final* times. Once we have set these values, we collect them in a so-called `ODEProblem`.
```@example intro_1
u0 = [:S => 99, :I => 1, :R => 0]
tspan = (0.0, 20.0)
ps = [:β => 0.025, :γ => 0.2]
oprob = ODEProblem(sir_model, u0, tspan, ps)
```
To simulate our `ODEProblem`, we simply input it to the `solve` command. Before we can do so, however, we also need to activate the [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) package (which is used for all ODE simulations in Julia).
```@example intro_1
using OrdinaryDiffEq
sol = solve(oprob)
```
Finally, we can plot the solution using the [Plots.jl](https://github.com/JuliaPlots/Plots.jl) package's `plot` function:
```@example intro_1
using Plots
plot(sol)
```

### [SDE simulations](@id introduction_to_catalyst_model_simulation_sde)
SDE simulations can be carried out in a very similar manner to ODE ones, but by creating an `SDEProblem` instead of an `ODEProblem`. However, we can re-use the problem inputs from previously:
```@example intro_1
sprob = SDEProblem(sir_model, u0, tspan, ps)
```
Next, we import the [StochasticDiffEq.jl](https://github.com/SciML/StochasticDiffEq.jl) package (which is used for all SDE simulations in Julia) and use it to simulate our SDE. 
```@example intro_1
using StochasticDiffEq
sol = solve(sprob, STrapezoid())
sol = solve(sprob, STrapezoid(); seed = 1234) # hide
nothing # hide
```
Note that we here have to provide a second argument to the `solver` command. This is our choice of [simulation method](@ref ref). While we can provide this for ODE simulations as well, OrdinaryDiffEq.jl is able to [automatically select a suitable solver for any problem](@ref ref). This is currently not possible for SDEs. For now, `STrapezoid` is often a passable default choice, however, for important applications it can be good to [closer study the available SDE solvers](@ref ref).

Next, we can again use `plot` to plot the solution.
```@example intro_1
plot(sol)
```
Since SDEs are *stochastic* (random), if we simulate it again we get another result:
```@example intro_1
sol = solve(sprob, ImplicitEM())
sol = solve(sprob, STrapezoid(); seed = 5678) # hide
plot(sol)
```

### [Jump simulations](@id introduction_to_catalyst_model_simulation_jumps)
Finally, Catalyst models can be simulated using jump simulations. These simulate the exact (and randomised) occurrence of system reactions, and their effect on the system's state. Jump simulations require the creation of `JumpProblem`. These are a bit different from `ODEProblem`s and `SDEProblem`s in that we first create a `DiscreteProblem` (to which we provide the model, initial condition, time span, and parameter values), and then use this (in addition to again providing our model) as input to `JumpProblem`.
```@example intro_1
dprob = DiscreteProblem(sir_model, u0, tspan, ps)
jprob = JumpProblem(sir_model, dprob)
nothing # hide
```
Finally, this problem can be simulated and plotted just like ODEs and SDEs (however, this time requiring the [JumpProcesses.jl](https://github.com/SciML/JumpProcesses.jl) simulation package):
```@example intro_1
using JumpProcesses
sol = solve(jprob)
sol = solve(jprob; seed = 12345) # hide
plot(sol)
```
Looking at the plot, we can actually distinguish the individual jumps of the simulation.

### [Basic model and simulation querying](@id introduction_to_catalyst_querying)
In this section, we will consider some basic functions for querying the content of Catalyst models. If we wish to list the species a model consists of, we can use the `species` function:
```@example intro_1
species(sir_model)
```
Next, to list the parameters we can use the `parameters` function:
```@example intro_1
parameters(sir_model)
```
!!! note
 The `species` and `parameters` functions return a model's species/parameters as so-called *symbolic variables*. This is a special type (which can be used to form *symbolic expressions) that is described in more detail [here](@ref ref).

Finally, to list the reactions we can use
```@example intro_1
reactions(sir_model)
```
A full list of similar and related functions can be found in [the API](@ref ref).

Simulation solutions can also be queried in various ways. To receive the value of species $P$ across the solution `sol` we simply call
```@example intro_1
sol[:P]
```
A more throughout tutorial on how to query solutions (and other relevant structures) for parameters and species values can be found [here](@ref simulation_structure_interfacing).

Finally, it is possible to print a model in [LaTeX format](https://en.wikipedia.org/wiki/LaTeX) using the [Latexify.jl](https://github.com/korsbo/Latexify.jl) package. To print it formatted as chemical reactions we simply call `latexify` with our model as input:
```@example intro_1
using Latexify
latexify(sir_model)
```
If we instead wish to print equations the model is converted to for ODE simulations, we have to add the `form = :ode` argument:
```@example intro_1
latexify(sir_model; form = :ode)
```
Printing of models using Latexify is described in more detail [here](@ref ref).

!!! note
    To have `latexify`'s output printed in LaTeX format you must use a programming environment that actually supports this (like a [JuPyteR](https://github.com/JuliaLang/IJulia.jl) or [Pluto](https://github.com/fonsp/Pluto.jl) notebook). Otherwise, you will simply print the LaTeX code which would generate the LaTeX print. 

## [Reaction rate laws used in ODE simulations](@id introduction_to_catalyst_rate_equations)
When generating mathematical models from Catalyst-generated [`ReactionSystem`](@ref) models, reaction rates are treated as *microscopic* rates. That is, for a general mass action reaction of the form 
```math
n_1 S_1 + n_2 S_2 + \dots n_M S_M \to \dots
```
 with stoichiometric substrate coefficients $\{n_i\}_{i=1}^M$ and rate $k$, the corresponding ODE and SDE rate laws are taken to be
```math
k \prod_{i=1}^M \frac{(S_i)^{n_i}}{n_i!},
```
while the jump process transition rate (i.e., the propensity or intensity function) is
```math
k \prod_{i=1}^M \frac{S_i (S_i-1) \dots (S_i-n_i+1)}{n_i!}.
```
For example, the rate law of the reaction $2X + 3Y \to Z$ with rate constant $k$ would be
```math
k \frac{X^2}{2!} \frac{Y^3}{3!} \\
```
giving the ODE model
```math
\begin{align*}
\frac{dX}{dt} &= -2 k \frac{X^2}{2!} \frac{Y^3}{3!}, &
\frac{dY}{dt} &= -3 k \frac{X^2}{2!} \frac{Y^3}{3!}, &
\frac{dZ}{dt} &= k \frac{X^2}{2!} \frac{Y^3}{3!}.
\end{align*}
```
This implicit rescaling of rate constants can be disabled through [explicit providing the `combinatoric_ratelaws=false` argument](@ref ref) when a set of equations is generated.

For the previous example using this keyword argument would give the rate law
```math
k X^2 Y^3
```
and the ODE model
```math
\begin{align*}
\frac{dX}{dt} &= -2 k X^2 Y^3, &
\frac{dY}{dt} &= -3 k X^2 Y^3, &
\frac{dZ}{dt} &= k X^2 Y^3.
\end{align*}
```

## [Next steps](@id introduction_to_catalyst_next_steps)
The above tutorial gives enough background to use create and simulate basic CRN models using Catalyst. The remaining documentation goes through additional features of the package. In some places, it will also provide introductions to basic theory, and general advice on how to carry out various workflows. While it would be possible to read it from start to finish, you might also select specific sections depending on your intended use of Catalyst:
- If you have not read it already, this documentation's [home](@ref ref) page provides some useful information.
- The [basic](@ref ref) and [advanced](@ref ref) modelling tutorials describe a range of options for model creation, and are useful to everyone who plans on using Catalyst extensively.
- The introduction to [simulations](@ref ref), [spatial modelling](@ref ref), and [inverse problem solving](@ref ref) are useful to anyone who plans on using Catalyst for any of these purposes.
- The [Advanced introduction to Catalyst](@ref ref) is not required to utilise any Catalyst feature. While it is not intended to be read at an early stage, it is still useful to anyone who has, or plans to, use Catalyst extensively.


---
## [References](@id introduction_to_catalyst_references)
[^1]: [Torkel E. Loman, Yingbo Ma, Vasily Ilin, Shashi Gowda, Niklas Korsbo, Nikhil Yewale, Chris Rackauckas, Samuel A. Isaacson, *Catalyst: Fast and flexible modeling of reaction networks*, PLOS Computational Biology (2023).](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011530)