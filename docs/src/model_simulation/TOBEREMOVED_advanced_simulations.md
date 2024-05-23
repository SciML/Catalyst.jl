# [Advanced Simulation Options](@id advanced_simulations)
Throughout the preceding tutorials, we have shown the basics of how to solve
ODE, SDE, and jump process models generated from Catalyst `ReactionSystem`s. In
this tutorial we'll illustrate some more advanced functionality that can be
useful in many modeling contexts, and that provide conveniences for common
workflows. For a comprehensive overview of solver properties, parameters, and
manipulating solution objects, please read the [documentation of the
DifferentialEquations package](https://docs.sciml.ai/DiffEqDocs/stable/), which
Catalyst uses for all simulations.


## [Monte Carlo simulations using `EnsembleProblem`s](@id advanced_simulations_ensemble_problems)
In many contexts one needs to run multiple simulations of a model, for example
to collect statistics of SDE or jump process solutions, or to systematically
vary parameter values within a model. While it is always possible to manually
run such ensembles of simulations via a `for` loop, DifferentialEquations.jl
provides the `EnsembleProblem` as a convenience to manage structured collections
of simulations. `EnsembleProblem`s provide a simple interface for modifying a
problem between individual simulations, and offers several options for batching
and/or parallelizing simulation runs. For a more thorough description, please
read [the Parallel Ensemble Simulations section of the DifferentialEquations
documentation](https://docs.sciml.ai/DiffEqDocs/stable/features/ensemble/#ensemble).
Here, we will give a brief introduction to the use of `EnsembleProblem`s from
Catalyst-generated models.

Let's look at a single-component bistable self-activation model:
```@example ex1
using Catalyst, DifferentialEquations, Plots

rn = @reaction_network begin
    v0 + hill(X,v,K,n), ∅ --> X
    deg, X --> ∅
end
u0 = [:X => 0.0]
tspan = (0.0,1000.0)
p = [:v0 => 0.1, :v => 2.5, :K => 75.0, :n => 2.0, :deg => 0.01];
sprob = SDEProblem(rn, u0, tspan, p)
nothing # hide
```
we can then use our `SDEProblem` as input to an `EnsembleProblem`:
```@example ex1
eprob = EnsembleProblem(sprob)
```
The `EnsembleProblem` can now be used as input to the `solve` command. It has
the same options as when simulating the `SDEProblem` directly, however, it has
an additional argument `trajectories` to determine how many simulations should
be performed.
```@example ex1
esol = solve(eprob; trajectories=5)
```
This simulation is automatically multithreaded over all available threads.
Please read [this
documentation](https://docs.sciml.ai/DiffEqDocs/stable/features/ensemble/#EnsembleAlgorithms)
for more information on parallelisation alternatives. The ensemble simulations
can be plotted using the `plot` function, which by default displays all
trajectories:
```@example ex1
plot(esol)
```

Sometimes when performing a large number of ensemble simulations, the plots get
very dense. In these cases, the plot argument `linealpha` (which sets trajectory
transparency) may be useful:
```@example ex1
esol = solve(eprob; trajectories = 100)
plot(esol; linealpha = 0.5)
```

Sometimes, one wishes to perform the same simulation a large number of times,
while making minor modifications to the problem each time. This can be done by
giving a problem function, `prob_func`, argument to the `EnsembleProblem`. Let
us consider ODE simulations of a simple birth/death process:
```@example ex1
rn = @reaction_network begin
    (b,1.0), ∅ <--> X
end
u0 = [:X => 1.0]
tspan = (0.0, 1.0)
p = [:b => 1.];
oprob = ODEProblem(rn, u0, tspan, p)
nothing # hide
```
We wish to simulate this model for a large number of values of `b`. We do this
by creating a `prob_func` that will make a modification to the problem at the
start of each Monte Carlo simulation:
```@example ex1
b_values = 1.0:0.1:2.0
function prob_func(prob, i, repeat)
    @unpack b = prob.f.sys    # Fetches the b parameter to be used in the local scope.
    remake(prob; p = [b => b_values[i]])
end
nothing # hide
```
Here, `prob_func` takes three arguments:
 - `prob`: The problem given to our `EnsembleProblem`, this is the problem that
   `prob_func` modifies in each iteration.
 - `i`: The number of this specific Monte Carlo iteration in the interval
   `1:trajectories`.
 - `repeat`: The repeat of this specific Monte Carlo simulation (We will ignore
this argument in this brief overview). In our case, for each Monte Carlo
simulation, our `prob_func` takes our original `ODEProblem` and uses the
`remake` function to change the parameter vector. Here, for the `i`th Monte
Carlo simulation, the value of `b` is also the `i`th value of our `b_values`
vector. Finally, we can simulate and plot our problem:
```@example ex1
eprob = EnsembleProblem(oprob; prob_func = prob_func)
esol = solve(eprob; trajectories = length(b_values))
plot(esol)
```

Note that plot legends are disabled when plotting ensemble solutions. These can
be re-enabled using the `legend` plotting keyword. However, when plotting a
large number of trajectories, each will generate a label. Sometimes the best
approach is to remove these and add a label manually:
```@example ex1
p = plot(esol; label = nothing)
plot!(p, Float64[], Float64[]; label = "X", legend = :topleft)
```

## [Event handling using callbacks](@id advanced_simulations_callbacks)
Sometimes one wishes to add discrete events during simulations. Examples could include:
 - A chemical system where an amount of some species is added at a time point
   after the simulation's initiation.
 - A simulation of a circadian rhythm, where light is turned on/off every 12 hours.
 - A cell divides when some size variable reaches a certain threshold, randomly
   allocating all species to two daughter cells.

In simple cases events such as these can be modelled symbolically, as described
in the [Constraint Equations and Events](@ref constraint_equations) tutorial. A
more flexible, but low-level, interface is also available via the callback
functionality of DifferentialEquations.jl. A callback is a function that is
passed to the `solve()` command, combing an `affect!` function (defining how the
callback changes the system) with a `condition` function (a condition for
triggering a callback). For a thorough introduction, please read [the section
about callbacks in the DifferentialEquations.jl
documentation](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/).

There exist three types of callbacks, `PresetTimeCallback`s `DiscreteCallback`s,
and `ContinuousCallback`s. Here, we will limit ourselves to introducing the
`PresetTimeCallback`. For our example, we are going to use a simple network
where a single component, `X`, degrades linearly.
```@example ex2
using Catalyst
degradation_model = @reaction_network begin
    d, X --> 0
end
nothing # hide
```
we can simulate the model without using a callback:
```@example ex2
using DifferentialEquations, Plots
u0 = [:X => 10.0]
tspan = (0.0, 10.0)
p = [:d => 1.0]

oprob = ODEProblem(degradation_model, u0, tspan, p)
sol = solve(oprob)
plot(sol)
```
We now wish to modify our simulation so that at the times `t = 3.0` and `t =
7.0` we add `5` units of `X` to the system. For this we create a
`PresetTimeCallback`:
```@example ex2
condition = [3.0, 7.0]
function affect!(integrator)
    integrator[:X] += 5.0
end
ps_cb = PresetTimeCallback(condition, affect!)
nothing # hide
```
Here, `condition` is simply a vector with all the time points during which we
want the callback to trigger. The `affect!` function determines what happens to
the simulation when the callback is triggered. It takes a single object, an
`integrator` and makes some modification to it (please read more about
integrators [here](https://docs.sciml.ai/DiffEqDocs/stable/basics/integrator/)).
Here, we access the system's unknown `X` through the `integrator`, and add
`5.0` to its amount. We can now simulate our system using the
callback:
```@example ex2
sol = solve(oprob, Tsit5(); callback = ps_cb)
plot(sol)
```

Next, we can also use a callback to change the parameters of a system. The following code
plots the concentration of a two-state system, as we change the equilibrium
constant between the two states:
```@example ex2
rn = @reaction_network begin
    (k,1), X1 <--> X2
end
u0 = [:X1 => 10.0, :X2 => 0.0]
tspan = (0.0, 20.0)
p = [:k => 1.0]
oprob = ODEProblem(rn, u0, tspan, p)

condition = [5.0]
affect!(integrator) = integrator.ps[:k] = 5.0
ps_cb = PresetTimeCallback(condition, affect!)

sol = solve(oprob, Tsit5(); callback = ps_cb)
plot(sol)
```
The result looks as expected. However, what happens if we attempt to run the
simulation again?
```@example ex2
sol = solve(oprob, Tsit5(); callback = ps_cb)
plot(sol)
```
The plot looks different, even though we simulate the same problem. Furthermore,
the callback does not seem to have any effect on the system. If we check our
`ODEProblem`
```@example ex2
oprob.p
```
we note that `k = 5.0`, rather than `k = 1.0` as we initially specified. This is
because the callback modifies our `ODEProblem` during the simulation, and this
modification remains during the second simulation. An improved workflow to avoid
this issue is:
```@example ex2
rn = @reaction_network begin
    (k,1), X1 <--> X2
end
u0 = [:X1 => 10.0,:X2 => 0.0]
tspan = (0.0, 20.0)
p = [:k => 1.0]
oprob = ODEProblem(rn, u0, tspan, p)

condition = [5.0]
affect!(integrator) = integrator.ps[:k] = 5.0
ps_cb = PresetTimeCallback(condition, affect!)

sol = solve(deepcopy(oprob), Tsit5(); callback = ps_cb)
plot(sol)
```
where we parse a copy of our `ODEProblem` to the solver (using `deepcopy`). We can now run
```@example ex2
sol = solve(deepcopy(oprob), Tsit5(); callback = ps_cb)
plot(sol)
```
and get the expected result.

It is possible to give several callbacks to the `solve()` command. To do so, one
has to bundle them together in a `CallbackSet`, here follows one example:
```@example ex2
rn = @reaction_network begin
    (k,1), X1 <--> X2
end
u0 = [:X1 => 10.0, :X2 => 0.0]
tspan = (0.0, 20.0)
p = [:k => 1.0]
oprob = ODEProblem(rn, u0, tspan, p)

ps_cb_1 = PresetTimeCallback([3.0, 7.0], integ -> integ[:X1] += 5.0)
ps_cb_2 = PresetTimeCallback([5.0], integ -> integ.ps[:k] = 5.0)

sol = solve(deepcopy(oprob), Tsit5(); callback=CallbackSet(ps_cb_1, ps_cb_2))
plot(sol)
```

The difference between the `PresetTimeCallback`s and the `DiscreteCallback`s and
`ContiniousCallback`s is that the latter two allow the condition to be a
function, permitting the user to give more general conditions for the callback
to be triggered. An example could be a callback that triggers whenever a species
surpasses some threshold value.

### [Callbacks during SSA simulations](@id advanced_simulations_ssa_callbacks)
An assumption of (most) SSA simulations is that the state of the system is unchanged between reaction events. However, callbacks that affect the system's state can violate this assumption. To prevent erroneous simulations, users must inform a SSA solver when the state has been updated in a callback. This allows the solver to reinitialize any internal state information that may have changed. This can be done through the `reset_aggregated_jumps!` function, see the following example:

```@example ex2
rn = @reaction_network begin
    (k,1), X1 <--> X2
end
u0 = [:X1 => 10.0,:X2 => 0.0]
tspan = (0.0, 20.0)
p = [:k => 1.0]
dprob = DiscreteProblem(rn, u0, tspan, p)
jprob = JumpProblem(rn, dprob, Direct())

condition = [5.0]
function affect!(integrator)
    integrator[:X1] += 5.0
    integrator.ps[:k] += 2.0
    reset_aggregated_jumps!(integrator)
    nothing
end
cb = PresetTimeCallback(condition, affect!)

sol = solve(deepcopy(jprob), SSAStepper(); callback=cb)
plot(sol)
```


## Scaling the noise magnitude in the chemical Langevin equations
When using the CLE to generate SDEs from a CRN, it can sometimes be desirable to
scale the magnitude of the noise terms. Here, each reaction of the system generates a separate noise term in the CLE. If you require identical scaling for all reactions, the `@default_noise_scaling` option can be used. Else, you can supply a `noise_scaling` metadata for each individual reaction, describing how to scale the noise for that reaction.

We begin with considering the first approach. First, we simulate a simple two-state CRN model using the
CLE:
```@example ex3
using Catalyst, StochasticDiffEq, Plots

rn_1 = @reaction_network begin
    (k1,k2), X1 <--> X2
end
u0 = [:X1 => 10.0, :X2 => 10.0]
tspan = (0.0, 10.0)
p_1 = [:k1 => 1.0, :k2 => 1.0]

sprob_1 = SDEProblem(rn_1, u0, tspan, p_1)
sol_1 = solve(sprob_1, ImplicitEM())
plot(sol_1; idxs = :X1, ylimit = (0.0, 20.0))
```
Here we can see that the $X$ concentration fluctuates around a steady state of $X≈10.0$.

Next, we wish increase the amount of noise by a factor 2. To do so, we use the `@default_noise_scaling` option, to which we provide the desired scaling 
```@example ex3
rn_2 = @reaction_network begin
    @default_noise_scaling 2
    (k1,k2), X1 <--> X2
end
```
If we re-simulate the system we see that the amount of noise have increased:
```@example ex3
sprob_1 = SDEProblem(rn_2, u0, tspan, p_1)
sol_1 = solve(sprob_1, ImplicitEM())
plot(sol_1; idxs = :X1, ylimit = (0.0, 20.0))
```

It is possible to scale the amount of noise using any expression. A common use of this is to set a parameter which determines the amount of noise. Here we create a parameter $η$, and uses its value to scale the noise.
```@example ex3
using Catalyst, StochasticDiffEq, Plots

rn_3 = @reaction_network begin
    @parameters η
    @default_noise_scaling η
    (k1,k2), X1 <--> X2
end
u0 = [:X1 => 10.0, :X2 => 10.0]
tspan = (0.0, 10.0)
p_3 = [:k1 => 1.0, :k2 => 1.0, :η => 0.2]

sprob_3 = SDEProblem(rn_3, u0, tspan, p_3)
sol_3 = solve(sprob_3, ImplicitEM())
plot(sol_3; idxs = :X1, ylimit = (0.0, 20.0))
```
Here we saw how, by setting a small $η$ value, the amount of noise was reduced.

It is possible to use a different noise scaling expression for each reaction. Here, each reaction's noise scaling expression is provided using the `noise_scaling` metadata. In the following example, we use this to tune the noise of for both reactions involving the species $Y$.

```@example ex3
rn_4 = @reaction_network begin
    (p, d), 0 <--> X
    (p, d), 0 <--> Y, ([noise_scaling=0.0], [noise_scaling=0.0])
end

u0_4 = [:X => 10.0, :Y => 10.0]
tspan = (0.0, 10.0)
p_4 = [:p => 10.0, :d => 1.]

sprob_4 = SDEProblem(rn_4, u0_4, tspan, p_4)
sol_4 = solve(sprob_4, ImplicitEM())
plot(sol_4; ylimit = (0.0, 20.0))
```
Here, we not that there is n fluctuation in the value of $Y$. If the `@default_noise_scaling` option is used, its value is used for all reactions for which the `noise_scaling` metadata is unused. If `@default_noise_scaling` is not used, the default noise scaling value is `1.0` (i.e. no scaling).

## Useful plotting options
Catalyst, just like DifferentialEquations, uses the Plots package for all
plotting. For a detailed description of differential equation plotting, see
[DifferentialEquations documentation on the
subject](https://docs.sciml.ai/DiffEqDocs/stable/basics/plot/). Furthermore, the
[Plots package documentation](https://docs.juliaplots.org/stable/) contains
additional information and describes [a large number of plotting
options](https://docs.juliaplots.org/stable/attributes/). Here follows a very
short tutorial with a few useful options.

Let us consider the Brusselator model:
```@example ex4
using Catalyst, DifferentialEquations, Plots

brusselator = @reaction_network begin
    A, ∅ --> X
    1, 2X + Y --> 3X
    B, X --> Y
    1, X --> ∅
end
u0 = [:X => 1.0, :Y => 0.0]
tspan = (0.0, 50.0)
p = [:A => 1.0, :B => 4.0]

oprob = ODEProblem(brusselator, u0, tspan, p)
sol = solve(oprob)
plot(sol)
```
If we want to plot only the `X` species, we can use the `idxs` command:
```@example ex4
plot(sol; idxs = [:X])
```
If we wish to plot a single species (such as we do in this case), vector notation is not required and
we could simply write `plot(sol; idxs=:X)`.

Next, if we wish to plot a solution in phase space (instead of across time) we
again use the `idxs` notation, but use `()` instead of `[]` when designating the
species we wish to plot. Here, we plot the solution in `(X,Y)` space:
```@example ex4
plot(sol; idxs=(:X, :Y))
```

---
## References
[^1]: [DifferentialEquations.jl online documentation.](https://docs.sciml.ai/DiffEqDocs/stable/)
[^2]: [Chris Rackauckas, Qing Nie, *DifferentialEquations.jl – A Performant and Feature-Rich Ecosystem for Solving Differential Equations in Julia*, Journal of Open Resource Software (2017).](https://openresearchsoftware.metajnl.com/articles/10.5334/jors.151)