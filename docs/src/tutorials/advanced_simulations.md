# [Advanced Simulation Options](@id advanced_simulations)
Throughout the preceding tutorials, we have shown the basics of how to solve ODE, SDE, and jump process models generated from Catalyst `ReactionSystem`s. In this tutorial we'll illustrate some more advanced functionality that can be useful in many modeling contexts, and that provide conveniences for common workflows. For a comprehensive overview of solver properties, parameters, and manipulating solution objects, please read the [documentation of the DifferentialEquations package](https://docs.sciml.ai/DiffEqDocs/stable/), which Catalyst uses for all simulations. 



### Monte Carlo simulations using `EnsembleProblem`s
In many contexts one needs to run multiple simulations of a model, for example to collect statistics of SDE or jump process solutions, or to systematically vary parameter values within a model. While it is always possible to manually run such ensembles of simulations via a `for` loop, DifferentialEquations.jl provides the `EnsembleProblem` as a convenience to manage structured collections of simulations. `EnsembleProblem`s provide a simple interface for modifying a problem between individual simulations, and offers several options for batching and/or parallelizing simulation runs. For a more thorough description, please read [the Parallel Ensemble Simulations section of the DifferentialEquations documentation](https://docs.sciml.ai/DiffEqDocs/stable/features/ensemble/#ensemble). Here, we will give a brief introduction to the use of `EnsembleProblem`s from Catalyst-generated models.

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
sprob = SDEProblem(rn,u0,tspan,p)
nothing # hide
```
we can then use our `SDEProblem` as input to an `EnsembleProblem`:
```@example ex1
eprob = EnsembleProblem(sprob)
```
The `EnsembleProblem` can now be used as input to the `solve` command. It has the same options as when simulating the `SDEProblem` directly, however, it has an additional argument `trajectories` to determine how many simulations should be performed. 
```@example ex1
esol = solve(eprob; trajectories=5)
```
This simulation is automatically multithreaded over all available threads. Please read [this documentation](https://docs.sciml.ai/DiffEqDocs/stable/features/ensemble/#EnsembleAlgorithms) for more information on parallelisation alternatives. The ensemble simulations can be plotted using the `plot` function, which by default displays all trajectories:
```@example ex1
plot(esol)
```

Sometimes when performing a large number of ensemble simulations, the plots get very dense. In these cases, the plot argument `linealpha` (which sets trajectory transparency) may be useful:
```@example ex1
esol = solve(eprob; trajectories=100)
plot(esol; linealpha=0.5)
```

Sometimes, one wishes to perform the same simulation a large number of times, while making minor modifications to the problem each time. This can be done by giving a problem function, `prob_func`, argument to the `EnsembleProblem`. Let us consider ODE simulations of a simple birth/death process:
```@example ex1
rn = @reaction_network begin
    (b,1.0), ∅ <--> X
end
u0 = [:X => 1.0]
tspan = (0.0,1.0)
p = [:b => 1.];
oprob = ODEProblem(rn,u0,tspan,p)
nothing # hide
```
We wish to simulate this model for a large number of values of `b`. We do this by creating a `prob_func` that will make a modification to the problem at the start of each Monte Carlo simulation:
```@example ex1
b_values = 1.0:0.1:2.0
function prob_func(prob,i,repeat)
    remake(prob; p=[b_values[i]])
end
nothing # hide
```
Here, `prob_func` takes three arguments:
 - `prob`: The problem given to our `EnsembleProblem`, this is the problem that `prob_func` modifies in each iteration.
 - `i`: The number of this specific Monte Carlo iteration in the interval `1:trajectories`.
 - `repeat`: The repeat of this specific Monte Carlo simulation (We will ignore this argument in this brief overview).
In our case, for each Monte Carlo simulation, our `prob_func` takes our original `ODEProblem` and uses the `remake` function to change the parameter vector. Here, for the `i`th Monte Carlo simulation, the value of `b` is also the `i`th value of our `b_values` vector. Finally, we can simulate and plot our problem:
```@example ex1
eprob = EnsembleProblem(oprob; prob_func=prob_func)
esol = solve(eprob; trajectories=length(b_values))
plot(esol)
```

Note that plot legends are disabled when plotting ensemble solutions. These can be re-enabled using the `legend` plotting keyword. However, when plotting a large number of trajectories, each will generate a label. Sometimes the best approach is to remove these and add a label manually:
```
plot(esol; label="")
plot!([],[],label="X", legend=:best)
```


### Event handling using callbacks
Sometimes one wishes to add discrete events during simulations. Examples could include:
 - A chemical system where an amount of some species is added at a time point after the simulation's initiation.
 - A simulation of a circadian rhythm, where light is turned on/off every 12 hours.
 - A cell divides when some size variable reaches a certain threshold, randomly allocating all species to two daughter cells.
Events such as these can be modelled using callbacks. A callback is a function that is parsed to a `solve()` command, combing an `affect!` function (defining how the callback changes the system) with a `condition` function (a condition for triggering a callback). For a throughout the introduction, please read [the section about callbacks in the DifferentialEquations.jl documentation](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/).

There exist three types of callbacks, `PresetTimeCallback`s `DiscreteCallback`s, and `ContinuousCallback`s. Here, we will limit ourselves to introducing the `PresetTimeCallback`. For our example, we are going to use a simple network where a single component, `X`, degrades linearly.
```@example ex2
using Catalyst
degradation_model = @reaction_network begin
    d, X --> 0
end
nothing # hide
```
we can simulate the model without using a callback:
```@example ex2
using DifferentialEquations
u0 = [:X => 10.0]
tspan = (0.0,10.0)
p = [:d => 1.0]

oprob = ODEProblem(degradation_model, u0, tspan, p)
sol = solve(oprob)
plot(sol)
```
We now wish to modify our simulation so that at the times `t=3.0` and `t=7.0` we add `5` units of `X` to the system. For this we create a `PresetTimeCallback`:
```@example ex2
condition = [3.0, 7.0]
function affect!(integrator)
    integrator.u[1] += 5.0
end
ps_cb = PresetTimeCallback(condition, affect!)
nothing # hide
```
Here, `condition` is simply a vector with all the time points during which we want the callback to trigger. The `affect!` function determines what happens to the simulation when the callback is triggered. It takes a single object, an `integrator` and makes some modification to it (please read more about integrators [here](https://docs.sciml.ai/DiffEqDocs/stable/basics/integrator/)). Here, we access the system's current state vector as `integrator.u`, and add `5.0` to the amount of `X` present. We can now simulate our system using the callback:
```@example ex2
sol = solve(oprob; callback=ps_cb)
plot(sol)
```

Next, we can also use a callback to change the parameters of a system. Now, instead of accessing `integrator.u` we access `integrator.p`. The following code plots the concentration of a two-state system, as we change the equilibrium constant between the two states:
```@example ex2
rn = @reaction_network begin
    (k,1), X1 <--> X2
end
u0 = [:X1 => 10.0,:X2 => 0.0]
tspan = (0.0,20.0)
p = [:k => 1.0]
oprob = ODEProblem(rn, u0, tspan, p)

condition = [5.0]
affect!(integrator) = integrator.p[1] = 5.0
ps_cb = PresetTimeCallback(condition, affect!)

sol = solve(oprob; callback=ps_cb)
plot(sol)
```
The result looks as expected. However, what happens if we attempt to run the simulation again?
```@example ex2
sol = solve(oprob; callback=ps_cb)
plot(sol)
```
The plot looks different, even though we simulate the same problem. Furthermore, the callback does not seem to have any effect on the system. If we check our `ODEProblem`
```@example ex2
oprob.p
```
we note that `k=5.0`, rather than `k=1.0` as we initially specify. This is because the callback modifies our `ODEProblem` during the simulation, and this modification remains during the second simulation. An improved workflow to avoid this issue is:
```@example ex2
rn = @reaction_network begin
    (k,1), X1 <--> X2
end
u0 = [:X1 => 10.0,:X2 => 0.0]
tspan = (0.0,20.0)
p = [:k => 1.0]
oprob = ODEProblem(rn, u0, tspan, p)

condition = [5.0]
affect!(integrator) = integrator.p[1] = 5.0
ps_cb = PresetTimeCallback(condition, affect!)

sol = solve(deepcopy(oprob); callback=ps_cb)
plot(sol)
```
where we parse a copy of our `ODEProblem` to the solver. We can now run
```@example ex2
sol = solve(deepcopy(oprob); callback=ps_cb)
plot(sol)
```
and get the expected result.

It is possible to give several callbacks to the `solve()` command. To do so, one has to bundle them together in a `CallbackSet`, here follows one example:
```@example ex2
rn = @reaction_network begin
    (k,1), X1 <--> X2
end
u0 = [:X1 => 10.0,:X2 => 0.0]
tspan = (0.0,20.0)
p = [:k => 1.0]
oprob = ODEProblem(rn, u0, tspan, p)

ps_cb_1 = PresetTimeCallback([3.0,7.0], integ -> integ.u[1] += 5.0)
ps_cb_2 = PresetTimeCallback([5.0], integ -> integ.p[1] = 5.0)

sol = solve(deepcopy(oprob); callback=CallbackSet(ps_cb_1, ps_cb_2))
plot(sol)
```

The difference between the `PresetTimeCallback`s and the `DiscreteCallback`s and `ContiniousCallback`s is that the latter two allow the condition to be a function, permitting the user to give more general conditions for the callback to be triggered. An example could be a callback that triggers whenever a species surpasses some threshold value.



### Scaling the noise magnitude in the chemical Langevin equations
When using the CLE to generate SDEs from a CRN, it can sometimes be desirable to scale the magnitude of the noise terms. This can be done by introducing a *noise scaling parameter*. First, we simulate a simple two-state CRN model using the CLE:
```@example ex3
using Catalyst, StochasticDiffEq, Plots

rn_1 = @reaction_network begin
    (k1,k2), X1 <--> X2
end
u0 = [:X1 => 10.0, :X2 => 10.0]
tspan = (0.0,10.0)
p_1 = [:k1 => 1.0, :k2 => 1.0]

sprob_1 = SDEProblem(rn_1,u0,tspan,p_1)
sol_1 = solve(sprob_1)
plot(sol_1; idxs=1, ylimit=(0.0,20.0))
```
Here we can see that the `X` concentration fluctuations around a steady state of *X≈10.0*. 

Next, we wish to introduce a noise scaling parameter ,`η`. This will scale the noise magnitude so that for *η≈0.0* the system lacks noise (and its SDE simulations are identical to its ODE simulations) and for *η≈1.0* noise is not scaled (and SDE simulations are identical to as if no noise scaling was used). Setting *η<1.0* will reduce noise and *η>1.0* will increase noise. The syntax for setting a noise scaling parameter `η` is
```@example ex3
rn_2 = @reaction_network begin
    @parameters η
    (k1,k2), X1 <--> X2
end
u0 = [:X1 => 10.0, :X2 => 10.0]
tspan = (0.0,10.0)
p_2 = [:k1 => 1.0, :k2 => 1.0, :η => 0.1]

sprob_2 = SDEProblem(rn_2,u0,tspan,p_2; noise_scaling=(@parameters η)[1])
```
Here, we first need to add `η` as a parameter to the system using the `@parameters η` option. Next, we pass the `noise_scaling=(@parameters η)[1]` argument to the `SDEProblem`. We can now simulate our system and confirm that noise is reduced:
```@example ex3
sol_2 = solve(sprob_2)
plot(sol_2; idxs=1, ylimit=(0.0,20.0))
```

Finally, it is possible to set individual noise scaling parameters for each reaction of the system. Our model has two reactions (`X1 --> X2` and `X2 --> X1`) so we will use two noise scaling parameters, `η1` and `η2`. We use the following syntax:
```@example ex3
rn_3 = @reaction_network begin
    @parameters η1 η2
    (k1,k2), X1 <--> X2
end
u0 = [:X1 => 10.0, :X2 => 10.0]
tspan = (0.0,10.0)
p_3 = [:k1 => 1.0, :k2 => 1.0, :η1 => 0.1, :η2 => 1.0]

sprob_3 = SDEProblem(rn_3,u0,tspan,p_3; noise_scaling=@parameters η1 η2)
```
plotting the results, we see that we have less fluctuation than for the first simulation, but more as compared to the second one (which is as expected):
```@example ex3
sol_3 = solve(sprob_3)
plot(sol_3; idxs=1, ylimit=(0.0,20.0))
```



### Useful plotting options
Catalyst, just like DifferentialEquations, uses the Plots package for all plotting. For a detailed description of differential equation plotting, see [DifferentialEquations documentation on the subject](https://docs.sciml.ai/DiffEqDocs/stable/basics/plot/). Furthermore, the [Plots package documentation](https://docs.juliaplots.org/stable/) contains additional information and describes [a large number of plotting options](https://docs.juliaplots.org/stable/attributes/). Here follows a very short tutorial with a few useful options.

Let us consider the Brusselator model:
```@example ex4
using Catalyst, DifferentialEquations, Plots

brusselator = @reaction_network begin
    A, ∅ → X
    1, 2X + Y → 3X
    B, X → Y
    1, X → ∅
end
u0 = [:X => 1.0, :Y => 0.0]
tspan = (0.0,50.0)
p = [:A => 1.0, :B => 4.0]

oprob = ODEProblem(brusselator, u0, tspan, p)
sol = solve(oprob)
plot(sol)
```
If we want to plot only the `X` species, we can use the `idxs` command:
```@example ex4
@unpack X = brusselator
plot(sol; idxs=[X])
```
Here we use the `brusselator.X` notation to denote that we wish to plot the `X` species. The input to `idxs` is a vector listing all the species we wish to plot. If we wish to plot a single species, vector notation is not required and we could simply write `plot(sol; idxs=brusselator.X)`. 
(The plotting feature that automatically sets the label when using this interface is currently not working optimally, hence we manually set the label using the `label` option.)

Next, if we wish to plot a solution in phase space (instead of across time) we again use the `idxs` notation, but use `()` instead of `[]` when designating the species we wish to plot. Here, we plot the solution in `(X,Y)` space:
```@example ex4
@unpack X,Y = brusselator
plot(sol; idxs=(X,Y))
```

