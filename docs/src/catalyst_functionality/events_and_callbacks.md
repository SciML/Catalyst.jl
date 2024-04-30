# [Simulation events](@id events)

Sometimes one wishes to model events, describing things that can happen to a system during its simulation. Examples can include:
 - A chemical system where an amount of some species is added at a time point after the simulation's initiation.
 - A simulation of a circadian rhythm, where light is turned on/off every 12 hours.
 - A cell divides when some size variable reaches a certain threshold, halving its size.
  
Generally, there are two approaches for creating events:
1. Including them as a part of the model when it is created.
2. Creating a separate *callback*, which is then supplied to the `solve` function.

Generally, the first approach is more convenient, and should be used if possible. However, the second approach is more potent, and can be used to create more general types of events. Below, we will describe both approaches.

!!! warn
    Independently of how they are created, there are some additional considerations when creating events. These are described [at the end of this tutorial](@ref ref). Please briefly check these cases, and if any apply to you, read the corresponding section thoroughly.

## [Creating models with incorporated events](@id events_model)
Catalyst's DSL provides two options, `@continuous_events` and `@discrete_events`, which permits users to add events to the models when they are constructed. Each corresponds to the creation of a specific type of events, and both are described both. The syntax for creating evens in the DSL is identical to that used for [programmatic] model creation (described in detail [here]) and which is used by the [ModelingToolkit.jl package](https://docs.sciml.ai/ModelingToolkit/stable/basics/Events/)

Both continuous and discrete events combine a *condition* (for triggering the event) with an *affect* (describing its effect on the system). They differ in the following ways:
- They use slightly different notation.
- Discrete events' conditions are checked at *the end of* each simulation time step. For continuous events, the simulation instead finds the *exact time point* when the event is triggered at.
- Continuous events cannot be supplied to jump simulations.

### [Continuous events](@id events_model_continuous)
Let us consider a simple system where species `X` degrades at a constant rate `d`. We wish to create an event that adds `2.0` units of `X` whenever `X` reaches a critical threshold `1.0`. This can be done in the following manner:
```@example dsl_advanced_events
using Catalyst # hide
rn = @reaction_network begin
  @continuous_events begin
    [X ~ 1.0] => [X ~ X + 2.0]
  end
  d, X --> 0
end
nothing # hide
```
Here, the `@continuous_events` option is followed by a `begin ... end` block. Next, each line corresponds to a separate event. Each continuous event is created in the following manner:
- It combines a set of *conditions* (describing when the event is triggered) with a set of *affects* (describing the event's effect on the system once triggered). 
- The conditions and affects are vectors (containing each containing any number of condition or affect equations).
- The conditions (left) and affects (right) are separated by `=>`.
- Each condition is an [equation](@ref ref). Here, the event is triggered whenever any of its conditions equations' hold (i.e. the equation's left and right ahnd sides are identical).
- Each affect is a single [equation](@ref ref) that describes how a parameter or species (or [variable](@ref ref)) is updated when the event is triggered.
- Each affect's equation's left-hand side must contain only the parameter/species(/variable) whose value should be updated.
- Each affect's equation's right-hand side is an expression describing its updated value.
- A single event's affects can updates species (and variables) only, or parameters only (but not parameters and species/variables).

The declared model can be simulated using standard syntax:
```@example dsl_advanced_events
using OrdinaryDiffEq, Plots
u0 = [:X => 2.0]
tspan = (0.0, 10.0)
ps = [:d => 1.0]

oprob = ODEProblem(rn, u0, tspan, ps)
sol = solve(ODEProblem, Tsit5())
plot(sol)
```
Inspecting the solution, we can confirm that whenever `X` reaches a value of `1.0`, `2.0` units of `X` is added to the system.

In our example, we can also denote the critical quantities using parameters:
```@example dsl_advanced_events
rn = @reaction_network begin
  @parameters X_thres X_add
  @continuous_events begin
    [X ~ X_thres] => [X ~ X + X_add]
  end
  d, X --> 0
end
nothing # hide
```
Here, since `X_thres` and `X_add` do not appear in any reactions, Catalyst cannot determine whether they are parameters or species. Hence, they must be [explicitly designated as parameters by using the `@parameters` option](@ref dsl_advanced_options_declaring_species_and_parameters). These values can now be designated when the `ODEProblem` is created:
```@example dsl_advanced_events
u0 = [:X => 2.0]
tspan = (0.0, 10.0)
ps = [:d => 1.0, :X_thres => 0.5, :X_add => 3.0]

oprob = ODEProblem(rn, u0, tspan, ps)
sol = solve(ODEProblem, Tsit5())
plot(sol)
```

As previously noted, each continuous event can have multiple affects. The following system has two components (`X` and `Y`, one being produced and one being degraded). When their concentrations are equal, a continuous events reduce the concentration of `X` while increasing the concentration of `Y`:
```@example dsl_advanced_events
rn = @reaction_network begin
  @continuous_events begin
    [X ~ Y] => [X ~ X - 1.0, Y ~ Y + 1.0]
  end
  p, 0 --> X
  d, Y --> 0
end

u0 = [:X => 1.0, :Y => 3.0]
tspan = (0.0, 10.0)
ps = [:p => 1.0, :d => 1.0]

oprob = ODEProblem(rn, u0, tspan, ps)
sol = solve(ODEProblem, Tsit5())
plot(sol)
```

In the above examples we have modelled a system with a single event. In these cases, the `begin end` block is not required, and the event can be provided on the same line as the `@continuous_events` option:
```@example dsl_advanced_events
rn = @reaction_network begin
  @continuous_events  X ~ Y => [X ~ X - 1.0, Y ~ Y + 1.0]
  p, 0 --> X
  d, Y --> 0
end
nothing # hide
```

### [Discrete events](@id events_model_discrete)
Discrete events are similar to continuous events in that they combine a condition for being triggered with an affect once triggered. The affect(s) of discrete events are declared, and works, identical to [those of continuous events](@ref events_model_continuous). However, the condition is different. There exist 3 different types of discrete events, each with a different type of condition. All three types are created using the `@discrete_events` option, and a single system can contain a mix of all types. The three types are:
- Preset-time discrete events.
- Periodic discrete events.
- Conditional discrete events.

#### [Preset-time discrete events](@id events_model_discrete_presettime)
*Present-time events* are events that happen at specific time points. Their conditions are vectors with all the time points at which the event is triggered. E.g. here we create a production/degradation loop, where `2.0` units of `X` is added at time points `3.0` and `7.0`
```@example dsl_advanced_events
rn = @reaction_network begin
  @discrete_events begin
    [3.0, 7.0] => [X ~ X + 2.0]
  end
  (p,d), 0 <--> X
end

u0 = [:X => 0.1, :Y => 3.0]
tspan = (0.0, 10.0)
ps = [:p => 1.0, :d => 0.5]

oprob = ODEProblem(rn, u0, tspan, ps)
sol = solve(ODEProblem, Tsit5())
plot(sol)
```

Currently, the time points of preset cannot be parameters. If you want to create preset time event with a callback, either use a [continuous callback with condition `[t ~ preset_t]`](@ref events_model_continuous) (note, not possible for jump simulations), or [create a conditional discrete event and supply a `tstops` to the `solve` command](@ref events_additional_considerations_tstops)

#### [Periodic discrete events](@id events_model_discrete_periodic)
When a discrete event's condition is a vector, a preset-time event is created. If it instead is a single value, a *periodic event* is created. These occur repeatedly throughout a simulation, with its period set by the affect term. E.g. here we create a system where `0.5` units of `X` is added every `1.0` time units.
```@example dsl_advanced_events
rn = @reaction_network begin
  @discrete_events begin
    1.0 => [X ~ X + 0.5]
  end
  d, X --> 0
end

u0 = [:X => 1.0]
tspan = (0.0, 10.0)
ps = [:d => 1.0]

oprob = ODEProblem(rn, u0, tspan, ps)
sol = solve(ODEProblem, Tsit5())
plot(sol)
```
Like for preset-time events, periodic events' conditions cannot contain parameters.

#### [Conditional discrete events](@id events_model_discrete_conditional)
Finally, discrete events' condition may be a boolean expression (consisting of any combinations of parameters, species, variables, numbers, and the time variable). Let's say that we want to create an event which, if the concentration of `X` is below a threshold `1.0`, adds `1.0` units of `X` to the system, then we can use the following discrete event:
```@example dsl_advanced_events
rn = @reaction_network begin
  @discrete_events begin
    X < 1.0 => [X ~ X + 2.0]
  end
  d, X --> 0
end
```
If we simulate the system using the same conditions as for our [similar, continuous, example](@ref dsl_advanced_options_events_continuous), the result is similar:
```@example dsl_advanced_events
u0 = [:X => 2.0]
tspan = (0.0, 10.0)
ps = [:d => 1.0]

oprob = ODEProblem(rn, u0, tspan, ps)
sol = solve(ODEProblem, Tsit5())
plot(sol)
```
So, how is modelling this event as a discrete or continuous event different? There are four differences:
1) For continuous events, the simulation method finds the exact time point when the condition triggers. Discrete events are triggered after a time step after which the condition holds.
2) This discrete event will be triggered whenever `X < 1.0` holds, not just when the concentration of `X` passes the threshold. E.g. it will be triggered if the initial concentration of `X` is less than `1.0`.
3) Only discrete event events can be used with jump simulations.
4) Discrete event can be used to create more advanced conditions.

E.g. using (4), we can modify our system so that the event is only triggered when time is less than `5.0` (after which `X` decays towards `0`):
```@example dsl_advanced_events
rn = @reaction_network begin
  @discrete_events begin
    (X < 1.0) & (t < 5.0) => [X ~ X + 2.0]
  end
  d, X --> 0
end

u0 = [:X => 2.0]
tspan = (0.0, 10.0)
ps = [:d => 1.0]

oprob = ODEProblem(rn, u0, tspan, ps)
sol = solve(ODEProblem, Tsit5())
plot(sol)
```
!!! note
    When we form composite boolean conditions for conditional discrete events, we use `&` to denote the AND operator (not `&&`, as this is currently not supported). Similarly, `|` is sued instead of `||`.

!!! warn
    Generally, discrete events including equality (`==`) will not be triggered. The reason is that the condition is only checked at the end of every time step. Hence, unless special precautions are taken to ensure that the simulator stops when the condition holds, it will unlikely be triggered.


## [Creating event through callbacks](@id events_model_callbacks)

### [Discrete callbacks](@id events_model_callbacks_discrete)

### [Present time callbacks](@id events_model_callbacks_presenttime)

### [Continuous callbacks](@id events_model_callbacks_continuous)

### [Callbacks during jump simulations](@id events_model_callbacks_jump_callbacks)


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
sol = solve(oprob; callback = ps_cb)
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
affect!(integrator) = integrator[:k] = 5.0
ps_cb = PresetTimeCallback(condition, affect!)

sol = solve(oprob; callback = ps_cb)
plot(sol)
```
The result looks as expected. However, what happens if we attempt to run the
simulation again?
```@example ex2
sol = solve(oprob; callback = ps_cb)
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
affect!(integrator) = integrator[:k] = 5.0
ps_cb = PresetTimeCallback(condition, affect!)

sol = solve(deepcopy(oprob); callback = ps_cb)
plot(sol)
```
where we parse a copy of our `ODEProblem` to the solver (using `deepcopy`). We can now run
```@example ex2
sol = solve(deepcopy(oprob); callback = ps_cb)
plot(sol)
```
and get the expected result.

It is possible to give several callbacks to the `solve()` command. To do so, one
has to bundle them together in a `CallbackSet`, here follows one example:
```@example ex2
rn = @reaction_network begin
    (k,1), X1 <--> X2
end
u0 = [:X1 => 10.0,:X2 => 0.0]
tspan = (0.0, 20.0)
p = [:k => 1.0]
oprob = ODEProblem(rn, u0, tspan, p)

ps_cb_1 = PresetTimeCallback([3.0, 7.0], integ -> integ[:X1] += 5.0)
ps_cb_2 = PresetTimeCallback([5.0], integ -> integ[:k] = 5.0)

sol = solve(deepcopy(oprob); callback=CallbackSet(ps_cb_1, ps_cb_2))
plot(sol)
```

The difference between the `PresetTimeCallback`s and the `DiscreteCallback`s and
`ContiniousCallback`s is that the latter two allow the condition to be a
function, permitting the user to give more general conditions for the callback
to be triggered. An example could be a callback that triggers whenever a species
surpasses some threshold value.

## [Additional consideration when using events (or callbacks)](@id events_additional_considerations)

### [Events affecting changing model parameter values](@id events_additional_considerations_parameters)

### [Events affecting changing model parameter values during jump simulations](@id events_additional_considerations_jump_parameters)

### [Events occurring at specific time points](@id events_additional_considerations_tstops)

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
    integrator[:k] += 2.0
    reset_aggregated_jumps!(integrator)
    nothing
end
cb = PresetTimeCallback(condition, affect!)

sol = solve(deepcopy(jprob), SSAStepper(); callback=cb)
plot(sol)
```