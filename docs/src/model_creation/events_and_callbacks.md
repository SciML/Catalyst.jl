# [Simulation events](@id events)

Sometimes one wishes to model events, describing things that can happen to a system during its simulation. Examples can include:
 - A chemical system where an amount of some species is added at a time point after the simulation's initiation.
 - A simulation of a circadian rhythm, where light is turned on/off every 12 hours.
 - A cell divides when a size variable reaches a certain threshold (subsequently halving the size).
  
Generally, there are two approaches for creating events:
1. Including them as a part of the model when it is created.
2. Creating a separate *callback*, which is then supplied to the `solve` function.

Generally, the first approach is more convenient and should be used if possible. However, the second approach is more potent and can be used to create more general types of events. Below, we will describe both approaches.

!!! warn
    Independently of how they are created, there are some additional considerations when creating events. These are described [at the end of this tutorial](@ref events_additional_considerations). Please briefly check these cases, and if any apply to you, read the corresponding section thoroughly.

## [Creating models with incorporated events](@id events_model)
Catalyst's DSL provides two options, `@continuous_events` and `@discrete_events`, which permit users to add events to the models when they are constructed. Each corresponds to the creation of a specific type of event, and both are described below. The syntax for creating events in the DSL is identical to that used for [programmatic model creation]() (described in detail [here](@ref events_programmatic)) and which is used by the [ModelingToolkit.jl package](https://docs.sciml.ai/ModelingToolkit/stable/basics/Events/)

Both continuous and discrete events combine a *condition* (for triggering the event) with an *affect* (describing its effect on the system). They differ in the following ways:
- They use slightly different notation.
- Discrete events' conditions are checked at *the end of* each simulation time step. For continuous events, the simulation instead finds the *exact time point* when the event is triggered at.
- Continuous events cannot be supplied to jump simulations.

### [Continuous events](@id events_model_continuous)
Let us consider a simple system where species `X` degrades at a constant rate `d`. We wish to create an event that adds `2.0` units of `X` whenever `X` reaches a critical threshold `1.0`. This can be done in the following manner:
```@example events_model
using Catalyst
decay_model = @reaction_network begin
  @continuous_events begin
    [X ~ 1.0] => [X ~ X + 2.0]
  end
  d, X --> 0
end
nothing # hide
```
Here, the `@continuous_events` option is followed by a `begin ... end` block. Next, each line corresponds to a separate event. Each continuous event is created in the following manner:
- It combines a set of *conditions* (describing when the event is triggered) with a set of *affects* (describing the event's effect on the system once triggered). 
- The conditions and affects are vectors (each containing any number of condition or affect equations).
- The conditions (left) and affects (right) are separated by `=>`.
- Each condition is an [equation](@ref ref). Here, the event is triggered whenever any of its condition equations' hold (i.e. the equation's left and right-hand sides are identical).
- Each affect is a single [equation](@ref ref) that describes how a parameter or species (or [variable](@ref ref)) is updated when the event is triggered.
- Each affect's equation's left-hand side must contain only the parameter/species(/variable) whose value should be updated.
- Each affect's equation's right-hand side is an expression describing its updated value.
- A single event's affects can update species (and variables) only, or parameters only (but not parameters and species/variables).

Models containing events can be simulated using the standard syntax:
```@example events_model
using OrdinaryDiffEq, Plots
u0 = [:X => 2.0]
tspan = (0.0, 10.0)
ps = [:d => 0.5]

oprob = ODEProblem(decay_model, u0, tspan, ps)
sol = solve(oprob)
plot(sol)
```
Inspecting the solution, we can confirm that whenever `X` reaches a value of `1.0`, `2.0` units of `X` are added to the system.

In our example, we can also denote the critical quantities using parameters:
```@example events_model
decay_model = @reaction_network begin
  @parameters X_thres X_add
  @continuous_events begin
    [X ~ X_thres] => [X ~ X + X_add]
  end
  d, X --> 0
end
nothing # hide
```
Here, since `X_thres` and `X_add` do not appear in any reactions, Catalyst cannot determine whether they are parameters or species. Hence, they must be [explicitly designated as parameters by using the `@parameters` option](@ref ref). These values can now be designated when the `ODEProblem` is created:
```@example events_model
u0 = [:X => 2.0]
tspan = (0.0, 10.0)
ps = [:d => 0.5, :X_thres => 0.5, :X_add => 3.0]

oprob = ODEProblem(decay_model, u0, tspan, ps)
sol = solve(oprob)
plot(sol)
```

As previously noted, each continuous event can have multiple affects (and/or multiple conditions). The following system has two components (`X` and `Y`, one being degraded and one being produced). When their concentrations are equal, a continuous event reduces the concentration of `X` while also increasing the concentration of `Y`:
```@example events_model
XY_model = @reaction_network begin
  @continuous_events begin
    [X ~ Y] => [X ~ X - 1.0, Y ~ Y + 1.0]
  end
  p, 0 --> X
  d, Y --> 0
end

u0 = [:X => 1.0, :Y => 3.0]
tspan = (0.0, 10.0)
ps = [:p => 1.0, :d => 1.0]

oprob = ODEProblem(XY_model, u0, tspan, ps)
sol = solve(oprob)
plot(sol)
```

In the above examples we have modelled a system with a single event. In these cases, the `begin ... end` block is not required, and the event can be provided on the same line as the `@continuous_events` option:
```@example events_model
XY_model = @reaction_network begin
  @continuous_events  [X ~ Y] => [X ~ X - 1.0, Y ~ Y + 1.0]
  p, 0 --> X
  d, Y --> 0
end
nothing # hide
```

### [Discrete events](@id events_model_discrete)
Discrete events are similar to continuous events in that they combine a condition for being triggered with one (or several) affect(s) once triggered. The affect(s) of discrete events are declared, and works, identically to [those of continuous events](@ref events_model_continuous). However, the condition is different. There exist 3 different types of discrete events, each with a different type of condition. All three types are created using the `@discrete_events` option, and a single system can contain a mix of all types. The three types are:
- Preset-time discrete events.
- Periodic discrete events.
- Conditional discrete events.

#### [Preset-time discrete events](@id events_model_discrete_presettime)
*Present-time events* are events that happen at specific time points. Their conditions are vectors with all the time points at which the event is triggered. E.g. in this model, `2.0` units of `X` are added at time points `3.0` and `7.0`
```@example events_model
decay_model = @reaction_network begin
  @discrete_events begin
    [3.0, 7.0] => [X ~ X + 2.0]
  end
  d, X --> 0
end

u0 = [:X => 0.1]
tspan = (0.0, 10.0)
ps = [:d => 0.5]

oprob = ODEProblem(decay_model, u0, tspan, ps)
sol = solve(oprob)
plot(sol)
```

Currently, preset-time events' condition time points cannot be parameters. If you want to create a preset-time event with parameters, use a [continuous event with condition `[t ~ preset_t]`](@ref events_model_continuous) (not possible for jump simulations).

#### [Periodic discrete events](@id events_model_discrete_periodic)
When a discrete event's condition is a vector, a preset-time event is created. If it instead is a single value, a *periodic event* is created. These occur repeatedly throughout a simulation, with the period determined by the condition term. E.g. here we create a system where `0.5` units of `X` is added every `2.0` time units.
```@example events_model
decay_model = @reaction_network begin
  @discrete_events begin
    1.0 => [X ~ X + 0.5]
  end
  d, X --> 0
end

u0 = [:X => 2.0]
tspan = (0.0, 10.0)
ps = [:d => 0.5]

oprob = ODEProblem(decay_model, u0, tspan, ps)
sol = solve(oprob)
plot(sol)
```
As for preset-time events, periodic events' conditions cannot contain parameters.

#### [Conditional discrete events](@id events_model_discrete_conditional)
Finally, a discrete event's condition may be a boolean expression (consisting of any combinations of parameters, species, variables, numbers, and the time variable). Let's say that we want to create an event which, if the concentration of `X` is below a threshold `1.0`, adds `1.0` units of `X` to the system, then we can use the following discrete event:
```@example events_model
decay_model = @reaction_network begin
  @discrete_events begin
    (X < 1.0) => [X ~ X + 2.0]
  end
  d, X --> 0
end
```
If we simulate the system using the same conditions as for our [similar, continuous, example](@ref dsl_advanced_options_events_continuous), the result is similar
```@example dsl_advanced_events
u0 = [:X => 2.0]
tspan = (0.0, 10.0)
ps = [:d => 0.5]

oprob = ODEProblem(decay_model, u0, tspan, ps)
sol = solve(oprob)
plot(sol)
```
but not identical. This is explained by continuous events' ability to trigger exactly when the condition holds, while this discrete event triggers at the first time step when the condition holds (we can see this in the simulation, where it is randomly triggered sometime after the critical threshold `1.0` have been passed). In addition to this, discrete and continuous events have the following differences:
1) This discrete event will be triggered whenever `X < 1.0` holds, not just when the concentration of `X` passes the threshold. E.g. it will be triggered if the initial concentration of `X` is less than `1.0`.
2) Only discrete events can be used with jump simulations.
3) Discrete events can depend on more advanced conditions.

E.g. using (3), we can modify our system so that the event is only triggered when time is less than `5.0` (after which `X` decays towards `0`):
```@example events_model
decay_model = @reaction_network begin
  @discrete_events begin
    (X < 1.0) & (t < 5.0) => [X ~ X + 2.0]
  end
  d, X --> 0
end

oprob = ODEProblem(decay_model, u0, tspan, ps)
sol = solve(oprob)
plot(sol)
```
Again, the solver does not find the exact time points when the event is triggered. A workaround is to also supply a continuous event that is triggered just as the conditions holds, but which have no effect on the system. This way, the continuous event will ensure that the solver stops at the correct time points, at which point the discrete event will also be triggered:
```@example events_model
decay_model = @reaction_network begin
  @continuous_events [X ~ 0.9999] => [X ~ X]
  @discrete_events begin
    (X <= 1.0) & (t < 5.0) => [X ~ X + 2.0]
  end
  d, X --> 0
end

oprob = ODEProblem(decay_model, u0, tspan, ps)
sol = solve(oprob)
plot(sol)
```

Another alternative is to [force non-adaptive time stepping with a small `dt`](@ref ref):
```@example events_model
decay_model = @reaction_network begin
  @discrete_events begin
    (X < 1.0) & (t < 5.0) => [X ~ X + 2.0]
  end
  d, X --> 0
end

oprob = ODEProblem(decay_model, u0, tspan, ps)
sol = solve(oprob; dt = 0.001, adaptive = false)
plot(sol)
```

!!! note
    When we form composite boolean conditions for conditional discrete events, we use `&` to denote the AND operator (not `&&`, which is currently not supported). Similarly, `|` is used instead of `||`.

!!! warn
    Generally, discrete events including equality (`==`) will not be triggered. The reason is that the condition is only checked at the end of every time step. Hence, unless special precautions are taken to ensure that the simulator stops when the condition holds, the condition is unlikely to hold.


## [Creating event through callbacks](@id events_model_callbacks)

An alternative approach to events is *callbacks*. Like events, these combine a condition (for triggering the callback) with an affect (determining the callback's effect on the system once triggered). Both conditions and affects are functions, both taking an *integrator* as their input. Integrators are Julia structures which keep track of the solution's current state during the simulation (and can be used to e.g. check whether a species concentration has surpassed a threshold, or to update the value of a parameter). Integrators are described in more detail [here](https://docs.sciml.ai/DiffEqDocs/stable/basics/integrator/), with [this section](@ref ref) providing some useful details on how to interact with integrators generated from Catalyst models. Finally, a more throughout description of callbacks can be found [here](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/), while we below will provide a briefer introduction.

Unlike events, which are declared as a part of a model, callbacks are declared separately and then passed to the `solve` command (through the `callback` option). To use multiple callbacks, a [`CallbackSet` must be used](@ref events_model_callbacks_callbacksets). Before we can create any callbacks, we must activate the [DiffEqCallbacks.jl](https://github.com/SciML/DiffEqCallbacks.jl) package:
```@example events_callbacks
using DiffEqCallbacks
```

### [Preset-time callbacks](@id events_model_callbacks_presettime)
The simplest type of callbacks are preset-time callbacks (which are similar to [preset-time discrete events](@ref events_model_discrete_presettime)). Let us again consider a decay model. First, we create a normal `ODEProblem` for it:
```@example events_callbacks
using Catalyst
decay_model = @reaction_network begin
  d, X --> 0
end

u0 = [:X => 2.0]
tspan = (0.0, 10.0)
ps = [:d => 0.5]
oprob = ODEProblem(decay_model, u0, tspan, ps)
nothing # hide
```
Next, we create our callback. First, we create its condition and affect. For preset-time callbacks, the condition is a vector with all the time points at which we wish to trigger it. The affect is a function which takes a single argument (the simulation `integrator`), and updates it in some way (typically changing species or parameter values). E.g. here we create a callback which adds `2.0` units of `X` at time points `3.0` and `7.0`. 
```@example events_callbacks
cond = [3.0, 7.0] 
function affect!(integrator)
    integrator[:X] += 2.0
end
cb = PresetTimeCallback(cond, affect!)
```
How to index into integrators to check/update parameter/species values is described [here](@ref ref). Next, we simulate our model while supplying the callback to `solve` using the `callback` optional argument:

```@example events_callbacks
using OrdinaryDiffEq, Plots
sol = solve(oprob; callback = cb)
plot(sol)
```

As for events, it is possible to make the amount of `X` added when the callback is triggered a parameter. To do so we must declare a parameter `X_add` as part of our system, and then access its value in the `affect` function:
```@example events_callbacks
decay_model = @reaction_network begin
    @parameters X_add
    d, X --> 0
end

u0 = [:X => 2.0]
tspan = (0.0, 10.0)
ps = [:d => 0.5, :X_add => 3.0]
oprob = ODEProblem(decay_model, u0, tspan, ps)

cond = [3.0, 7.0] 
affect!(integrator) = (integrator[:X] += integrator.ps[:X_add])
cb = PresetTimeCallback(cond, affect!)
sol = solve(oprob; callback = cb)
plot(sol)
```

### [Discrete callbacks](@id events_model_callbacks_discrete)
Discrete callbacks use identical affect functions to preset-time callbacks. However, rather than their condition being a vector of values, it is a function which takes three inputs:
- `u`: The current state of the system (for Catalyst models a system's state is most easily accessed through the `integrator`, hence, the `u` input is best ignored).
- `t`: The current time point of the system
- `integrator` The simulation integrator (the same one which is modified by the affect function).

Its output is a [`Bool`](https://docs.julialang.org/en/v1/base/numbers/#Core.Bool). After each time step of the simulation, the callback's condition function is evaluated, if it returns `true`, the affect is triggered.

Here we simulate the same decay model, but with a discrete callback that adds `2.0` units of `X` whenever the concentration of `X` goes below `1.0`.
```@example events_callbacks
condition(u, t, integrator) = (integrator[:X] < 1.0)
affect!(integrator) = (integrator[:X] += 2.0)
cb = DiscreteCallback(condition, affect!)
sol = solve(oprob; callback = cb)
plot(sol)
```
Please note that (like for discrete events)(@ref events_model_discrete_conditional), discrete callbacks do not trigger at the exact time points when the condition holds (and only of the first time step when it holds). To achieve exact trigger times, a continuous callback is required.

### [Continuous callbacks](@id events_model_callbacks_continuous)
Finally, continuous callbacks (like continuous events) are different from discrete ones in that, rather than being checked after each simulation time step, the simulation finds the exact moment when the callback is triggered. Continuous callbacks utilise identical affects to discrete ones, and their conditions take identical input. However, their conditions return a number. This number should be non-negative, with the event being triggered whenever it evaluates to zero (and only when it transitions from positive to zero, transitions from negative to zero do not trigger continuous events).

Here, we can re-implement our previous discrete callback, but using a continuous one. This time our condition returns `integrator[:X] - 1.0`, which will be zero exactly when $X = 1.0$ (triggering our callback).
```@example events_callbacks
condition(u, t, integrator) = (integrator[:X] - 1.0)
affect!(integrator) = (integrator[:X] += 2.0)
cb = ContinuousCallback(condition, affect!)
```
We can now supply this callback to our `solve` command, yielding a similar simulation:
```@example events_callbacks
sol = solve(oprob; callback = cb)
plot(sol)
```

For situations where you have a large number of continuous callbacks, you should consider [implementing a `VectorContinuousCallback`](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/#VectorContinuousCallback), as this improves performance.

### [Creating `CallbackSet`s with multiple callbacks](@id events_model_callbacks_callbacksets)
The `callback` option in `solve` only accepts a single callback. To supply several callbacks, first [bundle these into a `CallbackSet`](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/#CallbackSet). Here we simulate our birth-death model, using both our preset-time callback and our discrete callback:
```@example events_callbacks
condition1 = [3.0, 7.0] 
affect1!(integrator) = (integrator[:X] += 2.0)
cb1 = PresetTimeCallback(condition1, affect1!)

condition2(u, t, integrator) = (integrator[:X] < 1.0)
affect2!(integrator) = (integrator[:X] += 2.0)
cb2 = DiscreteCallback(condition2, affect2!)

cbs = CallbackSet(cb1, cb2)
sol = solve(oprob; callback = cbs)
plot(sol)
```

### [Simulation termination and other callback-specific features](@id events_model_callbacks_special)
Callbacks give the user much more control over the affect and condition functions than events. They can e.g. be used to:
- Print a statement with information of the simulation state under certain conditions (using [`println`](https://docs.julialang.org/en/v1/base/io-network/#Base.println)).
- Reset a system's state to its value of the previous time step (using the `integrator.uprev` field).
- Terminate the simulation once some condition is reached (using the [`terminate!` function](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/#Example-2:-Terminating-an-Integration)).

Here we simulate our decay model, but terminates the solution once $X < 1.0$:
```@example events_callbacks
condition(u, t, integrator) = (integrator[:X] - 1.0)
affect!(integrator) = terminate!(integrator)
cb = ContinuousCallback(condition, affect!)
sol = solve(oprob; callback = cb)
plot(sol)
```

The DiffEqCallbacks package implements a range of ready-made callbacks that are available for the user. A full list of these can be found [here](https://docs.sciml.ai/DiffEqCallbacks/stable/).

## [Additional consideration when using events (or callbacks)](@id events_additional_considerations)

### [Events affecting changing model parameter values](@id events_additional_considerations_parameters)
It is not uncommon to have events (or callbacks) which change a system parameter. When doing so one has to be careful, since the event (or callback) will also update the original `ODEProblem`. E.g. in the following [birth-death process model](@ref ref), an event updates `p`'s value at time `5.0`:
```@example events_considerations
using Catalyst, OrdinaryDiffEq, Plots # hide
bd_model = @reaction_network begin
  @discrete_events begin
    [5.0] => [p ~ 2.0]
  end
  (p,d), 0 <--> X
end

u0 = [:X => 0.1]
tspan = (0.0, 10.0)
ps = [:p => 0.2, :d => 0.1]

oprob = ODEProblem(bd_model, u0, tspan, ps)
sol = solve(oprob)
plot(sol)
```
If we check the original `ODEProblem`, we can confirm that `ps`'s value has been updated:
```@example events_considerations
oprob.ps[:p]
```
Indeed, if we simulate our model again, the result is different:
```@example events_considerations
sol = solve(oprob)
plot(sol)
```

To deal with these situations, it is recommended to not pass your `ODEProblem` directly into `solve`, but rather a [`deepcopy`](https://docs.julialang.org/en/v1/base/base/#Base.deepcopy) of it:
```@example events_considerations
oprob = ODEProblem(bd_model, u0, tspan, ps)
sol = solve(deepcopy(oprob), Tsit5())
nothing # hide
```

### [Callbacks during jump simulations](@id events_additional_considerations_jump_parameters)

An assumption of (most) jump simulations is that the system's state is unchanged between reaction events. However, callbacks that change the system's state violate this assumption. To prevent erroneous simulations, users must inform the jump solver when the state has been updated this way. This allows the solver to reinitialize any internal state information that may have changed. This can be done through the `reset_aggregated_jumps!` function.

Here, we create a birth-death process where we add `2.0` additional units of `X` through a preset-time callback. At the end of the `affect!` function, we also call `reset_aggregated_jumps!` on the integrator:
```@example events_considerations
using DiffEqCallbacks # hide
using JumpProcesses
rn = @reaction_network begin
  (p,d), 0 <--> X
end

u0 = [:X => 1]
tspan = (0.0, 10.0)
ps = [:p => 1.0, :d => 0.2]
dprob = DiscreteProblem(bd_model, u0, tspan, ps)
jprob = JumpProblem(bd_model, dprob, Direct())

condition = [3.0, 7.0]
function affect!(integrator)
    integrator[:X] += 2.0
    reset_aggregated_jumps!(integrator)
end
cb = PresetTimeCallback(condition, affect!)

sol = solve(jprob, SSAStepper(); callback = cb)
plot(sol)
```

### [Events occurring at specific time points](@id events_additional_considerations_tstops)
Typically, events that occur at designated time points should be implemented through [preset-time events](@ref events_model_discrete_presettime) (or [preset-time callbacks](@ref events_model_callbacks_presettime)). However, sometimes one might want to use a a discrete event (or callback) instead. This is a potential problem, as discrete events (and callbacks) are only checked at the end of each time step. If you have defined a discrete event which is triggered when `t == 5.0`, but the simulation time stepper stops at $t = ... 4.9007, 5.0214, 5.1094, ...$, the condition won't hold, and is never triggered. Here, you have to explicitly instruct the simulation to make a stop at the specific time using the `tstops` argument.

In this example we wish to introduce an event which is triggered at time `5.0`, but only if `X > 2.0`. First, we create an `ODEProblem`:
```@example events_considerations
bd_model = @reaction_network begin
  @discrete_events begin
    (X > 2.0) & (t == 5.0) => [X ~ 2.0]
  end
  (p,d), 0 <--> X
end

u0 = [:X => 0.1]
tspan = (0.0, 10.0)
ps = [:p => 2.0, :d => 0.2]

oprob = ODEProblem(bd_model, u0, tspan, ps)
nothing # hide
```
If we now simulate the model normally, the event is not triggered:
```@example events_considerations
sol = solve(oprob)
plot(sol)
```
However, if we instruct `solve` to specifically stop at time `5.0`, it *is* triggered:
```@example events_considerations
sol = solve(oprob; tstops = [5.0])
plot(sol)
```
Here, `tstops` takes a vector with all the time points we wish to stop at (here we only designate `5.0`).

Since continuous events (and callbacks) actively find the exact time point when they are triggered, supplying `tstops` is not required for these.

## [Programmatic creation of models containing events](@ref events_programmatic)
We have previously described how [programmatic modelling](@ref ref) provides an alternative approach (to the DSL) for creating models in Catalyst. Here, we briefly describe how to add discrete and continuous events to models created programmatically.

In the DSL, events are created through the `continuous_events` and `discrete_events` options. In programmatic modelling, events are created using the same syntax, however:
- Any parameter or species (or [variables](@ref ref)) that occurs in an event must first be declared using the `@parameter` or `@species` (or `@variables`) options. If these are only accessed in a Callback, these must also [be explicitly added to the `ReactionSystem` upon creation](@ref ref).
- The events are stored in vectors (one for continuous events and one for discrete events).
- The events are supplied to the `ReactionSystem` constructor through the `continuous_events` and `discrete_events` optional arguments.

E.e. below we create a birth-death process programmatically, including both a continuous event, a preset-time discrete event, and a conditional discrete event. Since we have two different discrete events, we put both in the `discrete_events` vector (the `continuous_events` vector contains the single continuous events):
```@example events_programmatic
using Catalyst, OrdinaryDiffEq, Plots

t = default_t()
@species X(t)
@parameters p d X_add X_thres

rxs = [
  Reaction(p, [], [X]),
  Reaction(d, [X], []),
]

continuous_events = [
    [X ~ X_thres] => [X ~ X + X_add]
]
discrete_events = [
    [3.0, 7.0] => [d ~ 2*d],
    (X > 0.3) & (t > 2.0) => [X ~ 2.0]
]

@named bd_model = ReactionSystem(rxs, t; continuous_events, discrete_events)
bd_model = complete(bd_model)

u0 = [X => 0.1]
tspan = (0.0, 10.0)
ps = [p => 0.5, d => 0.2, X_add => 1.0, X_thres => 1.0]
oprob = ODEProblem(bd_model, u0, tspan, ps)

sol = solve(oprob)
plot(sol)
```