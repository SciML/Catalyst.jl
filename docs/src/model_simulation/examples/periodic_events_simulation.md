# [Modelling a periodic event during ODE and jump simulations](@id periodic_event_simulation_example)
This tutorial will describe how to simulate systems with periodic events in ODE (straightforward) and jump (slightly less, but still fairly, straightforward) simulations (SDEs use identical syntax to ODEs). We will consider a model with a [circadian rhythm](https://en.wikipedia.org/wiki/Circadian_rhythm), where a parameter represents the level of light. While outdoor light varies smoothly, in experimental settings a lamp is often simply turned on/off every 12 hours. Here we will model this toggling of the light using a periodic event that is triggered every 12 hours.

## [Modelling a circadian periodic event in an ODE simulation](@id periodic_event_simulation_example_ode)
We will consider a simple circadian model, consisting of a single protein ($X$), which is phosphorylated ($X \to Xᴾ$) in the presence of light ($l$). Here, the light parameter can either be $0$ (night) or $1$ (day). We can model this using a simple periodic event which switches the value of $l$ every 12 hours (here, `%` is the [remainder operator](https://docs.julialang.org/en/v1/manual/mathematical-operations/#Arithmetic-Operators)).
```@example periodic_event_example
using Catalyst
circadian_model = @reaction_network begin
    @discrete_events begin
        12 => [l ~ (l + 1)%2]
    end
 (kₚ*l,kᵢ), X <--> Xᴾ
end
```
We can now simulate this model, observing how a 24-hour cycle is reached
```@example periodic_event_example
using OrdinaryDiffEq, Plots
u0 = [:X => 150.0, :Xᴾ => 50.0]
ps = [:kₚ => 0.1, :kᵢ => 0.1, :l => 1.0]
oprob = ODEProblem(circadian_model, u0, (0.0, 100.0), ps)
sol = solve(oprob)
plot(sol)
```

## [Modelling a circadian periodic event in a jump simulation](@id periodic_event_simulation_example_ode)
While defining periodic events is easy for ODE and SDE simulations, due to how events are internally implemented, these cannot currently be used for jump simulations. Instead, we create a new model (and its corresponding `JumpProblem`), without the event:
```@example periodic_event_example
circadian_model = @reaction_network begin
 (kₚ*l,kᵢ), X <--> Xᴾ
end

using JumpProcesses
u0 = [:X => 150, :Xᴾ => 50]
ps = [:kₚ => 0.1, :kᵢ => 0.1, :l => 1.0]
dprob = DiscreteProblem(circadian_model, u0, (0.0, 100.0), ps)
jprob = JumpProblem(circadian_model, dprob, Direct())
nothing # hide
```
Next, we implement our event through a `DiscreteCallback`. It triggers every 12 time units, and when it does, it flips the value of `l`. Since we will interface with $l$'s value multiple times, we [create specific `getl` and `setl` functions](@ref simulation_structure_interfacing_functions) to do this (in practice, due to this model's small size, the performance increase is minimal).
```@example periodic_event_example
getl = ModelingToolkit.getp(jprob, :l)
setl = ModelingToolkit.setp(jprob, :l)
function condition(u, t, integrator)
 (t % 12) == 0.0
end
function affect!(integrator)
    setl(integrator, (getl(integrator) + 1) % 2)
    reset_aggregated_jumps!(integrator)
end
callback = DiscreteCallback(condition, affect!)
nothing # hide
```
!!! danger
    Whenever the integrator of a jump simulator is modified, the `reset_aggregated_jumps!` function must be called (with the integrator as its input) afterwards. Omitting to do so will result in simulations (incorrectly) using the old values of parameters/states internally instead of the updated values.

Next, if we simulate our model, we note that the events do not seem to be triggered
```@example periodic_event_example
sol = solve(jprob, SSAStepper(); callback)
plot(sol)
```
The reason is that discrete callbacks' conditions are only checked at the end of each simulation time step, and the probability that these will coincide with the callback's trigger times ($t = 12, 24, 36, ...$) is infinitely small. Hence, we must directly instruct our simulation to stop at these time points. We can do this using the `tstops` argument:
```@example periodic_event_example
tstops = 12.0:12.0:dprob.tspan[2]
sol = solve(jprob, SSAStepper(); callback, tstops)
plot(sol)
```

## [Plotting the light level](@id periodic_event_simulation_plotting_light)
Sometimes when simulating models with periodic parameters, one would like to plot the parameter's value across the simulation. To do this, we must turn our parameter $l$ to a *variable* (so that its value is recorded throughout the simulation):
```@example periodic_event_example
circadian_model = @reaction_network begin
    @variables l(t)
    @discrete_events begin
        12 => [l ~ (l + 1)%2]
    end
 (kₚ*l,kᵢ), X <--> Xᴾ
end
nothing # hide
```
Next, we simulate our model like before (but providing $l$'s value as an initial condition):
```@example periodic_event_example
u0 = [:X => 150.0, :Xᴾ => 50.0, :l => 1.0]
ps = [:kₚ => 0.1, :kᵢ => 0.1]
oprob = ODEProblem(circadian_model, u0, (0.0, 100.0), ps)
sol = solve(oprob)
nothing # hide
```
If we directly plot $l$'s value, it will be too small (compared to $X$ and $Xᴾ$ to be discernible). We instead [`@unpack` our variables](@ref dsl_advanced_options_symbolics_and_DSL_unpack), and then plot a re-scaled version:
```@example periodic_event_example
@unpack X, Xᴾ, l = circadian_model
plot(sol; idxs = [X, Xᴾ, 150*l], labels = ["X" "Xᴾ" "Light amplitude"])
```

!!! note
    If you wish to reproduce this in a jump simulation, remember to make appropriate modifications (like using `setu` instead of `setp`).