# [Modelling a periodic event during ODE and jump simulations](@id periodic_event_simulation_example)
This tutorial will describe how to simulate systems with periodic events in ODEs and jump simulations (SDEs use identical syntax). We will consider a model with a [circadian rhythm](https://en.wikipedia.org/wiki/Circadian_rhythm), where a parameter represents the level of light. While outdoor light varies smoothly, in experimental settings a lamp is often simply turned on/off every 12 hours. Here we will model this toggling of the light using a periodic event that is triggered every 12 hours.

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
using OrdinaryDiffEqTsit5, Plots
u0 = [:X => 150.0, :Xᴾ => 50.0]
ps = [:kₚ => 0.1, :kᵢ => 0.1, :l => 1.0]
tspan = (0.0, 100.0)
oprob = ODEProblem(circadian_model, u0, tspan, ps)
sol = solve(oprob, Tsit5())
plot(sol)
```

## [Modelling a circadian periodic event in a jump simulation](@id periodic_event_simulation_example_ode)
We can define periodic events in an identical manner for jump simulations. Let's
reuse our previously defined network, but now simulate it as a stochastic
chemical kinetics jump process model
```@example periodic_event_example
using JumpProcesses
u0 = [:X => 150, :Xᴾ => 50]  # define u0 as integers now
jinput = JumpInputs(circadian_model, u0, tspan, ps)
jprob = JumpProblem(jinput)
jsol = solve(jprob)
plot(jsol)
Catalyst.PNG(plot(jsol; fmt = :png, dpi = 200)) # hide
```

Sometimes we might want to setup a more complicated condition than simply that
our event occurs at a fixed time or with a fixed frequency. For example, suppose
we want to skip the first event at time `t = 12` and then have the event be
periodic after that point every 12 units of time. We can do so using a more
general discrete callback as follows
```@example periodic_event_example
circadian_model = @reaction_network begin
    @discrete_events begin
        ((t % 12 == 0) & (t > 12)) => [l ~ (l + 1)%2]
    end
    (kₚ*l,kᵢ), X <--> Xᴾ
end
```
Here our condition `((t % 12 == 0) & (t > 12))` determines when the event
occurs, evaluating to `true` when `t` is a multiple of `12` that is also larger
than `12`. We now finish specifying our model
```@example periodic_event_example
using JumpProcesses
u0 = [:X => 150, :Xᴾ => 50]
ps = [:kₚ => 0.1, :kᵢ => 0.1, :l => 1.0]
tspan = (0.0, 100.0)
jinput = JumpInputs(circadian_model, u0, tspan, ps)
jprob = JumpProblem(jinput)
nothing # hide
```
Next, if we simulate our model, we note that the events do not seem to be
triggered
```@example periodic_event_example
sol = solve(jprob)
plot(sol)
Catalyst.PNG(plot(sol; fmt = :png, dpi = 200)) # hide
```
The reason is that general discrete callbacks' conditions are only checked at
the end of each simulation time step, and the probability that these will
coincide with the callback's trigger times ($t = 12, 24, 36, ...$) is infinitely
small. Hence, we must directly instruct our simulation to stop at these time
points. We can do this using the `tstops` argument:
```@example periodic_event_example
tstops = range(12.0, tspan[2]; step = 12.0)
sol = solve(jprob; tstops)
plot(sol)
Catalyst.PNG(plot(sol; fmt = :png, dpi = 200)) # hide
```

## [Plotting the light level](@id periodic_event_simulation_plotting_light)
Sometimes when simulating models with periodic parameters, one would like to plot the parameter's value across the simulation. For this, there are two potential strategies. One includes creating a [*saving callback*](https://docs.sciml.ai/DiffEqCallbacks/stable/output_saving/#DiffEqCallbacks.SavingCallback). The other one, which we will demonstrate here, includes turning the parameter $l$ to a *variable* (so that its value is recorded throughout the simulation):
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
sol = solve(oprob, Tsit5())
nothing # hide
```
If we directly plot $l$'s value, it will be too small (compared to $X$ and $Xᴾ$ to be discernible). We instead [`@unpack` our variables](@ref dsl_advanced_options_symbolics_and_DSL_unpack), and then plot a re-scaled version:
```@example periodic_event_example
@unpack X, Xᴾ, l = circadian_model
plot(sol; idxs = [X, Xᴾ, 150*l], labels = ["X" "Xᴾ" "Light amplitude"])
```
