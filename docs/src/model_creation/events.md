# [Modelling Events](@id events)
In many applications one needs to model events that can occur when a set
condition is reached, such as providing a drug treatment at a specified time, or
turning off production of cells once the population reaches a given level.
Catalyst supports the event representation provided by ModelingToolkit, see
[here](https://docs.sciml.ai/ModelingToolkit/stable/basics/Events/), allowing
for both continuous and discrete events.

In this tutorial we'll illustrate how to make use of events. We will continue on the model 
of a cell with volume $V(t)$ that grows at a rate $\lambda$ that was considered in the
[coupled equations tutorial](@ref coupled_models).

## [Adding events](@id events_base_example)
Our current model is unrealistic in assuming the cell will grow exponentially
forever. Let's modify it such that the cell divides in half every time its
volume reaches a size of `2`. We also assume we lose half of the protein upon
division. Note, we will only keep track of one cell, and hence follow a specific
 lineage of the system. To do this we can create a continuous event using the
ModelingToolkit symbolic event interface and attach it to our system. Please see
the associated [ModelingToolkit
tutorial](https://docs.sciml.ai/ModelingToolkit/stable/basics/Events/) for more
details on the types of events that can be represented symbolically. A
lower-level approach for creating events via the DifferentialEquations.jl
callback interface is illustrated [here](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/) tutorial.

```@example events1
using Catalyst, OrdinaryDiffEqTsit5, Plots

rn = @reaction_network growing_cell begin
    # the growth rate
    @parameters λ = 1.0

    # assume there is no protein initially
    @species P(t) = 0.0

    # set the initial volume to 1.0
    @variables V(t) = 1.0

    # the reactions
    V, 0 --> P
    1.0, P --> 0

    # the coupled ODE for V(t)
    @equations begin
        D(V) ~ λ * V
    end

    # every 1.0 time unit we half the volume of the cell and the number of proteins
    @continuous_events begin 
        [V ~ 2.0] => [V ~ V/2, P ~ P/2]
    end
end
```
We can now create and simulate our model
```@example events1
oprob = ODEProblem(rn, [], (0.0, 10.0))
sol = solve(oprob, Tsit5())
plot(sol)
Catalyst.PNG(plot(sol; fmt = :png, dpi = 200)) # hide
```
We can also model discrete events. Here at a time `switch_time` we will set the parameter `k_on` to be
zero:
```@example events1
rn = @reaction_network param_off_ex begin
    @parameters switch_time
    k_on, A --> B
    k_off, B --> A
    
    @discrete_events begin
        (t == switch_time) => [k_on ~ 0.0]
    end
end

u0 = [:A => 10.0, :B => 0.0]
tspan = (0.0, 4.0)
p = [:k_on => 100.0, :switch_time => 2.0, :k_off => 10.0]
oprob = ODEProblem(rn, u0, tspan, p)
sol = solve(oprob, Tsit5(); tstops = 2.0)
plot(sol)
```
Note that for discrete events we need to set a stop time via `tstops` so that
the ODE solver can step exactly to the specific time of our event. In the
previous example we just manually set the numeric value of the parameter in the
`tstops` kwarg to `solve`, however, it can often be convenient to instead get
the value of the parameter from `oprob` and pass this numeric value. This helps
ensure consistency between the value passed via `p` and/or symbolic defaults and
what we pass as a `tstop` to `solve`. We can do this as
```@example events1
oprob = ODEProblem(rn, u0, tspan, p)
switch_time_val = oprob.ps[:switch_time]
sol = solve(oprob, Tsit5(); tstops = switch_time_val)
plot(sol)
```
For a detailed discussion on how to directly use the lower-level but more
flexible DifferentialEquations.jl event/callback interface, see the
[tutorial](https://docs.sciml.ai/Catalyst/stable/catalyst_applications/advanced_simulations/#Event-handling-using-callbacks)
on event handling using callbacks.

## [Adding events via the symbolic interface](@id events_symbolic)
Let's repeat the previous two models using the symbolic interface. We first
create our equations and unknowns/species again
```@example events2
using Catalyst, OrdinaryDiffEqTsit5, Plots
t = default_t()
D = default_time_deriv()

@parameters λ = 1.0
@variables V(t) = 1.0
@species P(t) = 0.0
eq = D(V) ~ λ * V
rx1 = @reaction $V, 0 --> $P
rx2 = @reaction 1.0, $P --> 0
```
Before creating our `ReactionSystem` we make the event.
```@example events2
# every 1.0 time unit we half the volume of the cell and the number of proteins
continuous_events = [V ~ 2.0] => [V ~ Pre(V)/2, P ~ Pre(P)/2]
```
We can now create and simulate our model
```@example events2
@named rs = ReactionSystem([rx1, rx2, eq], t; continuous_events)
rs = complete(rs)

oprob = ODEProblem(rs, [], (0.0, 10.0))
sol = solve(oprob, Tsit5())
plot(sol)
Catalyst.PNG(plot(sol; fmt = :png, dpi = 200)) # hide
```
We can again also model discrete events. Similar to our example with continuous
events, we start by creating reaction equations, parameters, variables, and
unknowns.
```@example events2
t = default_t()
@parameters switch_time k_off
@discretes k_on(t)
@species A(t) B(t)

rxs = [(@reaction $k_on, A --> B), (@reaction k_off, B --> A)]
```
Now we add an event such that at time `t` (`switch_time`), `k_on` is set to zero.
```@example events2
discrete_events = ModelingToolkitBase.SymbolicDiscreteCallback((t == switch_time) => [k_on ~ 0.0]; discrete_parameters = [k_on])

u0 = [:A => 10.0, :B => 0.0]
tspan = (0.0, 4.0)
p = [k_on => 100.0, switch_time => 2.0, k_off => 10.0]
```
Simulating our model,
```@example events2
@named rs2 = ReactionSystem(rxs, t, [A, B], [k_on, k_off, switch_time]; discrete_events)
rs2 = complete(rs2)

oprob = ODEProblem(rs2, u0, tspan, p)
switch_time_val = oprob.ps[:switch_time]
sol = solve(oprob, Tsit5(); tstops = switch_time_val)
plot(sol)
```
