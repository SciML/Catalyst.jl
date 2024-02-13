# [Constraint Equations and Events](@id constraint_equations)
In many applications one has additional algebraic or differential equations for
non-chemical species that can be coupled to a chemical reaction network model.
Catalyst supports coupled differential and algebraic equations, and currently
allows conversion of such coupled systems to ModelingToolkit `ODESystem`s and
`NonlinearSystem`s. Likewise, one often needs events that can occur when a set
condition is reached, such as providing a drug treatment at a specified time, or
turning off production of cells once the population reaches a given level.
Catalyst supports the event representation provided by ModelingToolkit, see
[here](https://docs.sciml.ai/ModelingToolkit/stable/basics/Events/), allowing
for both continuous and discrete events (though only the latter are supported
when converting to `JumpSystem`s currently).

In this tutorial we'll illustrate how to make use of constraint equations and
events. Let's consider a model of a cell with volume $V(t)$ that grows at a rate
$\lambda$. For now we'll assume the cell can grow indefinitely. We'll also keep
track of one protein $P(t)$, which is produced at a rate proportional $V$ and
can be degraded.

## [Coupling ODE constraints via extending a system](@id constraint_equations_coupling_constratins)

There are several ways we can create our Catalyst model with the two reactions
and ODE for $V(t)$. One approach is to use compositional modeling, create
separate `ReactionSystem`s and `ODESystem`s with their respective components,
and then extend the `ReactionSystem` with the `ODESystem`. Let's begin by
creating these two systems. 

Here, to create differentials with respect to time (for our differential equations), we must import the time differential operator from Catalyst. We do this through `import Catalyst: D_nounits as D` (calling it `D`). Here, `D(V)` denotes the differential of the variable `V` with respect to time.

```@example ceq1
using Catalyst, DifferentialEquations, Plots
import Catalyst: t_nounits as t, D_nounits as D

# set the growth rate to 1.0
@parameters λ = 1.0

# set the initial volume to 1.0
@variables V(t) = 1.0

# build the ODESystem for dV/dt
eq = [D(V) ~ λ * V]
@named osys = ODESystem(eq, t)

# build the ReactionSystem with no protein initially
rn = @reaction_network begin
    @species P(t) = 0.0
    $V,   0 --> P
    1.0, P --> 0
end
```
Notice, here we interpolated the variable `V` with `$V` to ensure we use the
same symbolic unknown variable in the `rn` as we used in building `osys`. See the
doc section on [interpolation of variables](@ref
dsl_description_interpolation_of_variables) for more information.

We can now merge the two systems into one complete `ReactionSystem` model using
[`ModelingToolkit.extend`](@ref):
```@example ceq1
@named growing_cell = extend(osys, rn)
```

We see that the combined model now has both the reactions and ODEs as its
equations. To solve and plot the model we proceed like normal
```@example ceq1
oprob = ODEProblem(growing_cell, [], (0.0, 1.0))
sol = solve(oprob, Tsit5())
plot(sol)
```

## Coupling ODE constraints via directly building a `ReactionSystem`
As an alternative to the previous approach, we could have constructed our
`ReactionSystem` all at once by directly using the symbolic interface:
```@example ceq2
using Catalyst, DifferentialEquations, Plots
import Catalyst: t_nounits as t, D_nounits as D

@parameters λ = 1.0
@variables V(t) = 1.0
eq = D(V) ~ λ * V
rx1 = @reaction $V, 0 --> P
rx2 = @reaction 1.0, P --> 0
@named growing_cell = ReactionSystem([rx1, rx2, eq], t)
setdefaults!(growing_cell, [:P => 0.0])

oprob = ODEProblem(growing_cell, [], (0.0, 1.0))
sol = solve(oprob, Tsit5())
plot(sol)
```

## [Adding events](@id constraint_equations_events)
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
callback interface is illustrated in the [Advanced Simulation Options](@ref
advanced_simulations) tutorial.

Let's first create our equations and unknowns/species again
```@example ceq3
using Catalyst, DifferentialEquations, Plots
import Catalyst: t_nounits as t, D_nounits as D

@parameters λ = 1.0
@variables V(t) = 1.0
@species P(t) = 0.0
eq = D(V) ~ λ * V
rx1 = @reaction $V, 0 --> $P
rx2 = @reaction 1.0, $P --> 0
```
Before creating our `ReactionSystem` we make the event.
```@example ceq3
# every 1.0 time unit we half the volume of the cell and the number of proteins
continuous_events = [V ~ 2.0] => [V ~ V/2, P ~ P/2]
```
We can now create and simulate our model
```@example ceq3
@named rs = ReactionSystem([rx1, rx2, eq], t; continuous_events)

oprob = ODEProblem(rs, [], (0.0, 10.0))
sol = solve(oprob, Tsit5())
plot(sol; plotdensity = 1000)
```
We can also model discrete events. Similar to our example with continuous events, we start by creating reaction equations, parameters, variables, and unknowns. 
```@example ceq3
import Catalyst: t_nounits as t
@parameters k_on switch_time k_off
@species A(t) B(t)

rxs = [(@reaction k_on, A --> B), (@reaction k_off, B --> A)]
```
Now we add an event such that at time `t` (`switch_time`), `k_on` is set to zero. 
```@example ceq3
discrete_events = (t == switch_time) => [k_on ~ 0.0]

u0 = [:A => 10.0, :B => 0.0]
tspan = (0.0, 4.0)
p = [k_on => 100.0, switch_time => 2.0, k_off => 10.0]
```
Simulating our model, 
```@example ceq3
@named osys = ReactionSystem(rxs, t, [A, B], [k_on, k_off, switch_time]; discrete_events)

oprob = ODEProblem(osys, u0, tspan, p)
sol = solve(oprob, Tsit5(); tstops = 2.0)
plot(sol)
```
Note that for discrete events we need to set a stop time, `tstops`, so that the ODE solver can step exactly to the specific time of our event. For a detailed discussion on how to directly use the lower-level but more flexible DifferentialEquations.jl event/callback interface, see the [tutorial](https://docs.sciml.ai/Catalyst/stable/catalyst_applications/advanced_simulations/#Event-handling-using-callbacks) on event handling using callbacks. 