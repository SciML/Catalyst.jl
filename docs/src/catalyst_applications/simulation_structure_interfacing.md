# [Interfacing problems, integrators, and solutions](@id simulation_structure_interfacing)
When simulating a model, one begins with creating a [problem](https://docs.sciml.ai/DiffEqDocs/stable/basics/problem/). Next, a simulation is performed on a problem, during which the state of the simulation is recorded through an [integrator](https://docs.sciml.ai/DiffEqDocs/stable/basics/integrator/). Finally, the simulation output is returned as a [solution](https://docs.sciml.ai/DiffEqDocs/stable/basics/solution/). This tutorial describes how to access, or modify the state, or parameter, values of problems, integrators, and solutions structures.

Generally, when we have a structure `simulation_struct` and want to interface with the state (or parameter) `G`, we use `simulation_struct[:G]` to access the value, and `simulation_struct[:G] = 5.0` to set it to a new value. However, see the following examples for full details.

## Interfacing problem objects

We begin by demonstrating how we can interface with problem objects. We will demonstrate using a `ODEProblem`, however, it works similarly for other problem types.
```@example ex1
using Catalyst
rn = @reaction_network begin
    (k1,k2), X1 <--> X2
end
```

To later enable us to index solutions and problems using symbolic variables we
also need to indicate that we are finished constructing our system using `complete`:
```@example ex1
rn = complete(rn)
nothing  # hide
```

We can then finish creating our `ODEProblem`
```@example ex1
u0 = [:X1 => 1.0, :X2 => 5.0]
p = [:k1 => 5.0, :k2 => 2.0]
tspan = (0.0, 10.0)
oprob = ODEProblem(rn, u0, tspan, p)
nothing    # hide
```

Note that we could have alternatively specified the `u0` and `p` mappings using symbolic variables like
```@julia
@unpack X1, X2, k1, k2 = rn
u0 = [X1 => 1.0, X2 => 5.0]
p = [k1 => 5.0, k2 => 2.0]
```
or if we did not want to call `@unpack`, and have to list the variables, by
```@example ex1
u0 = [rn.X1 => 1.0, rn.X2 => 5.0]
p = [rn.k1 => 5.0, rn.k2 => 2.0]
```
Note that this last form of construction requires that we called `rn =
complete(rn)` earlier.

Given our problem, `oprob`, we can access the value of the initial condition via
```@example ex1
oprob.u0
```
and the parameter values via
```@example ex1
oprob.p
```
Note that these are ordered consistent to `states(rn)` and `parameters(rn)`
respectively. Individual species can be indexed like

```@example ex1
# load the symbolic variable
@unpack X1 = rn
oprob[X1]
```
Alternatively, we could access this value via
```@example ex1```
# these are all the same
oprob[X1] == oprob[rn.X1] == oprob[:X1]
```

To ensure type stability, parameters must be accessed using the interface provided by [SymbolicIndexingInterface.jl](https://docs.sciml.ai/SymbolicIndexingInterface/stable/). For example, we can create a function to return the value of `k1` and then evaluate it via

```@example ex1
# `getp` is used to make a function to return parameter values
using ModelingToolkit: getp

@unpack k1, k2 = rn

# this creates a function that returns k1
getk1 = getp(rn, k1)
getk1(oprob)
```
Equivalent to this would be
```@example ex1
g2 = getp(rn, rn.k1)
g3 = getp(rn, :k1)
getk1(oprob) == g2(oprob) == g3(oprob)
```

Finally, the time span to solve over is given by
```@example ex1
oprob.tspan
```

#### Remaking problems using the `remake` function
To modify these fields, it is recommended to use the `remake` function. `remake` creates a new problem object, taking an existing problem as input along with a mapping that indicates any fields you wish to modify (and their new values). Thus, we can do:
```@example ex1
using DifferentialEquations
@unpack X1, X2, k1, k2 = rn
oprob2 = remake(oprob; u0 = [X1 => 10.0, X2 => 50.0], tspan = (0.0, 100.0),
                        p = [k1 => 50.0, k2 => 20.0])
nothing    # hide
```
and we can now check the fields of `oprob2`
```@example ex1
oprob2.u0
```
```@example ex1
oprob2.tspan
```
```@example ex1
oprob2.p
```

When using `remake`, we only have to provide the fields that we actually wish to change, e.g.
```@example ex1
oprob3 = remake(oprob; u0 = [X1 => 10.0])
nothing    # hide
```
will only update `X1`'s value in the initial conditions.

Note that if we didn't want to use `@unpack`, we could instead access the symbolic variables and parameters via the `ReactionSystem`, `rn`, i.e. like
```@example ex1
oprob2 = remake(oprob; u0 = [rn.X1 => 10.0, rn.X2 => 50.0], tspan = (0.0, 100.0),
                        p = [rn.k1 => 50.0, rn.k2 => 20.0])
```
giving

```@example ex1
oprob2.u0
```
```@example ex1
oprob2.tspan
```
```@example ex1
oprob2.p
```

Note that using `Symbol`s like `:k1` for parameters is not currently supported with `remake`.

## Interfacing integrator objects

During a simulation, the solution is stored in an integrator object, we will here describe how to interface with these. The primary circumstance under which a user may wish to do so is when using [callbacks](@ref advanced_simulations_callbacks). We can create an integrator by calling `init` on our problem ([while circumstances where the user might want to use `init` function exist](https://docs.sciml.ai/DiffEqDocs/stable/basics/integrator/#Initialization-and-Stepping), since integrators are automatically created during simulations, these are rare):
```@example ex1
integrator = init(oprob)
```
As for problems, `integrator.p` stores the vector of numeric parameter values.
During an ODE simulation, `integrator.u` stores the current solution value, i.e.
the vector of `[X1(integrator.t), X2(integrator.t)]`.

Using a similar syntax to problems, we can get the current values of a state within the integrator via symbolic indexing
```@example ex1
integrator[X1]
```
or a parameter:
```@example ex1
getk1(integrator)
```
Similarly, we can update the current values of state variables via
```@example ex1
integrator[X1] = 10.0
```
This is equivalent to manually setting `integrator.u[1] = 10.0` (since `X1` is the first entry in `states(rn)`).

Updating parameters requires create a function for setting their value, similar to getting parameter values:
```@example ex1
using ModelingToolkit: setp
setk1 = setp(rn, k1)
setk1(integrator, 2.5)
getk1(integrator) == 2.5
```

Note that parameter getting and setting functions can be built to work for multiple parameters by passing a vector or tuple (depending on the type of `oprob.p`):
```@example ex1
getpvec = getp(rn, [k1, k2])
getpvec(integrator)
```
and
```@example ex1
setpvec = setp(rn, [k1, k2])
setpvec(integrator, [1.0, 9.0])
integrator.p == [1.0, 9.0]
```

As for problems, we can also use `rn.X1` or `:X1` to access state values and `rn.k1` or `:k1` when accessing parameter values.

*Note, when updating integrators for Gillespie type simulations, i.e. using `JumpProblem`s, users need to call a function to tell the solver if a parameter has been updated (as this requires changing internal solver state too). Please see [this](@ref advanced_simulations_ssa_callbacks) for details.*
