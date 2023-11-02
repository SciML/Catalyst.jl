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

u0 = [:X1 => 1.0, :X2 => 5.0]
p = [:k1 => 5.0, :k2 => 2.0]
oprob = ODEProblem(rn, u0, (0.0,10.0), p)
nothing    # hide
```

We can find the value of a state simply by interfacing with the corresponding symbol:
```@example ex1
oprob[:X1]
```
with the notation being identical for parameters:
```@example ex1
oprob[:k1]
```

If we want to change a state's initial condition value, we use the following notation
```@example ex1
oprob[:X1] = 10.0
```
with parameters using the same notation.

#### [Remaking problems using the `remake` function](@od simulation_structure_interfacing_remake)
Typically, when modifying problems, it is recommended to use the `remake` function. Unlike when we do `oprob[:X1] = 10.0` (which modifies the problem in question), `remake` creates a new problem object. The `remake` function takes a problem as input, and any fields you wish to modify (and their new values) as optional inputs. Thus, we can do:
```@example ex1
using DifferentialEquations
@unpack X1, X2, k1, k2 = rn
oprob1 = ODEProblem(rn, u0, (0.0,10.0), p)
oprob2 = remake(oprob1; u0=[X1 => 10.0, X2 => 50.0], tspan=(0.0,100.0), p=[k1 => 50.0,k2 => 20.0])
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
Please note that, currently, `remake` does not work while giving `Symbol`s as input (e.g `[:X1 => 10.0, :X2 => 50.0]`), but we need to unpack the symbolic variables and use them instead (please see the end of this tutorial for more information on using symbolic variables rather than `Symbol`s).

When using `remake`, we only have to provide the fields that we actually wish to change, e.g.
```@example ex1
oprob3 = remake(oprob1; u0=[X1 => 10.0, X2 => 50.0])
nothing    # hide
```
will only update the initial conditions.


## Interfacing integrator objects

During a simulation, the solution is stored in an integrator object, we will here describe how to interface with these. The primary circumstance under which a user may wish to do so is when using [callbacks](@ref advanced_simulations_callbacks). We can create an integrator by calling `init` on our problem ([while circumstances where the user might want to use `init` function exist](https://docs.sciml.ai/DiffEqDocs/stable/basics/integrator/#Initialization-and-Stepping), since integrators are automatically created during simulations, these are rare):
```@example ex1
integrator = init(oprob)
```
Using a similar syntax to problems, we can get the current values of a state within the integrator:
```@example ex1
integrator[:X1]
```
or a parameter:
```@example ex1
integrator[:k1]
```
Similarly, we can update their values using:
```@example ex1
integrator[:X1] = 10.0
```
Please read [this](@ref advanced_simulations_ssa_callbacks) with regards to updating integrators of `JumpProblem`s.


## [Interfacing solution objects](@id simulation_structure_interfacing_solutions)

Finally, we consider solution objects. First, we simulate our problem:
```@example ex1
sol = solve(oprob)
```
For solutions, when we access a state, we get its whole simulation vector:
```@example ex1
sol[:X1]
```
while when we access a parameter we only get a single value:
```@example ex1
sol[:k1]
```
Finally, we note that we cannot change the values of solution states or parameters (i.e. both `sol[:X1] = 0.0` and `sol[:k1] = 0.0` generate errors).

## [Interfacing using symbolic representation](@id simulation_structure_interfacing_symbolic_representation)

Catalyst is built on an *intermediary representation* implemented by (ModelingToolkit.jl)[https://github.com/SciML/ModelingToolkit.jl]. ModelingToolkit is a modelling framework where one first declares a set of symbolic variables and parameters using e.g.
```@example ex2
using ModelingToolkit
@parameters σ ρ β
@variables t x(t) y(t) z(t)
nothing # hide
```
and then uses these to build systems of equations. Here, these symbolic variables (`x`, `y`, and `z`) and parameters (`σ`, `ρ`, and `β`) can be used to interface a `problem`, `integrator`, and `solution` object (like we did previously, but using Symbols, e.g. `:X`). Since Catalyst models are built on ModelingToolkit, these models also contain similar symbolic variables and parameters.
```@example ex2
using Catalyst
rn = @reaction_network begin
    (k1,k2), X1 <--> X2
end

@unpack k1,k2,X1,X2 = rn
```
Here, we first list the parameters and variables (for reaction systems the latter are typically species) we wish to import (in this case we select all, but we could select only a subset), next we denote from which model (here `rn`) from which we wish to import from. Next, these values can be used directly to interface with e.g. an `ODEProblem`:
```@example ex2
u0 = [X1 => 1.0, X2 => 5.0]
p = [:k1 => 5.0, :k2 => 2.0]
oprob = ODEProblem(rn, u0, (0.0,10.0), p)

oprob[k1]
```
To interface with integrators and solutions we use a similar syntax.

Finally, instead of using `@unpack` to access a symbolic variable or parameter, we can access it directly using `rn.X1`, and thus access a state of our `ODEProblem` using
```@example ex2
oprob[rn.X1]
```
