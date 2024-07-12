# [Interfacing problems, integrators, and solutions](@id simulation_structure_interfacing)
When simulating a model, one begins with creating a [problem](https://docs.sciml.ai/DiffEqDocs/stable/basics/problem/). Next, a simulation is performed on the problem, during which the simulation's state is recorded through an [integrator](https://docs.sciml.ai/DiffEqDocs/stable/basics/integrator/). Finally, the simulation output is returned as a [solution](https://docs.sciml.ai/DiffEqDocs/stable/basics/solution/). This tutorial describes how to access (or modify) the state (or parameter) values of problem, integrator, and solution structures.

Generally, when we have a structure `simulation_struct` and want to interface with the unknown (or parameter) `x`, we use `simulation_struct[:x]` to access the value. For situations where a value is accessed (or changed) a large number of times, it can *improve performance* to first create a [specialised getter/setter function](@ref simulation_structure_interfacing_functions).

## [Interfacing problem objects](@id simulation_structure_interfacing_problems)

We begin by demonstrating how we can interface with problem objects. First, we create an `ODEProblem` representation of a [chemical cross-coupling model](@ref basic_CRN_library_cc) (where a catalyst, $C$, couples two substrates, $S₁$ and $S₂$, to form a product, $P$).
```@example structure_indexing
using Catalyst
cc_system = @reaction_network begin
    k₁, S₁ + C --> S₁C
    k₂, S₁C + S₂ --> CP
    k₃, CP --> C + P
end

u0 = [:S₁ => 1.0, :C => 0.05, :S₂ => 1.2, :S₁C => 0.0, :CP => 0.0, :P => 0.0]
tspan = (0., 10.0)
ps = [:k₁ => 5.0, :k₂ => 5.0, :k₃ => 100.0] 
oprob = ODEProblem(cc_system, u0, tspan, ps)
nothing    # hide
```

We can find a species's (or [variable's](@ref constraint_equations)) initial condition value by simply indexing with the species of interest as input. Here we check the initial condition value of $C$:
```@example structure_indexing
oprob[:C]
```
An almost identical notation can be used for parameters, however, here we use `oprob.ps` (rather than `oprob`):
with the notation being identical for parameters:
```@example structure_indexing
oprob.ps[:k₁]
```
To retrieve several species initial condition (or parameter) values, simply give a vector input. Here we check the values of the two substrates ($S₁$ and $S₂$):
```@example structure_indexing
oprob[[:S₁, :S₂]]
```

A problem's time span can be accessed through the `tspan` field:
```@example structure_indexing
oprob.tspan
```

!!! note
    Here we have used an `ODEProblem`to demonstrate all interfacing functionality. However, identical workflows work for the other problem types.

### [Remaking problems using the `remake` function](@id simulation_structure_interfacing_problems_remake)
To modify a problem, the `remake` function should be used. It takes an already created problem, and returns a new, updated, one (the input problem is unchanged). The `remake` function takes the following inputs:
- The problem that it remakes.
- (optionally) `u0`: A vector with initial conditions that should be updated. The vector takes the same form as normal initial condition vectors, but does not need to be complete (in which case only a subset of the initial conditions are updated).
- (optionally) `tspan`: An updated time span (using the same format as time spans normally are given in).
- (optionally) `p`: A vector with parameters that should be updated. The vector takes the same form as normal parameter vectors, but does not need to be complete (in which case only a subset of the parameters are updated).

Here we modify our problem to increase the initial condition concentrations of the two substrates ($S₁$ and $S₂$), and also confirm that the new problem is different from the old (unchanged) one:
```@example structure_indexing
using OrdinaryDiffEq
oprob_new = remake(oprob; u0 = [:S₁ => 5.0, :S₂ => 2.5])
oprob_new != oprob
```
Here, we instead use `remake` to simultaneously update all three fields:
```@example structure_indexing
oprob_new_2 = remake(oprob; u0 = [:C => 0.2], tspan = (0.0, 20.0), p = [:k₁ => 2.0, :k₂ => 2.0])
nothing # hide
```

Typically, when using `remake` to update a problem, the common workflow is to overwrite the old one with the output. E.g. to set the value of `k₁` to `5.0` in `oprob`, you would do:
```@example structure_indexing
oprob = remake(oprob; p = [:k₁ => 5.0])
nothing # hide
```

## [Interfacing integrator objects](@id simulation_structure_interfacing_integrators)

During a simulation, the solution is stored in an integrator object. Here, we will describe how to interface with these. The almost exclusive circumstance when integrator-interfacing is relevant is when simulation events are implemented through [callbacks](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/). However, to demonstrate integrator indexing in this tutorial, we will create one through the `init` function (while circumstances where one might [want to use `init` function exist](https://docs.sciml.ai/DiffEqDocs/stable/basics/integrator/#Initialization-and-Stepping), since integrators are automatically created during simulations, these are rare).
```@example structure_indexing
integrator = init(oprob)
nothing # hide
```
We can interface with our integrator using an identical syntax as [was used for problems](@ref simulation_structure_interfacing_problems). The primary exception is that there is no `remake` function for integrators. Instead, we can update species and parameter values using normal indexing. Here we update, and then check the values of, first the species $C$ and then the parameter $k₁$:
```@example structure_indexing
integrator[:C] = 0.0
integrator[:C]
```
or a parameter:
```@example structure_indexing
integrator.ps[:k₂] = 1.0
integrator.ps[:k₂]
```
Note that here, species-interfacing yields (or changes) a simulation's current value for a species, not its initial condition.

If you are interfacing with jump simulation integrators, you must always call `reset_aggregated_jumps!(integrator)` afterwards.

## [Interfacing solution objects](@id simulation_structure_interfacing_solutions)

Finally, we consider solution objects. First, we simulate our problem:
```@example structure_indexing
sol = solve(oprob)
nothing # hide
```
Next, we can access the simulation's values using the same notation as previously. When we access a species's, its values across the full simulation is returned as a vector:
```@example structure_indexing
sol[:P]
```
Parameter values can also be accessed (however, here we only get a single value):
```@example structure_indexing
sol.ps[:k₃]
```
Unlike for problems and integrators, species or parameter values of solutions cannot be changed. 

A vector with the time values for all simulation time steps can be retrieved using
```@example structure_indexing
sol.t
```

To find simulation values at a specific time point, simply use this time point as input to your solution object (treating it as a function). I.e. here we get our simulation's values at time $t = 1.0$
```@example structure_indexing
sol(1.0)
```
This works whenever the simulations actually stopped at time $t = 1.0$ (if not, an interpolated value is returned). To get the simulation's values for a specific subset of species, we can use the `idxs` optional argument. I.e. here we get the value of $C$ at time $t = 1.0$
```@example structure_indexing
sol(1.0; idxs = [:C])
``` 

## [Interfacing using specialised getter/setter functions](@id simulation_structure_interfacing_functions)
Internally, species and parameter values are stored in vectors. Whenever e.g. `oprob[:C]` is called, Julia must first find which index in the storage vector $C$ is stored in. Next, its value can be retrieved. If `oprob[:C]` is called a large number of times, this index must be found in each call. If a large number of such accesses are carried out, and performance is essential, it can be worthwhile to pre-compute a function to carry this out.

There exist four different functions, each returning a function for performing a specific type of interfacing:
- `ModelingToolkit.getu`: For accessing species values.
- `ModelingToolkit.getp`: For accessing parameter values.
- `ModelingToolkit.setu`: For changing species values.
- `ModelingToolkit.setp`: For changing parameter values.

For each species (or parameter) we wish to interface with, a new interfacing function must be created. Here we first creates a function for retrieving the value of $C$, and then use it for this purpose:
```@example structure_indexing
get_C = ModelingToolkit.getu(oprob, :C)
get_C(oprob)
```
Here, `getu` (as well as `getp`, `setu`, and `setp`) first takes the structure we wish to interface with, and then the target quantity. When using `setu` and `setp`, in the second step, we must also provide the update value:
```@example structure_indexing
set_C = ModelingToolkit.setu(oprob, :C)
set_C(oprob, 0.2)
get_C(oprob)
```

Like when indexing-based interfacing is used, these functions also work with vectors:
```@example structure_indexing
get_S = ModelingToolkit.getu(oprob, [:S₁, :S₂])
get_S(oprob)
```

## [Interfacing using symbolic representations](@id simulation_structure_interfacing_symbolic_representation)
When e.g. [programmatic modelling is used](@ref programmatic_CRN_construction), species and parameters can be represented as *symbolic variables*. These can be used to index a problem, just like symbol-based representations can. Here we create a simple [two-state model](@ref basic_CRN_library_two_states) programmatically, and use its symbolic variables to check, and update, an initial condition:
```@example structure_indexing_symbolic_variables
using Catalyst, OrdinaryDiffEq
t = default_t()
@species X1(t) X2(t)
@parameters k1 k2
rxs = [
    Reaction(k1, [X1], [X2]),
    Reaction(k2, [X2], [X1])
]
@named two_state_model = ReactionSystem(rxs, t)
two_state_model = complete(two_state_model)

u0 = [X1 => 2.0, X2 => 0.0]
tspan = (0.0, 1.0)
ps = [k1 => 1.0, k2 => 2.0]
oprob = ODEProblem(two_state_model, u0, tspan, ps)

oprob = remake(oprob; u0 = [X1 => 5.0])
oprob[X1]
```
Symbolic variables can be used to access or update species or parameters for all the cases when `Symbol`s can (including when using `remake` or e.g. `getu`).

An advantage when quantities are represented as symbolic variables is that symbolic expressions can be formed and used to index a structure. E.g. here we check the combined initial concentration of $X$ ($X1 + X2$) in our two-state problem:
```@example structure_indexing_symbolic_variables
oprob[X1 + X2]
```

Just like symbolic variables can be used to directly interface with a structure, symbolic variables stored in `ReactionSystem` models can be used:
```@example structure_indexing_symbolic_variables
oprob[two_state_model.X1 + two_state_model.X2]
```
This can be used to form symbolic expressions using model quantities when a model has been created using the DSL (as an alternative to @unpack). Alternatively, [creating an observable](@ref dsl_advanced_options_observables), and then interface using its `Symbol` representation, is also possible.

!!! warn
    When accessing a simulation structure using symbolic variables from a `ReactionSystem` model, such as `rn.A` for `rn` a `ReactionSystem` and `A` a species within it, ensure that the model is complete.
