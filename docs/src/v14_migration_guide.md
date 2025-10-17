# Version 14 Migration Guide

Catalyst is built on the [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) modelling language. A recent update of ModelingToolkit from version 8 to version 9 has required a corresponding update to Catalyst (from version 13 to 14). This update has introduced a couple of breaking changes, all of which will be detailed below.

!!! note
    Catalyst version 14 also introduces several new features. These will not be discussed here, however, they are described in Catalyst's [history file](https://github.com/SciML/Catalyst.jl/blob/master/HISTORY.md).

## System completeness

In ModelingToolkit v9 (and thus also Catalyst v14) all systems (e.g. `ReactionSystem`s and `ODESystem`s) are either *complete* or *incomplete*. Complete and incomplete systems differ in that

- Only complete systems can be used as inputs to simulations or most tools for model analysis.
- Only incomplete systems can be [composed with other systems to form hierarchical models](@ref compositional_modeling).

A model's completeness depends on how it was created:

- Models created programmatically (using the `ReactionSystem` constructor) are *not marked as complete* by default.
- Models created using the `@reaction_network` DSL are *automatically marked as complete*.
- To *use the DSL to create models that are not marked as complete*, use the `@network_component` macro (which in all other aspects is identical to `@reaction_network`).
- Models generated through the `compose` and `extend` functions are *not marked as complete*.

Furthermore, any systems generated through e.g. `convert(ODESystem, rs)`  are *not marked as complete*.

Complete models can be generated from incomplete models through the `complete` function. Here is a workflow where we take completeness into account in the simulation of a simple birth-death process.

```@example v14_migration_1
using Catalyst
t = default_t()
@species X(t)
@parameters p d
rxs = [
    Reaction(p, [], [X]),
    Reaction(d, [X], [])
]
@named rs = ReactionSystem(rxs, t)
```

Here we have created a model that is not marked as complete. If our model is ready (i.e. we do not wish to compose it with additional models) we mark it as complete:

```@example v14_migration_1
rs = complete(rs)
```

Here, `complete` does not change the input model, but simply creates a new model that is tagged as complete. We hence overwrite our model variable (`rs`) with `complete`'s output. We can confirm that our model is complete using the `Catalyst.iscomplete` function:

```@example v14_migration_1
Catalyst.iscomplete(rs)
```

We can now go on and use our model for e.g. simulations:

```@example v14_migration_1
using OrdinaryDiffEqDefault, Plots
u0 = [X => 0.1]
tspan = (0.0, 10.0)
ps = [p => 1.0, d => 0.2]
oprob = ODEProblem(rs, u0, tspan, ps)
sol = solve(oprob)
plot(sol)
```

If we wish to first manually convert our `ReactionSystem` to an `ODESystem`, the generated `ODESystem` will *not* be marked as complete

```@example v14_migration_1
osys = convert(ODESystem, rs)
Catalyst.iscomplete(osys)
```

(note that `rs` must be complete before it can be converted to an `ODESystem` or any other system type)

If we now wish to create an `ODEProblem` from our `ODESystem`, we must first mark it as complete (using similar syntax as for our `ReactionSystem`):

```@example v14_migration_1
osys = complete(osys)
oprob = ODEProblem(osys, u0, tspan, ps)
sol = solve(oprob)
plot(sol)
```

Note, if we had instead used the [`@reaction_network`](@ref) DSL macro to build
our model, i.e.

```@example v14_migration_1
rs2 = @reaction_network rs begin
    p, ∅ --> X
    d, X --> ∅
end
```

then the model is automatically marked as complete

```@example v14_migration_1
Catalyst.iscomplete(rs2)
```

In contrast, if we used the [`@network_component`](@ref) DSL macro to build our
model it is not marked as complete, and is equivalent to our original definition of `rs`

```@example v14_migration_1
rs3 = @network_component rs begin
    p, ∅ --> X
    d, X --> ∅
end
Catalyst.iscomplete(rs3)
```

## Unknowns instead of states

Previously, "states" was used as a term for system variables (both species and non-species variables). MTKv9 has switched to using the term "unknowns" instead. This means that there have been a number of changes to function names (e.g. `states` => `unknowns` and `get_states` => `get_unknowns`).

E.g. here we declare a `ReactionSystem` model containing both species and non-species unknowns:

```@example v14_migration_2
using Catalyst
t = default_t()
D = default_time_deriv()
@species X(t)
@variables V(t)
@parameters p d Vmax

eqs = [
    Reaction(p, [], [X]),
    Reaction(d, [X], []),
    D(V) ~ Vmax - V*X*d/p
]
@named rs = ReactionSystem(eqs, t)
```

We can now use `unknowns` to retrieve all unknowns

```@example v14_migration_2
unknowns(rs)
```

Meanwhile, `species` and `nonspecies` (like previously) returns all species or non-species unknowns, respectively:

```@example v14_migration_2
species(rs)
```

```@example v14_migration_2
nonspecies(rs)
```

## Lost support for most units

As part of its v9 update, ModelingToolkit changed how units were handled. This includes using the package [DynamicQuantities.jl](https://github.com/SymbolicML/DynamicQuantities.jl) to manage units (instead of [Unitful.jl](https://github.com/PainterQubits/Unitful.jl), like previously).

While this should lead to long-term improvements, unfortunately, as part of the process support for most units was removed. Currently, only the main SI units are supported (`s`, `m`, `kg`, `A`, `K`, `mol`, and `cd`). Composite units (e.g. `N = kg/(m^2)`) are no longer supported. Furthermore, prefix units (e.g. `mm = m/1000`) are not supported either. This means that most units relevant to Catalyst (such as `µM`) cannot be used directly. While composite units can still be written out in full and used (e.g. `kg/(m^2)`) this is hardly user-friendly.

The maintainers of ModelingToolkit have been notified of this issue. We are unsure when this will be fixed, however, we do not think it will be a permanent change.

## Removed support for system-mutating functions

According to the ModelingToolkit system API, systems should not be mutable. In accordance with this, the following functions have been deprecated and removed: `addparam!`, `addreaction!`, `addspecies!`, `@add_reactions`, and `merge!`. Please use `ModelingToolkit.extend` and `ModelingToolkit.compose` to generate new merged and/or composed `ReactionSystems` from multiple component systems.

It is still possible to add default values to a created `ReactionSystem`, i.e. the `setdefaults!` function is still supported.

## New interface for creating time variable (`t`) and its differential (`D`)

Previously, the time-independent variable (typically called `t`) was declared using

```@example v14_migration_3
using Catalyst
@variables t
nothing # hide
```

MTKv9 has introduced a standard global time variable, and as such a new interface has been developed:

```@example v14_migration_3
t = default_t()
nothing # hide
```

Note that internally MTK9 actually represents `t` as a parameter, so the old approach should never be used. As the type of `t` is now considered internal to MTK, it should always be declared using the `default_t()` function.

Similarly, the time differential (primarily relevant when creating combined reaction-ODE models) used to be declared through

```@example v14_migration_3
D = Differential(t)
nothing # hide
```

where now one must use

```@example v14_migration_3
D = default_time_deriv()
nothing # hide
```

!!! note
    If you look at ModelingToolkit documentation, these defaults are instead retrieved using `using ModelingToolkit: t_nounits as t, D_nounits as D`. This will also work, however, in Catalyst we have opted to instead use the functions `default_t()` and `default_time_deriv()` as our main approach.

## New interface for accessing problem/integrator/solution parameter (and species) values

Previously, it was possible to directly index problems to query them for their parameter values. e.g.

```@example v14_migration_4
using Catalyst
rn = @reaction_network begin
    (p,d), 0 <--> X
end
u0 = [:X => 1.0]
ps = [:p => 1.0, :d => 0.2]
oprob = ODEProblem(rn, u0, (0.0, 1.0), ps)
nothing # hide
```

```julia
oprob[:p]
```

This is *no longer supported*. When you wish to query a problem (or integrator or solution) for a parameter value (or to update a parameter value), you must append `.ps` to the problem variable name:

```@example v14_migration_4
oprob.ps[:p]
```

Furthermore, a few new functions (`getp`, `getu`, `setp`, `setu`) have been introduced from [SymbolicIndexingInterface](https://github.com/SciML/SymbolicIndexingInterface.jl) to support efficient and systematic querying and/or updating of symbolic unknown/parameter values. Using these can *significantly* improve performance when querying or updating a value multiple times, for example within a callback. These are described in more detail [here](@ref simulation_structure_interfacing_functions).

For more details on how to query various structures for parameter and species values, please read [this documentation page](@ref simulation_structure_interfacing).

## Other changes

#### Modification of problems with conservation laws broken

While it is possible to update e.g. `ODEProblem`s using the [`remake`](@ref simulation_structure_interfacing_problems_remake) function, this is currently not possible if the `remove_conserved = true` option was used. E.g. while

```@example v14_migration_5
using Catalyst, OrdinaryDiffEqDefault
rn = @reaction_network begin
    (k1,k2), X1 <--> X2
end
u0 = [:X1 => 1.0, :X2 => 2.0]
ps = [:k1 => 0.5, :k2 => 3.0]
oprob = ODEProblem(rn, u0, (0.0, 10.0), ps; remove_conserved = true)
solve(oprob)
# hide
```

is perfectly fine, attempting to then modify any initial conditions or the value of the conservation constant in `oprob` will likely silently fail:

```@example v14_migration_5
oprob_remade = remake(oprob; u0 = [:X1 => 5.0]) # NEVER do this.
solve(oprob_remade)
# hide
```

This might generate a silent error, where the remade problem is different from the intended one (the value of the conserved constant will not be updated correctly).

This bug was likely present on earlier versions as well, but was only recently discovered. While we hope it will be fixed soon, the issue is in ModelingToolkit, and will not be fixed until its maintainers find the time to do so.

#### Depending on parameter order is even more dangerous than before

In early versions of Catalyst, parameters and species were provided as vectors (e.g. `[1.0, 2.0]`) rather than maps (e.g. `[p => 1.0, d => 2.0]`). While we previously *strongly* recommended users to use the map form (or they might produce unintended results), the vector form was still supported (technically). Due to recent internal ModelingToolkit updates, the purely numeric form is no longer supported and should never be used -- it will potentially lead to incorrect values for parameters and/or initial conditions. Note that if `rn` is a complete `ReactionSystem` you can now specify such mappings via `[rn.p => 1.0, rn.d => 2.0]`.

*Users should never use vector-forms to represent parameter and species values.*  <!-- markdownlint-disable-line MD036 -->

#### Additional deprecated functions

The `reactionparams`, `numreactionparams`, and `reactionparamsmap` functions have been deprecated.
