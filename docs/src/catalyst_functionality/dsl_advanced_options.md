# [The Catalyst DSL - Advanced Features and Options](@id dsl_advanced_options)
Within the Catalyst DSL, each line can represent either *a reaction* or *an option*. The [previous DSL tutorial](@ref dsl_description) described how to create reactions. This one will focus on options. These are typically used to supplant a model with additional information. Examples include the declaration of initial condition/parameter default values, or the creation observables or events. 

All options designations begins with declaration starting with `@`, followed by its input. E.g. the `@observables` options allows for the generation of observables. Each option can only be used once within each use of `@reaction_network`. A full list of options can be found [here](@ref ref), with most (but not all) being described in more detail below. This tutorial will also describe some additional advanced DSL features that does not directly include using an option. 

As a first step, we import Catalyst (which is required to run the tutorial):
```@example dsl_advanced_1
using Catalyst
```

## [Explicit specification of network species and parameters](@id dsl_advanced_options_declaring_species_and_parameters)
[Previously](@ref ref), we mentioned that the DSL automatically determines which symbols corresponds to species and which to parameters. This is done by designating everything that appear as either a substrate or a product as a species, and all remaining quantities as parameters (i.e. those only appearing within rates or [stoichiometric constants](@ref ref)). Sometimes, one might want to manually override this default behaviour for a given symbol. I.e. consider the following model, where the conversion of a protein `P` from its inactive form (`Pᵢ`) to its active form (`Pₐ`) is catalysed by an enzyme (`E`). Using the most natural description:
```@example dsl_advanced_1
catalysis_sys = @reaction_network begin
  k*E, Pᵢ --> Pₐ
end
```
`E` (as well as `k`) will be considered a parameter, something we can confirm directly:
```@example dsl_advanced_1
parameters(catalysis_sys)
```
If we want `E` to be considered a species, we can designate this using the `@species` option:
```@example dsl_advanced_1
catalysis_sys = @reaction_network begin
  @species E(t)
  k*E, Pᵢ --> Pₐ
end
parameters(catalysis_sys)
```
!!! note
    When declaring species using the `@species` option, the species symbol must be followed by `(t)`. The reason is that species are time-dependent variables, and this time-dependency must be explicitly specified ([designation of non-`t` dependant species is also possible](@ref ref)).

Similarly, the `@parameters` option can be used to explicitly designate something as a parameter:
```@example dsl_advanced_1
catalysis_sys = @reaction_network begin
  @parameters k
  k*E, Pᵢ --> Pₐ
end
```
Here, while `k` is explicitly defined as a parameter, no information is provided about `E`. Hence, the default case will be used (setting `E` to a parameter). The `@species` and `@parameter` options can be used simultaneously (although a quantity cannot be declared *both* as a species and a parameter). They may be followed by a full list of all species/parameters, or just a subset.

While designating something which would default to a parameter as a species is straightforward, the reverse (creating a parameter which occur as a substrate or product) is more involved. This is, however, possible, and described [here](@ref dsl_advanced_options_constant_species).

Rather than listing all species/parameters on a single line after the options, a `begin ... end` block can be used (listing one species/parameter on each line). E.g. in the following example we use this notation to explicitly designate all species and parameters of the system:
```@example dsl_advanced_1
catalysis_sys = @reaction_network begin
  @species begin 
    E(t)
    Pᵢ(t)
    Pₐ(t)
  end
  @parameters begin
    k
  end
  k*E, Pᵢ --> Pₐ
end
```

A side-effect of using the `@species` and `@parameter` options is that they specify *the order in which the species and parameters are stored*. I.e. lets check the order of the parameters in the parameters in the following dimerisation model:
```@example dsl_advanced_1
dimerisation = @reaction_network begin
  (p,d), 0 <--> X
  (kB,kD), 2X <--> X2
end
parameters(dimerisation)
```
The default order is typically equal to the order with which the parameters (or species) are encountered in the DSL (this is, however, not guaranteed). If we specify the parameters using `@parameters`, the order used within the option is used instead:
```@example dsl_advanced_1
dimerisation = @reaction_network begin
  @parameters kB kD p d
  (p,d), 0 <--> X
  (kB,kD), 2X <--> X2
end
parameters(dimerisation)
```
!!! danger
    Generally, Catalyst and the SciML ecosystem *does not* guarantee that parameter and species order is preserved throughout various operations on the model. Writing programs that depend on these orders is *strongly discouraged*. There are, however, some legacy packages which still depends on order (one example is provided [here](@ref ref)). In these situations, this might be useful. However, in these cases, it is recommended that the user is extra wary, and also checks the order manually. 

!!! note
    The syntax of the `@species` and `@parameters` options is identical to that used by the `@species` and `@parameters` macros [used in programmatic modelling in Catalyst](@ref programmatic_CRN_construction) (for e.g. designating metadata or initial conditions). Hence, if one have learnt how to specify species/parameters using either approach, that knowledge can be transferred to the other one.

Generally, there are four main reasons for specifying species/parameters using the `@species` and `@parameters` option:
1. To designate a quantity, that would otherwise have defaulted to a parameter, as a species.
2. To designate default values for parameters/species initial conditions (described [here](@ref dsl_advanced_options_default_vals)).
3. To designate metadata for species/parameters (described [here](@ref dsl_advanced_options_species_and_parameters_metadata)).
4. To designate a species or parameters that does not occur in reactions, but are still part of the model (e.g a [parametric initial condition](@ref dsl_advanced_options_parametric_initial_conditions))

!!!! warn
    Catalyst's DSL automatically infer species and parameters from the input. However, it only does so for *quantities that appear in reactions*. Until now this has not been relevant. However, this tutorial will demonstrate cases where species/parameters that are not part of reactions are used. These *must* be designated using either the `@species` or `@parameters` options (or the `@variables` option, which is described [later](@ref dsl_advanced_options_variables)).

### [Setting default values for species and parameters](@id dsl_advanced_options_default_vals)
When declaring species/parameters using the `@species` and `@parameters` options, one can also assign them default values (by appending them with `=` followed by the desired default value). E.g here we set `X`'s default initial condition value to $1.0$, and `p` and `d`'s default values to $1.0$ and $0.2$, respectively:
```@example dsl_advanced_defaults
using Catalyst # hide
rn = @reaction_network begin
  @species X(t)=1.0
  @parameters p=1.0 d=0.1
  (p,d), 0 <--> X
end
```
Next, if we simulate the model, we do not need to provide values for species or parameters that have default values. In this case all have default values, so both `u0` and `ps` can be empty vectors:
```@example dsl_advanced_defaults
using OrdinaryDiffEq, Plots
u0 = []
tspan = (0.0, 10.0)
p = []
oprob = ODEProblem(rn, u0, tspan, p)
sol = solve(oprob, Tsit5())
plot(sol)
```
It is still possible to provide values for some (or all) initial conditions/parameters in `u0`/`ps` (in which case these overrides the default values):
```@example dsl_advanced_defaults
u0 = [:X => 4.0]
p = [:d => 0.5]
oprob = ODEProblem(rn, u0, tspan, p)
sol = solve(oprob, Tsit5())
plot(sol)
```
It is also possible to declare a model with default values for only some initial conditions/parameters:
```@example dsl_advanced_defaults
using Catalyst # hide
rn = @reaction_network begin
  @species X(t)=1.0
  (p,d), 0 <--> X
end

tspan = (0.0, 10.0)
p = [:p => 1.0, :D => 0.2]
oprob = ODEProblem(rn, u0, tspan, p)
sol = solve(oprob, Tsit5())
plot(sol)
```
API for checking the default values of species and parameters can be found [here](@ref ref).

### [Setting parametric initial conditions](@id dsl_advanced_options_parametric_initial_conditions)
In the previous section, we designated default values for species' initial conditions and parameters. However, the right-hand side of the designation accepts any valid expression (not only numeric values). While this can be used to set up some advanced default values, the most common use-case is to designate a species's initial condition as a parameter. E.e. in the following example we represent the initial condition of `X` using the parameter `X₀`. 
```@example dsl_advanced_defaults
rn = @reaction_network begin
  @species X(t)=X₀
  @parameters X₀
  (p,d), 0 <--> X
end
```
Please note that as the parameter `X₀` does not occur as part of any reactions, Catalyst's DSL cannot infer whether it is a species or a parameter. This must hence be explicitly declared. We can now simulate our model while providing `X`'s value through the `X₀` parameter:
```@example dsl_advanced_defaults
u0 = []
p = [:X₀ => 1.0, :p => 1.0, :d => 0.5]
oprob = ODEProblem(rn, u0, tspan, p)
sol = solve(oprob, Tsit5())
plot(sol)
```
It is still possible to designate $X$'s value in `u0`, in which case this overrides the default value.
```@example dsl_advanced_defaults
u0 = [:X => 0.5]
p = [:X₀ => 1.0, :p => 1.0, :d => 0.5]
oprob = ODEProblem(rn, u0, tspan, p)
sol = solve(oprob, Tsit5())
plot(sol)
```
Please note that `X₀` is still a parameter of the system, and as such its value must still be designated to simulate the model (even if it is not actually used).

### [Designating metadata for species and parameters](@id dsl_advanced_options_species_and_parameters_metadata)
Catalyst permits the user to define *metadata* for species and parameters. This permits the user to assign additional information to these, which can be used for a variety of purposes. Some Catalyst features depend on using metadata (with each such case describing specifically how this is done). Here we will introduce how to set metadata, and describe some common metadata types. 

Whenever a species/parameter is declared using the `@species`/`@parameters` options, it can be followed by a `[]` within which the metadata is given. Each metadata entry consists of the metadata's name, followed by a `=`, followed by its value. E.g. the `description` metadata allows you to attach a [`String`](https://docs.julialang.org/en/v1/base/strings/) to a species/parameter. Here we create a simple model where we add descriptions to all species and parameters.
```@example dsl_advanced_metadata
using Catalyst # hide
two_state_system = @reaction_network begin
  @species Xi(t) [description="The X's inactive form"] Xa(t) [description="The X's active form"]
  @parameters kA [description="X's activation rate"] kD [description="X's deactivation rate"]
  (ka,kD), Xi <--> Xa
end
```
A metadata can be given to only a subset of a system's species/parameters, and a quantity can be given several metadata entries. To give several metadata, separate each by a `,`. Here we only provide a description for `kA`, for which we also provide a [bounds metadata](@ref https://docs.sciml.ai/ModelingToolkit/dev/basics/Variable_metadata/#Bounds),
```@example dsl_advanced_metadata
two_state_system = @reaction_network begin
  @parameters kA [description="X's activation rate", bound=(0.01,10.0)]
  (ka,kD), Xi <--> Xa
end
```

It is possible to add both default values and metadata to a parameter/species. In this case, first provide the default value, next the metadata. I.e. to in the above example set $kA$'s default value to $1.0$ we use
```@example dsl_advanced_metadata
two_state_system = @reaction_network begin
  @parameters kA=1.0 [description="X's activation rate", bound=(0.01,10.0)]
  (ka,kD), Xi <--> Xa
end
```

When designating metadata for species/parameters in `begin ... end` blocks the syntax changes slight. Here, a `,` must be inserted before the metadata (but after any potential default value). I.e. a version of the previous example can be written as
```@example dsl_advanced_metadata
two_state_system = @reaction_network begin
  @parameters begin
    kA, [description="X's activation rate", bound=(0.01,10.0)]
    kD = 1.0, [description="X's deactivation rate"]
  end
  (ka,kD), Xi <--> Xa
end
```

Each metadata has its own getter functions. E.g. we can get the description of the parameter `kA` using `getdescription` (here we use [system indexing](@ref ref) to access the parameter):
```@example dsl_advanced_metadata
getdescription(two_state_system.kA)
```

It is not possible for the user to directly designate their own metadata. These have to first be added to Catalyst. Doing so is somewhat involved, and described in detail [here](@ref ref). A full list of metadata that can be used for species and/or parameters can be found [here](@ref ref).

### [Designating constant-valued/fixed species parameters](@id dsl_advanced_options_constant_species)

Catalyst enables the designation of parameters as `constantspecies`. These parameters can be used as species in reactions, however, their values are not changed by the reaction and remain constant throughout the simulation (unless changed by e.g. the [occurrence of an event]@ref ref). Practically, this is done by setting the parameter's `isconstantspecies` metadata to `true`. Here, we create a simple reaction where the species `X` is converted to `Xᴾ` at rate `k`. By designating `X` as a constant species parameter, we ensure that its quantity is unchanged by the occurrence of the reaction.
```@example dsl_advanced_constant_species
using Catalyst # hide
rn = @reaction_network begin
  @parameters X [isconstantspecies=true]
  k, X --> Xᴾ
end
```
We can confirm that $Xᴾ$ is the only species of the system:
```@example dsl_advanced_constant_species
species(rn)
```
Here, the produced model is actually identical to if $X$ had simply been a parameter in the reaction's rate:
```@example dsl_advanced_constant_species
rn = @reaction_network begin
  k*X, 0 --> Xᴾ
end
```

A common use-case for constant species are when modelling systems where some species are present in such surplus that their amounts the reactions' effect on it is negligible. A system which is commonly modelled this way is the [Brusselator](https://en.wikipedia.org/wiki/Brusselator).

### [Specifying non-species variables](@id dsl_advanced_options_variables)
--- MOVE THIS TO HYBRID SECTION ---
Chemical reaction network (CRN) models (which Catalyst creates) described how *species* are affected by the occurrence of *reaction events*. When they are converted to ODEs, the species become unknowns of the system. However, Catalyst permits the creation of [hybrid CRN models](@ref ref). These describe phenomenons which can only partially be modelled using CRNs. An example may be a bacterium. Here, we can use a CRN to model some internal system (e.g. controlling its growth rate). However, we might also want to model the bacterium's volume. Here, the volume cannot be considered a species (as it cannot plausibly be a reaction reactant). Instead, we should model it as a normal variable. Here, Catalyst provides the `@variables` option for adding non-species variables to the system. E.g. to create a model where a single growth factor ($G$) is produced and degraded, and where we also have a single volume variables ($V$) we can use:
```@example dsl_advanced_variables
using Catalyst # hide
rn = @reaction_network begin
  @variables V(t)
  (p,d), 0 <--> G
end
```
Note that $V$ (like species) is time-dependant, and (like species) must be declared as such when the `@variables` option is used. We can now simulate our model (remembering to provide a value for $V$ as well as $G$):
```@example dsl_advanced_variables
using OrdinaryDiffEq, Plots
u0 = [:G => 0.1, :V => 1.0]
tspan = (0.0, 10.0)
p = [:p => 1.0, :d => 0.5]
oprob = ODEProblem(rn, u0, tspan, p)
sol = solve(oprob, Tsit5())
plot(sol)
```
Here, we have not actually described how $V$ interacting with our model, or how its value may change. Primarily variables are declared as part of hybrid CRN/equation modelling, which is described in more detail [here](@ref ref).  

You can set metadata and default initial condition values for variables using the same syntax as used for parameters and species.

You can use the `variables` and `species` functions to retrieve a model's variables and species, respectively. The `unknowns` function retrieves both.

## [Setting reaction metadata](@id dsl_advanced_options_reaction_metadata)
--- DISCUSS WHAT TO DO WITH THIS SECTION ---
Reactions can also have metadata. This is described in detail [here](@ref ref).

## [Naming reaction networks](@id dsl_advanced_options_naming)
Each reaction network model has a name. It can be accessed using the `nameof` function. By default, some generic name is used:
```@example dsl_advanced_names
using Catalyst # hide
rn = @reaction_network begin
  (p,d), 0 <--> X
end
nameof(rn)
```
A specific name can be given as an argument between the `@reaction_network` and the `begin`. E.g. to name a network `my_network` we can use:
```@example dsl_advanced_names
rn = @reaction_network my_network begin
  (p,d), 0 <--> X
end
nameof(rn)
```

A consequence of generic names being used by default is that networks, even if seemingly identical, by default are not. E.g.
```@example dsl_advanced_names
rn1 = @reaction_network begin
  (p,d), 0 <--> X
end
rn2 = @reaction_network begin
  (p,d), 0 <--> X
end
rn1 == rn2
```
The reason can be confirmed by checking that their respective (randomly generated) names are different:
```@example dsl_advanced_names
nameof(rn1) == nameof(rn2)
```
By designating the networks to have the same name, however, identity is achieved.
```@example dsl_advanced_names
rn1 = @reaction_network my_network begin
  (p,d), 0 <--> X
end
rn2 = @reaction_network my_network begin
  (p,d), 0 <--> X
end
rn1 == rn2
```

Setting model names is primarily useful for [hierarchical modelling](@ref ref), where network names are appended to the display name of subnetworks' species and parameters.

## [Creating observables](@id dsl_advanced_options_observables)
Sometimes one might want to use observable variables. These are variables with values that can be computed directly from a system's state (rather than having their values implicitly given by reactions or equations). Observables can be designated using the `@observables` option. Here, the `@observables` option is followed by a `begin ... end` block with one line for each observable. Each line first gives the observable, followed by a `~` (*not* a `=`!), followed by an expression describing how to compute it.

Let us consider a model where two species (`X` and `Y`) can bind to form a complex (`XY`, which also can dissociate back into `X` and `Y`). If we wish to create a representation for the total amount of `X` and `Y` in the system, we can do this by creating observables `Xtot` and `Ytot`:
```@example dsl_advanced_observables
using Catalyst # hide
rn = @reaction_network begin
  @observables begin
    Xtot ~ X + XY
    Ytot ~ Y + XY
  end
  (kB,kD), X + Y <--> XY
end
```
We can now simulate our model using normal syntax (initial condition values for observables should not, and can not, be provided):
```@example dsl_advanced_observables
using OrdinaryDiffEq
u0 = [:X => 1.0, :Y => 2.0, :XY => 0.0]
tspan = (0.0, 10.0)
ps = [:kB => 1.0, :kD => 1.5]
oprob = ODEProblem(rn, u0, tspan, ps)
sol = solve(oprob, Tsit5())
nothing # hide
```

Next, we can use [symbolic indexing](@ref simulation_structure_interfacing) of our solution object, but with the observable as input. E.g. we can use 
```@example dsl_advanced_observables
sol[:Xtot]
```
to get a vector with `Xtot`'s value throughout the simulation. We can also use
```@example dsl_advanced_observables
using Plots
plot(sol; idxs = [:Xtot, :Ytot])
```
to plot the observables (rather than the species).

Observables can be defined using complicated expression containing species, parameters, and [variables](@ref ref) (but not other observables). In the following example (which uses a [parametric stoichiometry](@ref ref)) `X` polymerises to form a complex `Xn` containing `n` copies of `X`. Here, we create an observable describing the total number of `X` molecules in the system:
```@example dsl_advanced_observables
rn = @reaction_network begin
  @observables Xtot ~ X + n*Xn
  (kB,kD), n*X <--> Xn
end
nothing # hide
```
!!!
    If only a single observable is declared, the `begin .. end` block is not required and the observable can be declared directly after the `@observables` option.

[Metadata](@ref dsl_advanced_options_species_and_parameters_metadata) can be supplied to an observable directly after the its declaration (but before its formula). If so, the metadata must be separated from the observable with a `,`, and the observable plus the metadata encapsulated by `()`. E.g. to add a [description metadata](@ref dsl_advanced_options_species_and_parameters_metadata) to our observable we can use
```@example dsl_advanced_observables
rn = @reaction_network begin
  @observables (Xtot, [description="The total amount of X in the system."]) ~ X + n*Xn
  (kB,kD), n*X <--> Xn
end
nothing # hide
```

Observables are by default considered [variables](@ref ref) (not species). To designate them as a species, they can be pre-declared using the `@species` option. I.e. Here `Xtot` becomes a species:
```@example dsl_advanced_observables
rn = @reaction_network begin
  @species Xtot(t)
  @observables (Xtot, [description="The total amount of X in the system."]) ~ X + n*XnXY  
  (kB,kD), n*X <--> Xn
end
nothing # hide
```

Some final notes regarding observables:
- The left-hand side of the observable declaration must contain a single symbol only (with the exception of metadata, which can also be supplied).
- All quantities appearing on the right-hand side must be declared elsewhere within the `@reaction_network` call (either by being part of a reaction, or through the `@species`, `@parameters`, or `@variables` options).
- Observables may not depend on other observables.
- Observables have their [dependent variable(s)](@ref ref) automatically assigned as the union of the dependent variables of the species and variables on which it depends.

## [Creating events](@id dsl_advanced_options_events)
--- MOVE TO SEPARATE DOC PAGE ---
Sometimes one wishes to model events, describing things that can happen to it during a simulation.
 - A chemical system where an amount of some species is added at a time point after the simulation's initiation.
 - A simulation of a circadian rhythm, where light is turned on/off every 12 hours.
 - A cell divides when some size variable reaches a certain threshold, halving the amount of each species in the system.

Events are divided into *continuous* and *discrete* events, and these can be added directly to a system using the `continuous_events` and `discrete_events` options. Events can also be modelled through *callbacks*. These are different in that they are supplied in the simulation step (rather than on system creation), and generally provide more flexibility in how they may affect the system. Callbacks are described on a separate [page](@ref advanced_simulations_callbacks).

The notation described below for creating continuous and discrete events is the same which is used in [ModelingToolkit to create events](https://docs.sciml.ai/ModelingToolkit/stable/basics/Events/), and which is used for [events for programmatic model creation](@ref ref).

### [Continuous vs discrete events](@id dsl_advanced_options_events_continuous_vs_discrete)
Both continuous and discrete events combine some condition (for triggering the event) with some affect (describing their effects on the system). They differ in the following ways:
- They use slightly different notation.
- Discrete events' conditions are checked at *the end of* each simulation time step. For continuous events, the simulation instead finds the *exact time point* when the event is triggered at.
- Continuous events cannot be supplied to jump simulations.

### [Continuous events](@id dsl_advanced_options_events_continuous)
Let us consider a simple system where species `X` degraded at a constant rate `d`. Next, we wish to add an event which adds `2.0` units of `X` whenever `X` reaches a critical threshold `1.0`. This can be done in the following manner:
```@example dsl_advanced_events
using Catalyst # hide
rn = @reaction_network begin
  @continuous_events begin
    X ~ 1.0 => [X ~ X + 2.0]
  end
  d, X --> 0
end
nothing # hide
```
Here, the `@continuous_events` option is followed by a `begin ... end` block. Next, each line corresponds to a separate event. Each event is created in the following manner:
- It combines a *condition* (denoting when the event will happen) with one (or several) *affects* (denoting what the event does to the system when it happens). 
- The condition (left) and affect (right) are separated by a `=>`.
- The condition takes the form of an [equation](). Here, the event is triggered whenever the equation's two sides (separated by a `~`) are equal.
- The affect(s) are enclosed within `[]`. If there are multiple affects, these are separated by `,` (the example above contains a single affect).
- Each affect is a single equation that describes how a parameter, species, or [variable](@ref dsl_advanced_options_variables) is updated when the event is triggered.
- Each affect's equation's left-hand side must contain only the parameter/species/variable whose value should be updated.
- Each affect's equation's right-hand side is an expression describing its updated value.

We can simulate the model we declared, just like any other model:
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
    X ~ X_thres => [X ~ X + X_add]
  end
  d, X --> 0
end
nothing # hide
```
Here, since `X_thres` and `X_add` do not appear in any reactions, Catalyst cannot determine whether they are parameters or species. hence, they must be [explicitly designated as parameters by using the `@parameters` option](@ref dsl_advanced_options_declaring_species_and_parameters). Next, these can be designated as any value, and supplied to the `ODEProblem`:
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
    X ~ Y => [X ~ X - 1.0, Y ~ Y + 1.0]
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

!!!warn
    A single event (continuous or discrete) can update the value of either (one or several) species (and variables), or of (one or several) parameters. It is not possible for an event to update the values of both species/variables and parameters.

In the above examples we have modelled a system with a single event. In these cases, the `begin end` block is not required, and the event can be provided on the same line as the `@continuous_events` option:
```@example dsl_advanced_events
rn = @reaction_network begin
  @continuous_events  X ~ Y => [X ~ X - 1.0, Y ~ Y + 1.0]
  p, 0 --> X
  d, Y --> 0
end
nothing # hide
```

### [Discrete events](@id dsl_advanced_options_events_discrete)
Just like [continuous events](dsl_advanced_options_events_continuous), discrete events combine a condition with one or more affect statements. Here, discrete events' affects are created identically to those for continuous events. Discrete events' conditions are different. There exist 3 different types of discrete events, each with a different type of condition. All three types are created using the `@discrete_events` option, and a single system can contain a mix of all types. The three types are:
- Preset-time discrete events.
- Periodic discrete events.
- Conditional discrete events.

#### [Preset-time discrete events](@id dsl_advanced_options_events_discrete_presettime)
*Present-time events* are events that happen at specific time points. Here, the condition is a vector with all the time points at which an event is triggered. E.g. here we create a production/degradation loop, where `2.0` units of `X` is added at time points `3.0` and `7.0`
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

The preset time points can also be parameters (in which case, they have to be [designated as such using the `@parameters` option](@ref dsl_advanced_options_declaring_species_and_parameters)):
```@example dsl_advanced_events
rn = @reaction_network begin
  @parameters t1 t2
  @discrete_events begin
    [t1, t2] => [X ~ X + 2.0]
  end
  (p,d), 0 <--> X
end
nothing
```

#### [Periodic discrete events](@id dsl_advanced_options_events_discrete_periodic)
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
Like for preset-time events, periodic events' affects may contain parameters.

#### [Conditional discrete events](@id dsl_advanced_options_events_discrete_conditional)
Finally, discrete events' condition may be a boolean expression (consisting of parameters, species, variables, and the time variable). Let's say that we want to create an event which, if the concentration of `X` is below a threshold `1.0`, adds `1.0` units of `X` to the system, then we can use the following discrete event:
```@example dsl_advanced_events
rn = @reaction_network begin
  @discrete_events begin
    X < 1.0 => [X ~ X + 2.0]
  end
  d, X --> 0
end
```
If we simulate the system using the same conditions as for our [similar, continuous, example](@ref dsl_advanced_options_events_continuous) we get the same result:
```@example dsl_advanced_events
u0 = [:X => 2.0]
tspan = (0.0, 10.0)
ps = [:d => 1.0]

oprob = ODEProblem(rn, u0, tspan, ps)
sol = solve(ODEProblem, Tsit5())
plot(sol)
```
So, how is modelling this event as a discrete or continuous event different? There are four differences:
1) For continuous events, the simulation method finds the exact time point when the condition triggers. Discrete events are triggered at the first time step when the condition holds.
2) This discrete event will be triggered whenever `X < 1.0` holds, not just when the concentration of `X` passes the threshold. E.g. it will be triggered if the initial concentration of `X` is less than `1.0`.
3) Only the discrete event event can be used with jump simulations.
4) The discrete event can be used to create more advanced conditions.

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
!!!note
    When we form composite boolean conditions for conditional discrete events, we use `&` to denote the AND operator (not `&&`, as this is currently not supported).

!!!warn
    Generally, discrete events including equality (`==`) will not be triggered. The reason is that the condition is only checked at the end of every time step. Hence, unless special precautions are taken to ensure that the simulator stops when the condition holds, it will unlikely be triggered.


## [Specifying non-time independent variables](@id dsl_advanced_options_ivs)

As [described elsewhere](@ref ref), Catalyst's `ReactionSystem` models depends on a *time independent variable*, and potentially one or more *spatial independent variables*. By default, the independent variable `t` is used. We can declare another independent variable (which is automatically used as the default one) using the `@ivs` option. E.g. to use `τ` instead of `t` we can use
```@example dsl_advanced_ivs
using Catalyst # hide
rn = @reaction_network begin
  @ivs τ
  (ka,kD), Xi <--> Xa
end
nothing # hide
```
We can confirm that `Xi` and `Xa` depend on `τ` (and not `t`):
```@example dsl_advanced_ivs
species(rn)
```

It is possible to designate several independent variables using `@ivs`. If so, the first one is considered the default (time) independent variable, while the following one(s) are considered spatial independent variable(s). If we want some species to depend on a non-default independent variable, this has to be explicitly declared:
```@example dsl_advanced_ivs
rn = @reaction_network begin
  @ivs τ x
  @species X(x) Y(x)
  (p1,d1), 0 <--> X
  (p2,d2), 0 <--> Y
end
species(rn)
```
It is also possible to have species which depends on several independent variables:
```@example dsl_advanced_ivs
rn = @reaction_network begin
  @ivs t x
  @species Xi(t,x) Xa(t,x)
  (ka,kD), Xi <--> Xa
end
species(rn)
```

!!! note
    Setting spatial independent variables is primarily intended for modelling of spatial systems on continuous domains.  Catalyst's support for this is currently under development. Hence, the utility of specifying spatial independent variables is currently limited.