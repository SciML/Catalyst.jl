# [The Catalyst DSL - Advanced Features and Options](@id dsl_advanced_options)
Within the Catalyst DSL, each line can represent either *a reaction* or *an option*. The [previous DSL tutorial](@ref dsl_description) described how to create reactions. This one will focus on options. These are typically used to supply a model with additional information. Examples include the declaration of initial condition/parameter default values, or the creation of observables. 

All option designations begin with a declaration starting with `@`, followed by its input. E.g. the `@observables` option allows for the generation of observables. Each option can only be used once within each use of `@reaction_network`. A full list of options can be found [here](@ref ref), with most (but not all) being described in more detail below. This tutorial will also describe some additional advanced DSL features that do not involve using an option. 

As a first step, we import Catalyst (which is required to run the tutorial):
```@example dsl_advanced_explicit_definitions
using Catalyst
```

## [Explicit specification of network species and parameters](@id dsl_advanced_options_declaring_species_and_parameters)
[Previously](@ref ref), we mentioned that the DSL automatically determines which symbols correspond to species and which to parameters. This is done by designating everything that appears as either a substrate or a product as a species, and all remaining quantities as parameters (i.e. those only appearing within rates or [stoichiometric constants](@ref ref)). Sometimes, one might want to manually override this default behaviour for a given symbol. I.e. consider the following model, where the conversion of a protein `P` from its inactive form (`Pᵢ`) to its active form (`Pₐ`) is catalysed by an enzyme (`E`). Using the most natural description:
```@example dsl_advanced_explicit_definitions
catalysis_sys = @reaction_network begin
  k*E, Pᵢ --> Pₐ
end
```
`E` (as well as `k`) will be considered a parameter, something we can confirm directly:
```@example dsl_advanced_explicit_definitions
parameters(catalysis_sys)
```
If we want `E` to be considered a species, we can designate this using the `@species` option:
```@example dsl_advanced_explicit_definitions
catalysis_sys = @reaction_network begin
  @species E(t)
  k*E, Pᵢ --> Pₐ
end
parameters(catalysis_sys)
```
!!! note
    When declaring species using the `@species` option, the species symbol must be followed by `(t)`. The reason is that species are time-dependent variables, and this time-dependency must be explicitly specified ([designation of non-`t` dependant species is also possible](@ref ref)).

Similarly, the `@parameters` option can be used to explicitly designate something as a parameter:
```@example dsl_advanced_explicit_definitions
catalysis_sys = @reaction_network begin
  @parameters k
  k*E, Pᵢ --> Pₐ
end
```
Here, while `k` is explicitly defined as a parameter, no information is provided about `E`. Hence, the default case will be used (setting `E` to a parameter). The `@species` and `@parameter` options can be used simultaneously (although a quantity cannot be declared *both* as a species and a parameter). They may be followed by a full list of all species/parameters, or just a subset.

While designating something which would default to a parameter as a species is straightforward, the reverse (creating a parameter which occurs as a substrate or product) is more involved. This is, however, possible, and described [here](@ref dsl_advanced_options_constant_species).

Rather than listing all species/parameters on a single line after the options, a `begin ... end` block can be used (listing one species/parameter on each line). E.g. in the following example we use this notation to explicitly designate all species and parameters of the system:
```@example dsl_advanced_explicit_definitions
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
```@example dsl_advanced_explicit_definitions
dimerisation = @reaction_network begin
  (p,d), 0 <--> X
  (kB,kD), 2X <--> X2
end
parameters(dimerisation)
```
The default order is typically equal to the order with which the parameters (or species) are encountered in the DSL (this is, however, not guaranteed). If we specify the parameters using `@parameters`, the order used within the option is used instead:
```@example dsl_advanced_explicit_definitions
dimerisation = @reaction_network begin
  @parameters kB kD p d
  (p,d), 0 <--> X
  (kB,kD), 2X <--> X2
end
parameters(dimerisation)
```
!!! danger
    Generally, Catalyst and the SciML ecosystem *do not* guarantee that parameter and species order are preserved throughout various operations on a model. Writing programs that depend on these orders is *strongly discouraged*. There are, however, some legacy packages which still depend on order (one example is provided [here](@ref ref)). In these situations, this might be useful. However, in these cases, it is recommended that the user is extra wary, and also checks the order manually. 

!!! note
    The syntax of the `@species` and `@parameters` options is identical to that used by the `@species` and `@parameters` macros [used in programmatic modelling in Catalyst](@ref programmatic_CRN_construction) (for e.g. designating metadata or initial conditions). Hence, if one has learnt how to specify species/parameters using either approach, that knowledge can be transferred to the other one.

Generally, there are four main reasons for specifying species/parameters using the `@species` and `@parameters` options:
1. To designate a quantity, that would otherwise have defaulted to a parameter, as a species.
2. To designate default values for parameters/species initial conditions (described [here](@ref dsl_advanced_options_default_vals)).
3. To designate metadata for species/parameters (described [here](@ref dsl_advanced_options_species_and_parameters_metadata)).
4. To designate a species or parameters that do not occur in reactions, but are still part of the model (e.g a [parametric initial condition](@ref dsl_advanced_options_parametric_initial_conditions))

!!!! warn
    Catalyst's DSL automatically infer species and parameters from the input. However, it only does so for *quantities that appear in reactions*. Until now this has not been relevant. However, this tutorial will demonstrate cases where species/parameters that are not part of reactions are used. These *must* be designated using either the `@species` or `@parameters` options (or the `@variables` option, which is described [later](@ref ref)).

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
sol = solve(oprob)
plot(sol)
```
It is still possible to provide values for some (or all) initial conditions/parameters in `u0`/`ps` (in which case these overrides the default values):
```@example dsl_advanced_defaults
u0 = [:X => 4.0]
p = [:d => 0.5]
oprob = ODEProblem(rn, u0, tspan, p)
sol = solve(oprob)
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
sol = solve(oprob)
plot(sol)
```
API for checking the default values of species and parameters can be found [here](@ref ref).

### [Setting parametric initial conditions](@id dsl_advanced_options_parametric_initial_conditions)
In the previous section, we designated default values for initial conditions and parameters. However, the right-hand side of the designation accepts any valid expression (not only numeric values). While this can be used to set up some advanced default values, the most common use case is to designate a species's initial condition as a parameter. E.g. in the following example we represent the initial condition of `X` using the parameter `X₀`. 
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

It is possible to add both default values and metadata to a parameter/species. In this case, first provide the default value, and next the metadata. I.e. to in the above example set $kA$'s default value to $1.0$ we use
```@example dsl_advanced_metadata
two_state_system = @reaction_network begin
  @parameters kA=1.0 [description="X's activation rate", bound=(0.01,10.0)]
  (ka,kD), Xi <--> Xa
end
```

When designating metadata for species/parameters in `begin ... end` blocks the syntax changes slightly. Here, a `,` must be inserted before the metadata (but after any potential default value). I.e. a version of the previous example can be written as
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

A common use-case for constant species is when modelling systems where some species are present in such surplus that their amounts the reactions' effect on it is negligible. A system which is commonly modelled this way is the [Brusselator](https://en.wikipedia.org/wiki/Brusselator).

### [Designating parameter types](@id dsl_advanced_options_parameter_types)
Sometimes it is desired to designate that a parameter should have a specific [type](@ref ref). When supplying this parameter's value to e.g. an `ODEProblem`, that parameter will then be restricted to that specific type. Designating a type is done by appending the parameter with `::` followed by its type. E.g. in the following example we specify that the parameter `n` (the number of `X` molecules in the `Xn` polymer) must be an integer (`Int64`)
```@example dsl_advanced_parameter_types
using Catalyst # hide
polymerisation_network = @reaction_network begin
    @parameters n::Int64
    (kB,kD), n*X <--> Xn
end
nothing # hide
```
Generally, when simulating models with mixed parameter types, it is recommended to [declare parameter values as tuples, rather than vectors](@ref ref), e.g.:
```@example dsl_advanced_parameter_types
ps = (:kB => 0.2, :kD => 1.0, :n => 2)
nothing # hide
```

If a parameter has a type, metadata, and a default value, they are designated in the following order:
```@example dsl_advanced_parameter_types
polymerisation_network = @reaction_network begin
    @parameters n::Int64 = 2 [description="Parameter n, which is an integer and defaults to the value 2."]
    (kB,kD), n*X <--> Xn
end
nothing # hide
```

### [Vector-valued species or parameters](@id dsl_advanced_options_vector_variables)
Sometimes, one wishes to declare a large number of similar parameters or species. This can be done by *creating them as vectors*. E.g. below we create a [two-state system](@ref ref). However, instead of declaring `X1` and `X2` (and `k1` and `k2`) as separate entities, we declare them as vectors:
```@example dsl_advanced_vector_variables
using Catalyst # hide
two_state_model = @reaction_network begin
    @parameters k[1:2]
    @species X(t)[1:2]
    (k[1],k[2]), X[1] <--> X[2]
end
```
Now, we can also declare our initial conditions and parameter values as vectors as well:
```@example dsl_advanced_vector_variables
using OrdinaryDiffEq, Plots # hide
u0 = [:X => [0.0, 2.0]]
tspan = (0.0, 1.0)
ps = [:k => [1.0, 2.0]]
oprob = ODEProblem(two_state_model, u0, tspan, ps)
sol = solve(oprob)
plot(sol)
```

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

Setting model names is primarily useful for [hierarchical modelling](@ref ref), where network names are appended to the display names of subnetworks' species and parameters.

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

Observables can be defined using complicated expressions containing species, parameters, and [variables](@ref ref) (but not other observables). In the following example (which uses a [parametric stoichiometry](@ref ref)) `X` polymerises to form a complex `Xn` containing `n` copies of `X`. Here, we create an observable describing the total number of `X` molecules in the system:
```@example dsl_advanced_observables
rn = @reaction_network begin
  @observables Xtot ~ X + n*Xn
  (kB,kD), n*X <--> Xn
end
nothing # hide
```
!!!
    If only a single observable is declared, the `begin .. end` block is not required and the observable can be declared directly after the `@observables` option.

[Metadata](@ref dsl_advanced_options_species_and_parameters_metadata) can be supplied to an observable directly after its declaration (but before its formula). If so, the metadata must be separated from the observable with a `,`, and the observable plus the metadata encapsulated by `()`. E.g. to add a [description metadata](@ref dsl_advanced_options_species_and_parameters_metadata) to our observable we can use
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
  @observables Xtot ~ X + n*XnXY  
  (kB,kD), n*X <--> Xn
end
nothing # hide
```

Some final notes regarding observables:
- The left-hand side of the observable declaration must contain a single symbol only (with the exception of metadata, which can also be supplied).
- All quantities appearing on the right-hand side must be declared elsewhere within the `@reaction_network` call (either by being part of a reaction, or through the `@species`, `@parameters`, or `@variables` options).
- Observables may not depend on other observables.
- Observables have their [dependent variable(s)](@ref ref) automatically assigned as the union of the dependent variables of the species and variables on which it depends.

## [Specifying non-time independent variables](@id dsl_advanced_options_ivs)

As [described elsewhere](@ref ref), Catalyst's `ReactionSystem` models depend on a *time independent variable*, and potentially one or more *spatial independent variables*. By default, the independent variable `t` is used. We can declare another independent variable (which is automatically used as the default one) using the `@ivs` option. E.g. to use `τ` instead of `t` we can use
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
  @species X(τ) Y(x)
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
    Setting spatial independent variables is primarily intended for the modelling of spatial systems on continuous domains. Catalyst's support for this is currently under development. Hence, the utility of specifying spatial independent variables is limited.


## [Setting reaction metadata](@id dsl_advanced_options_reaction_metadata)
It is possible to supply reactions with *metadata*, containing some additional information of the reaction. A reaction's metadata follows after its declaration (first using the metadata's name, then a `=`, then its value) and is encapsulated by `[]` (where individual entries are separated by `,`). Here, we add a `description` metadata to the reactions of a birth-death process:
```@example dsl_advanced_reaction_metadata
using Catalyst # hide
bd_model = @reaction_network begin
  p, 0 --> X, [description="A production reaction"]
  d, X --> 0, [description="A degradation reaction"]
end
nothing # hide
```

When [bundling reactions](@ref dsl_description_reaction_bundling), reaction metadata can be bundled using the same rules as rates. Bellow we re-declare our birth-death process, but on a single line:
```@example dsl_advanced_reaction_metadata
bd_model = @reaction_network begin
  (p,d), 0 --> X, ([description="A production reaction"], [description="A degradation reaction"])
end
nothing # hide
```

Here we declare a model where we also provide a `misc` metadata (which can hold any quantity we require) to our birth reaction:
```@example dsl_advanced_reaction_metadata
bd_model = @reaction_network begin
  p, 0 --> X, [description="A production reaction", misc=:value]
  d, X --> 0, [description="A degradation reaction"]
end
nothing # hide
```

A reaction's metadata can be accessed using specific functions, e.g. `Catalyst.hasdescription` and `Catalyst.getdescription` can be used to check if a reaction have a description metadata, and to retrieve it, respectively:
```@example dsl_advanced_reaction_metadata
rx = @reaction p, 0 --> X, [description="A production reaction"]
Catalyst.getdescription(rx)
```

A list of all available reaction metadata can be found [here](@ref ref).