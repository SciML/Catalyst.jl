# [The Catalyst DSL - Advanced Features and Options](@id dsl_advanced_options)

Within the Catalyst DSL, each line can represent either *a reaction* or *an option*. The [previous tutorial](@ref dsl_description) described how to create reactions. This one will instead describe how to use options. These are typically used to supplant the model with additional information. Examples include the declaration of initial condition/parameter default values, or the creation observables or events. 

All options are designated begins with a symbol starting with `@`, followed by its input. E.g. the `@observables` options allows for the generation of observables. Each option can only be used once within each declaration of `@reaction_network`. A full list of options can be found [here], with most (but not all) being described in more detail below.

This tutorial will also describe some additional advanced DSL features that does not include using an option. As a first step, we import Catalyst (which is required to run the tutorial):
```@example dsl_advanced_1
using Catalyst
```

## [Explicit specification of network species and parameters](@id dsl_advanced_options_declaring_species_and_parameters)
[Previously](@ref ref), we mentioned that the DSL automatically determines which symbols corresponds to species and which to parameters. This is done by designating everything that appear as either a substrate or a product as a species, and all remaining symbols as parameters (i.e. those only appearing within rates or [stoichiometric constants](@ref ref)). Sometimes, one might want to manually override this default
behaviour for a given symbol. I.e. consider the following model, where the conversion of a protein $P$ from its inactive ($Pᵢ$) to its active ($Pₐ$) for for is catalysed by an enzyme $E$. Using the most natural description:
```@example dsl_advanced_1
catalysis_sys = @reaction_network begin
  k*E, Pᵢ --> Pₐ
end
```
`X` (as well as `k`) will be considered a parameter (this can be checked using e.g. `parameters(catalysis_sys)`). If we want $E$ to be considered a species, we can designate this using the `@species` option:
```@example dsl_advanced_1
catalysis_sys = @reaction_network begin
  @species E(t)
  k*E, Pᵢ --> Pₐ
end
```
!!! note
    When declaring species using the `@species` option, the species symbol must be followed by `(t)`. The reason is that species are time-dependent variables, and this time-dependency must be explicitly specified ([designation of non-time dependant species is also possible](@ref ref)).

Similarly, the `@parameters` option can be used to  
```@example dsl_advanced_1
catalysis_sys = @reaction_network begin
  @parameters k
  k*E, Pᵢ --> Pₐ
end
```
Here, while `k` is explicitly defined as a parameter, not information is provided about $E$. Hence, the default case will be used (setting $E$ to a parameter). The `@species` and `@parameter` options can be used simultaneously (although a symbol cannot be declared *both* as a species and a parameter). They may be followed by an extensive list of all species/parameters, or just a subset.

While designating something which would default to a parameter as a species is straightforward, the other direction (creating a parameter which occur as a substrate or product) is more involved. This is, however, possible, and described [here](@ref dsl_advanced_options_constant_species).

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
    Generally, Catalyst and the SciML ecosystem *does not* guarantee that parameter and species order is preserved throughout various operations on the model. Writing programs that depend on these orders is *strongly discouraged*. There are, however, some legacy packages which still depends on order (one example is provided [here](@ref ref)). In these situations, this might be useful. However, in these cases, it is recommended that the user is extra wary.

The syntax of the `@species` and `@parameters` options is identical to that used by the `@species` and `@parameters` macros [used in programmatic modelling in Catalyst](@ref programmatic_CRN_construction) (for e.g. designating metadata or initial conditions). Hence, if one have learnt how to specify species/parameters using either system, that knowledge can be transferred to the other one.

Generally, there are dour main reasons for specifying species/parameters using the `@species` and `@parameters` option:
1. To designate a symbol, that would otherwise have defaulted to a parameter, as a species.
2. To designate default values for parameters/species initial conditions (described [here](@ref dsl_advanced_options_default_vals)).
3. To designate metadata for species/parameters (described [here](@ref dsl_advanced_options_species_and_parameters_metadata)).
1. To designate a species or parameters that does not occur in reactions, but are still part of the model (e.g a [parametric initial condition](@ref dsl_advanced_options_parametric_initial_conditions))

!!!! warn
    Catalyst's DSL automatically infer species and parameters from the input. However, it only does so for *symbols that appear in reactions*. Until now this has not been relevant. However, this tutorial will demosntrate use cases where species/parameters, that are not part of reactions, are used. These *must* be designated using either the `@species` or `@parameters` options (or the `@variables` option, which is described [later](@ref dsl_advanced_options_variables)).

### [Setting default values for species and parameters](@id dsl_advanced_options_default_vals)
When declaring species/parameters using the `@species` and `@parameters` options, one can also assign them default values (by appending them with `=` and the desired default value). E.g here we set $X$'s default initial condition value to $1.0$, and $p$ and $d$'s default values to $1.0$ and $0.2$, respectively:
```@example dsl_advanced_defaults
using Catalyst # hide
rn = @reaction_network begin
  @species X(t)=1.0
  @parameters p=1.0 d=0.1
  (p,d), 0 <--> X
end
```
Next, if we simulate the model, we do not need to provide values for species/parameters which have default values. In this case all have default values, so both `u0` and `ps` are set to empty vectors:
```@example dsl_advanced_defaults
using OrdinaryDiffEq, Plots
u0 = []
tspan = (0.0, 10.0)
p = []
oprob = ODEProblem(rn, u0, tspan, p)
sol = solve(oprob, Tsit5())
plot(sol)
```
It is still possible to provide values for some (or all) initial conditions/parameters in `u0` and `ps` (in which case these overrides the default values):
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
API for checking the default values of a species or parameters can be found [here](@ref ref).

### [Setting parametric initial conditions](@id dsl_advanced_options_parametric_initial_conditions)
In the previous section, we designated default values for species' initial conditions and parameters. However, the right-hand side of the designation accepts any valid expression (not only numeric values). While this can be used to set up some advanced default values, the most common use-case is to designate a species's initial condition as a parameter. E.e. in the following example we represent the initial condition of $X$ using the parameter $X₀$. 
```@example dsl_advanced_defaults
rn = @reaction_network begin
  @species X(t)=X₀
  @parameters X₀
  (p,d), 0 <--> X
end
```
Please note that as the parameter $X₀$ does not occur as part of any reactions, Catalyst's DSL cannot infer whether it is a species or a parameter. This must hence be explicitly declared. We can now simulate our model while providing $X$'s value through the $X₀$ parameter:
```@example dsl_advanced_defaults
u0 = []
p = [:X₀ => 1.0, :p => 1.0, :d => 0.5]
oprob = ODEProblem(rn, u0, tspan, p)
sol = solve(oprob, Tsit5())
plot(sol)
```
It is still possible to designate $X$'s value in `u0`, in which case this overrides the default value. Please note that $X₀$ is still a parameter of the system, and its value must still be designated to simulate the model.
```@example dsl_advanced_defaults
u0 = [:X => 0.5]
p = [:X₀ => 1.0, :p => 1.0, :d => 0.5]
oprob = ODEProblem(rn, u0, tspan, p)
sol = solve(oprob, Tsit5())
plot(sol)
```

### [Designating metadata for species and parameters](@id dsl_advanced_options_species_and_parameters_metadata)
Catalyst permits the user to define *metadata* for species and parameters. This permits the user to assign additional information to these, which can be used for a variety of purposes. Some Catalyst features depend on using metadata (with each such case describing specifically ho this is done). Here we will introduce how to set metadata, and describe some common metadata types. 

Whenever a species/parameter is declared using the `@species`/`@parameters` options, it can be followed by a `[]` within which the metadata is given. Each metadata entry consists of teh metadata's name, followed by a `=`, followed by its value. E.g. the `description` metadata allows you to attach a `String` to a species/parameter. Here we create a simple model where we add descriptions to all species and parameters.
```@example dsl_advanced_metadata
using Catalyst # hide
two_state_system = @reaction_network begin
  @species Xi(t) [description="The species X's inactive form"] Xa(t) [description="The species X's active form"]
  @parameters kA [description="X's activation rate"] kD [description="X's deactivation rate"]
  (ka,kD), Xi <--> Xa
end
```
A metadata can be given to only a subset of a system's species/parameters. To give several metadata, separate each by a `,`. Here we remove some of the descriptions, and also add a [bounds metadata](@ref ref) to $kA$,
```@example dsl_advanced_metadata
two_state_system = @reaction_network begin
  @parameters kA [description="X's activation rate", bound=(0.01,10.0)] kD [description="X's deactivation rate"]
  (ka,kD), Xi <--> Xa
end
```

It is possible to add both default values and metadata to a parameter/species. In this case, first provide the default value, next the metadata. I.e. to in the above example set $kD$'s default value to $1.0$ we use
```@example dsl_advanced_metadata
two_state_system = @reaction_network begin
  @parameters kA [description="X's activation rate", bound=(0.01,10.0)] kD = 1.0 [description="X's deactivation rate"]
  (ka,kD), Xi <--> Xa
end
```

When designating metadata for species/parameters in `begin ... end` blocks the syntax changes slight. Here, a `,` must be inserted before the metadata (but after any potential default value). I.e. the previous example is rewritten as
```@example dsl_advanced_metadata
two_state_system = @reaction_network begin
  @parameters begin
    kA, [description="X's activation rate", bound=(0.01,10.0)]
    kD = 1.0, [description="X's deactivation rate"]
  end
  (ka,kD), Xi <--> Xa
end
```

Each metadata has its own getter functions. E.g. we can get the description of the parameter $pA$ using `getdescription` (here we use [system indexing](@ref ref) to access the parameter):
```@example dsl_advanced_metadata
getdescription(two_state_system.kA)
```

It is not possible for the user to directly designate their own metadata. These have to first be added to Catalyst. Doing so is somewhat involved, and described in detail [here](). A full list of metadata that can be used for species and/or parameters can be found [here](@ref ref).

### [Designating constant-valued/fixed species-parameters](@id dsl_advanced_options_constant_species)

Catalyst enables the designation of parameters as `constantspecies`. These can be used as species in reactions, however, their values are not changed by the reaction and remain constant throughout the simulation (unless changed by e.g. the [occurrence of a callback]@ref advanced_simulations_callbacks). Practically, this is done by setting the parameter's `isconstantspecies` metadata to `true`. Here, we create a simple reaction where the species $X$ is converted to $Xᴾ$ at rate $k$. By designating $X$ as a constant species parameter, we ensure that its quantity is unchanged by the occurrence of the reaction.
```@example dsl_advanced_constant_species
using Catalyst # hide
rn = @reaction_network begin
  @parameters X [isconstantspecies=true]
  k, X --> Xᴾ
end
```
We can confirm that $X$ is the only species of the system:
```@example dsl_advanced_constant_species
species(rn)
```
Here, the produced model is actually identical to if $X$ had simply been put as a parameter in the reaction's rate:
```@example dsl_advanced_constant_species
rn = @reaction_network begin
  k*X, 0 --> Xᴾ
end
```

A common use-case for constant species are when modelling systems where some species are present in such surplus that their amounts the reactions' effect on it is negligible. A system which is commonly modelled this way is [the Brusselator](https://en.wikipedia.org/wiki/Brusselator).

### [Specifying non-species variables](@id dsl_advanced_options_variables)
Chemical reaction network (CRN) models (which Catalyst creates) described how *species* are affected by the occurrence of reaction events. When they are converted to ODEs, the species are the variables of the system. However, Catalyst permits the creation of hybrid CRN models. These describe phenomenons which can only partially be modelled using CRNs. An example may be a bacterium. Here, we can use a CRN to model some internal system (e.g. controlling its growth rate). However, we might also want to model the bacterium's volume. Here, the volume cannot be considered a species (as it does not participate in reactions). Instead, we should model it as a normal variable. Here, Catalyst provides the `@variables` option for adding non-species variables to the system. E.g. to create a model where a single growth factor ($G$) is produced and degraded, and where we also have a single volume variables ($V$) we can use:
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

You can use the `variables` and `species` functions to retrieve a model's variables and species, respectively. The `unknown` function can be used to return both.

## [Setting reaction metadata](@id dsl_advanced_options_reaction_metadata)
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
nothing # hide
```
Here, since the networks' names are different:
```@example dsl_advanced_names
nameof(rn1) == nameof(rn2)
```
they are different
```@example dsl_advanced_names
rn1 == rn2
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

## [Creating observables](@id dsl_advanced_options_observables)

Sometimes, one might want to include observable variables. These are variables which values that can be computed directly from the systems species, parameters, and variables (rather than having their values implicitly given by an equation). Observables can be designated using the `@observables` option. Here, the `@observables` option is followed by a `begin ... end` block with one line for each observable. Each line first given the observable, followed by a `~` (*not* a `=`!), followed by an expression describing how to compute it.

Let us consider a model where two species ($X$ and $Y$) can bind to form a complex ($XY$, which also can dissociate back into $X$ and $Y$). If we wish to create a representation for the total amount of $X$ and $Y$ in the system, we can do this by creating observables $Xtot$ and $Ytot$:
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

Next, we can use [symbolic indexing](@ref ref) of our solution object, but with the observable as input. E.g. we can use 
```@example dsl_advanced_observables
sol[:Xtot]
```
to get a vector with teh value of $Xtot$ throughout the simulation. We can also use
```@example dsl_advanced_observables
using Plots
plot(sol; idxs = [:Xtot, :Ytot])
```
to plot the observables (rather than the species).

Observables can be defined using complicated expression containing species, parameters, and variables (but not other observables). in the following example (which uses a [parametric stoichiometry](@ref ref)) $X$ polymerises to form a complex $Xn$ containing $n$ copies of $X$. Here, we create an observable describing the total number of $X$ molecules in the system:
```@example dsl_advanced_observables
using Catalyst # hide
rn = @reaction_network begin
  @observables Xtot ~ X + n*XnXY  
  (kB,kD), n*X <--> Xn
end
nothing # hide
```
Note that, since we only have a single observable, the `begin .. end` block is not required and the observable can be declared directly after the `@observables` option.

[Metadata](@ref ref) can be supplied to an observable directly after the its declaration (but before its formula). If so, the metadata must be separated from the observable with a `,`, and the observable plus the metadata encapsulated by `()`. E.g. to add a [description metadata](@ref ref) to our observable we can do
```@example dsl_advanced_observables
using Catalyst # hide
rn = @reaction_network begin
  @observables (Xtot, [description="The total amount of X in the system."]) ~ X + n*XnXY  
  (kB,kD), n*X <--> Xn
end
nothing # hide
```

Observables are by default considered [variables](@ref dsl_advanced_options_variables) (not species). To designate them as a species, they can be pre-declared using the `@species` option. I.e. Here $Xtot$ becomes a species:
```@example dsl_advanced_observables
using Catalyst # hide
rn = @reaction_network begin
  @species Xtot(t)
  @observables (Xtot, [description="The total amount of X in the system."]) ~ X + n*XnXY  
  (kB,kD), n*X <--> Xn
end
nothing # hide
```

Some final notes regarding observables:
- The right-hand side of the observable declaration must contain a single symbol only (with the exception of metadata, which can also be supplied).
- All symbols appearing on the left-hand side must be declared elsewhere within the `@reaction_network` call (either by being part of a reaction, or through the `@species`, `@parameters`, or `@variables` options).
- Observables may not depend on other observables.
- Observables have their [dependent variables](@ref ref) automatically assigned as the union of the dependent variables of the species and variables on which it depends.

## [Creating events](@id dsl_advanced_options_events)

## [Specifying non-time independent variables](@id dsl_advanced_options_ivs)

As [described elsewhere](@ref ref), Catalyst's `ReactionSystem` models depends on on time independent variable, and potentially one or more spatial independent variables. By default, the independent variable `t` is used. We can declare another independent variable (which is automatically used as the default one) using teh `@ivs` option. E.g. to use `s` instead of `t` we can use
```@example dsl_advanced_ivs
using Catalyst # hide
rn = @reaction_network my_network begin
  @ivs s
  (ka,kD), Xi <--> Xa
end
nothing # hide
```
We can confirm that $Xi$ and $Xa$ depend on `s` (and not `t`):
```@example dsl_advanced_ivs
species(rn)
```

It is possible to designate several independent variables using `@ivs`. If so, the first one is considered the default, time, variable, while the following one is considered spatial variables. If we want some species to  depend on a non-time independent variable, this has to be explicitly declared:
```@example dsl_advanced_ivs
using Catalyst # hide
rn = @reaction_network my_network begin
  @ivs t x
  @species Y(s)
  (p1,d1), 0 <--> X
  (p2,d2), 0 <--> Y
end
species(rn)
```
Finally, it is possible to have species which depends on several independent variables:
```@example dsl_advanced_ivs
using Catalyst # hide
rn = @reaction_network my_network begin
  @ivs t x
  @species Xi(t,x) Xa(t,x)
  (ka,kD), Xi <--> Xa
end
species(rn)
```

!!! note
    Setting spatial independent variables is primarily intended for modelling of spatial systems on continuous domains.  Catalyst's support for this is currently under development. Hence, the utility of specifying indepdent variables is currently limited.