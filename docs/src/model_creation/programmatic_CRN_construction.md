# [Programmatic Model Construction](@id programmatic_CRN_construction)
While the [Catalyst DSL](@ref ref) provides a simple interface for creating Catalyst [`ReactionSystem`](@ref) models using the `@reaction_network` macro, Catalyst also permits models to be created *programmatically*. Here, instead of declaring the model in one go, the models components are instead declared one by one, and then assembled to a full model. While all Catalyst features are accessible through both approaches, programmatic model creation can sometimes give the user more control over the modelling procedure. The DSL instead provides a more concise and human-readable notation. 

Briefly, programmatic model creation consists of three steps:
1. Declaration of *symbolic variables* representing model species and parameters.
2. The assembly of these into model reactions.
3. The creation of a `ReactionSystem` model from the reactions.

The creation of more complicated models (containing e.g. [events](@ref ref) or [non-reaction-based equations](@ref ref)) is also possible, and will be described later on.

Finally, this tutorial will discuss the concepts of model *names* and *completeness*. While also relevant to DSL-created models, only programmatic modelling requires the users to handle these explicitly, and hence they will be given additional attention here.

!!! note
       Programmatic creation of `ReactionSystem`s is highly related to [programmatic creation of ModelingToolkit.jl models](@ref ref) (from which this approach is derived). Users familiar with ModelingToolkit will be able to reuse this knowledge when declaring Catalyst models programmatically (but should still read the below tutorial, to take note of a few crucial differences).

### [Basic example](@id programmatic_CRN_construction_basic_example)
Before describing the three steps one by one, we will give a simple example where we create a [birth-death model](@ref ref). First we declare the [time independent variable](@ref ref), the species, and the parameters:
```@example programmatic_1
using Catalyst
t = default_t()
@species X(t)
@parameters p d
nothing # hide
```
Next we assembly the production and degradation reactions:
```@example programmatic_1
rxs = [
       Reaction(p, [], [X]),
       Reaction(d, [X], []),
]
nothing # hide
```
Finally, we use the reactions and time independent variable as input to the `ReactionSystem` structure constructor:
```@example programmatic_1
@named bd_model = ReactionSystem(rxs, t)
bd_model = complete(bd_model)
```
He we also pre-append the model declaration with the `@named` macro (to [give the model a name](@ref ref)) and use `complete` to [mark our model as complete](@ref ref).

## [Declaration of symbolic variables](@id programmatic_CRN_construction_symbolic_variables)
Internally, all Catalyst quantities (parameters, species, and [variables](@ref ref)) are represented as *symbolic variables*. This enables them to [form algebraic expressions](@ref ref), something Catalyst uses to e.g. build differential equations.

In practise, symbolic variables are declared using the `@parameters` and `@species` macros. In this tutorial we will create a simple [two-state model](@ref ref), for which we need to declare two parameters ($k1$ and $k2$). We do this by simply listing these after the `@parameters` macro:
```@example programmatic_2
using Catalyst # hide
@parameters k1 k2
```
Here, `@parameters` create two [new Julia variables](@ref ref), `k1` and `k2`, and stores the symbolic variables `k1` and `k2` in these. This creates the variables `k1` and `k2` in the current scope, something we can confirm:
```@example programmatic_2
k1
```
!!! warn
       Since `@parameters` creates the variables `k1` and `k2` in the current scope, this will override any previous variables with these names. Furthermore, if you store new values in the variables (e.g. `k1 = 1.0`), this will override the symbolic variables stored within these, preventing them from being used to e.g. build `ReactionSystem`s.

Next, we declare the species. Species are not constants, but [functions of the time independent variable](@ref ref). Hence, this must be declared first:
```@example programmatic_2
t = default_t()
```
While [non-time independent variables is possible](@ref), and independent variables can be declared without using the `default_t` function, the approach is the simplest (and by far most common) one. Next, the species can be declared using similar notation as the parameters (however, by using the `@species` macro, and by appending the time dependency `(t)` to each species):
```@example programmatic_2
@species X1 X2
```

Here we have described the basics of symbolic variable declarations. A more throughout overview (describing additional options) can be found [here](@ref ref).

!!! note
       For users familiar with [ModelingToolkit](@ref ref), the `@species` macro have a very similar functionality to the `@variables` macro. However, `@species` appends additional metadata that the variable is a *species*, which is required for it to be used as a reaction reactant.

!!! note
       While symbolic variables are not explicitly created when models are created via the DSL, these are still declared and stored internally. When calling the [`parameters`](@ref) or [`species](@ref) functions on a DSL-created `ReactionSystem`, you will symbolic variables of the same form as we here have created programmatically.

## [Creation of `Reaction`s](@id programmatic_CRN_construction_reactions)
In the next step, we can assembly our models reactions. In DSL-based modelling, these are listed within the `@reaction_network` macro. For programmatic modelling, these are instead declared using the `Reaction` structure's constructor:
```@example programmatic_2
rxs = [
       Reaction(k1, [X1], [X2]),
       Reaction(k2, [X2], [X1])
]
```
Here, `Reaction` takes three arguments:
1. The rate (in these cases a single parameter, however, [other types of rates are possible](@ref ref)).
2. A vector with the reaction's substrates (in these case, both reactions have a single substrate).
3. A vector with the reaction's products (in these case, both reactions have a single product). 

Just like [when the DSL is used], more complicated reactions (featuring e.g. [](@ref ref), [](@ref ref), [](@ref ref), and [](@ref ref)) are possible. How to create such reactions is described [here](@ref ref).

!!! note
       While `Reaction`s are not explicitly created when models are created via the DSL, these are still declared and stored internally. When calling the [`reactions`](@ref) function on a DSL-created `REactionSystem`, you will receive `Reaction`s of the same form as we here have created programmatically.

## [Creation of `ReactionSystem`s](@id programmatic_CRN_construction_reactionsystems)
Finally, we can use our `Reaction` vector as input to the Catalyst's `ReactionSystem` constructor. In addition to these, we need two additional arguments:
1. The independent variable (typically time) must be explicitly provided.
2. A [model name](@ref ref) bust be specified.

```@example programmatic_2
two_state_model = ReactionSystem(rxs, t; name = :two_state_model)
```
Here, while the [DSL-created models also have names](@ref ref), these must not be explicitly declared on creation. Model declaration can be pre-appending with the `@named` macro: 
```@example programmatic_2
@named two_state_model = ReactionSystem(rxs, t)
nothing # hide
```
This will automatically take name of the variable in which the model is stored, and use that as the model's name (here `two_state_model`). We can use the `getname` function to confirm this:
```@example programmatic_2
getname(two_state_model)
```
There exists a few additional options that can be supplied to the `ReactionSystem` constructor. These are described in the `ReactionSystem`s API entry, which can be found [here](@ref ref).

While we now have programmatically created a `ReactionSystem`, there are two final points (described in the next two sections) we should consider before using it for e.g. simulations.

## [System completeness](@id programmatic_CRN_construction_completeness)
`ReactionSystem` models created in Catalyst can either be *complete* or *incomplete*. This is primarily important for two reasons:
- Only complete models can be used as inputs to simulations or certain tools for model analysis.
- Only incomplete models can be [composed with other models to form hierarchical models](@ref ref).

A model's completeness depends on how it was created:
- Models created using the `@reaction_network` DSL are *complete*.
- To create *incomplete models using the DSL*, use the [`@network_component` macro](@ref ref).
- Models created programmatically are *incomplete*.
- Models generated through the `compose` (and `extend`) functions are *incomplete*.

Here, even if we do not intend to use our two-state model for hierarchical modelling or not, since it was created programmatically it is incomplete. We can confirm this using the `Catalyst.iscomplete` function:
```@example programmatic_2
Catalyst.iscomplete(two_state_model)
```
Since `two_state_model` is an incomplete model, we cannot simulate it. To do so, we need to use it to generate a new, complete, model. We do so through the `complete` function (and then store the new model in the same variable in which we stored the incomplete one). Finally, we confirm that the new model is complete
```@example programmatic_2
two_state_model = complete(two_state_model)
Catalyst.iscomplete(two_state_model)
```

We have now created a two-state model programmatically. For all purposes, there is not difference between this model and the same model as created with the DSL. Below we confirm this:
```@example programmatic_2
two_state_model_dsl = @reaction_network two_state_model begin
       k1, X1 --> X2
       k2, X2 --> X1
end
two_state_model == two_state_model_dsl
```

## [Symbolic designation of model quantities and simulation of programmatic models](@id programmatic_CRN_construction_symbolic_representation)
Previously, we have described how to [simulate Catalyst models declared via the DSL](@ref ref). Here, we used `Symbol`s (e.g. `:k1` and `:X1`) to designate parameters and species' values. However, in programmatic model creation we have explicitly [declared symbolic variables](@ref programmatic_CRN_construction_symbolic_variables) corresponding to our parameters and species. Here, it is instead possible to use these to designate our simulation conditions. Below we utilise this to declare initial conditions and parameters, use these as input to a `ODEProblem`, and then simulate it and plot the result:
```@example programmatic_2
using OrdinaryDiffEq, Plots
u0 = [X1 => 0.0, X2 => 2.0]
tspan = (0.0, 1.0)
ps = [k1 => 1.0, ks => 2.0]
oprob = ODEProblem(two_state_model, u0, tspan, ps)
sol = solve(oprob)
plot(sol)
```
While programmatically created models can also have their parameters and species designated using `Symbol`s, the reverse is not possible for DSL-created models. Here, the symbolic variables are never explicitly declared, and thus not available for designating their values. Symbolic designation can be enabled for DSL-created models [by using `@unpack`](@ref programmatic_CRN_construction_symbolics_and_DSL_unpack).

Elsewhere, we also describe how e.g. `ODEProblem`s and [simulations solutions can be queried for the values of model quantities](@ref simulation_structure_interfacing). There, we use use `Symbol`s to represent model quantities, however, symbolic variables (when available) can again be used. E.g. we can use
```@example programmatic_2
sol[X1]
```
to retrieve $X$'s value across the simulation.

## [Additional options for declaration of symbolic variables](@id programmatic_CRN_construction_symbolic_variables_options)
The declaration of symbolic variables for programmatic Catalyst modelling uses identical syntax as when [parameters/species are explicitly declared within the DSL](@ref ref), or as used within [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl). Here we will provide a brief summary iterating the same information.

### [Designating default values](@id programmatic_CRN_construction_symbolic_variables_options_defaults)
Species and parameters can be assigned default values. These permits the omission of their values from the initial condition and parameter values vector (in which case the default value is used instead). A default value is specified by simply a quantity's declaration by `=` followed by its default value. E.g. here we give the two-state model's species default initial condition values:
```@example programmatic_4
using Catalyst # hide
@species X1(t)=0.0 X2(t)=2.0
nothing # hide
```
A symbolic variables default value can be retrieved using the [`getdefault`](@ref ref) function:
```@example programmatic_4
getdefault(X1)
```

A more throughout description on how to work with default values, including how to e.g. set parametric initial conditions, is provided [when defaults are introduced in the context of DSL-based modelling](@ref ref).

### [Designating metadata](@id programmatic_CRN_construction_symbolic_variables_options_metadata)
Sometimes one might want to attach additional information to a symbolic variable. This is done through the addition of *metadata*. To attach metadata to a symbolic variable, simply follow its declaration by brackets (`[]`) enclosing a list of metadata and their values (separated by a `=`). E.g. here we attach [`description`](@ref ref) metadata to the two-state model's species:
```@example programmatic_4
@species X1(t) [description="Species X1"] X2(t) [description="Species X2"]
nothing # hide
```
A symbolic variable can have several metadata, and different sets of metadata can be supplied to each declared symbolic variable:
```@example programmatic_4
@parameters k1 [description="Parameter X1"]  k2 [description="Parameter X1", bounds=(1e-2, 1e2)]
nothing # hide
```

If a symbolic variable have both metadata and a default value, the default value is designated first, and the metadata second:
```@example programmatic_4
@species X1(t)=0.0 [description="Species X1"] X2(t)=2.0 [description="Species X2"]
nothing # hide
```

A list of commonly used metadata, and how to access them, can be found [here](https://docs.sciml.ai/ModelingToolkit/stable/basics/Variable_metadata/).

### [Designating parameter types](@id programmatic_CRN_construction_symbolic_variables_options_types)
Sometimes it is desired to designate that a parameter should have a specific [type](@ref ref). When supplying this parameters value to e.g. an `ODEProblem`, that parameter will then be restricted to that specific type. Designating a type is done by appending the parameter with `::` followed by its type. E.g. to specify that a parameter `n` should be an `Int64` we do:
```@example programmatic_4
@parameters n(t)::Int64
nothing # hide
```

If a parameter have a type, a metadata, and a default value, they are designated in the following order:
```@example programmatic_4
@parameters n(t)::Int64 = 2 ["Parameter n, which is an integer and defaults to the value 2."]
nothing # hide
```

### [Multi-line declarations](@id programmatic_CRN_construction_symbolic_variables_options_multiline)
Sometimes, when declaring a large number of quantities (or providing extensive additional information), readability can be improved through a multi-line statement. Here, the `@parameters`/`@species` macros are followed by a `begin ... end` block. Each line within the block contains the declaration of a single parameter/species. E.g. the parameters and species of or two-state model could have been declared using
```@example programmatic_4
@parameters begin
    k1
    k2
end
@species begin
    X1(t)
    X2(t)
end
nothing # hide
```

When using multi-line statements, default values, metadata, and types is appended *almost* identically as for the single-lines one. The one exception is that when metadata is added, this must be separated from the quantity by a `,`:
```@example programmatic_4
@parameters begin
    k1, [description = "Parameter k1"]
    k2, [description = "Parameter k2"]
end
nothing # hide
```

### [Vector-valued symbolic variables](@id programmatic_CRN_construction_symbolic_variables_options_vectors)
Sometimes, one wish to declare a large number of similar symbolic variables. E.g. if we have a system with ten species, each being produced at different rates, we could declare ten separate production parameters:
```@example programmatic_5
using Catalyst # hide
@parameters p1 p2 p3 p4 p5 p6 p7 p8 p9 p10
```
However, it is also possible to *declare a vector parameter* with ten different values:
```@example programmatic_5
@parameters p[1:10]
```

We can use this to create our two-state model, but instead of declaring `X1`, `X2`, `k1`, and `k2` as separate entities, we create them as two length-two vectors:
```@example programmatic_5
t = default_t() 
@species X[1:2](t)
@parameters k[1:2]
nothing # hide
```
Next, we can designate the individual parameters using e.g. `X[1]`
```@example programmatic_5
rxs = [
    Reaction(k[1], [X[1]], [X[2]]),
    Reaction(k[2], [X[2]], [X[1]]),
]
@named two_state_model = ReactionSystem(rxs, t)
```
Now, while we still can provide individual values to `X[1]`, `X[2]`, `k[1]`, and `k[2]`, we can also declare their values directly as vectors:
```@example programmatic_5
using OrdinaryDiffEq, Plots # hide
u0 = [X => [0.0, 2.0]]
tspan = (0.0, 1.0)
ps = [k => [1.0, 2.0]]
oprob = ODEProblem(two_state_model, u0, tspan, ps)
sol = solve(oprob)
plot(sol)
```


## [Additional options for declaration of `Reaction`s](@id programmatic_CRN_construction_reactions_options)
When describing the DSL, we also describe a range of [options for declaring various types of reactions](@ref ref). Each type of reaction that can be created using the DSL can also be created programmatically. Below, we briefly describe each case.

### [Reactions with non-unitary stoichiometries](@id programmatic_CRN_construction_reactions_options_stoichiometries)
Previously, we assumed that all reactions' substrates and products had stoichiometry 1. Other stoichiometries ([including decimal, parametric, and distributed, ones](@ref ref)) are possible. To designate this we provide two additional arguments to the `Reaction` constructor, designating the substrates and products stoichiometries respectively.

E.g. to model a simple dimerisation system (where two copies of $X$ dimerise to form $X2$, which then maybe dissociate back into two copies of $X$) we use
```@example programmatic_6
using Catalyst # hide
t = default_t() 
@species X(t) X2(t)
@parameters kB kD
rxs = [
    Reaction(kB, [X], [X2], [2], [1]),
    Reaction(kD, [X2], [X], [1], [2])
]
@named dimerisation_model = ReactionSystem(rxs, t)
```
Here, `Reaction(k, [X], [X2], [2], [1])` indicates that there are `2` copies of $X$ and one copy of $X2$ involved in the reaction.

If there are multiple substrates and/or products, the order of the stoichiometries must correspond to the order with which these occur. E.g. to create a reaction `k, 2X + Y --> X + 2Y` we would use
```@example programmatic_6
@parameters k
@species X(t) Y(t)
Reaction(k, [X, Y], [X, Y], [2, 1], [1, 2]),
```

### [Production and degradation reactions](@id programmatic_CRN_construction_reactions_options_production_and_degradation)
To designate [an absence of substrates and/or products](@ref ref), we simply give an empty vector. E.g. to create a [birth-death model](@ref ref) we use
```@example programmatic_6
@parameters p d
@species X(t)
rxs = [
    Reaction(p, [], [X]),
    Reaction(d, [X], [])
]
@named bd_model = ReactionSystem(rxs, t)
nothing # hide
```
If there also are non-unitary stoichiometries, we must supply empty stoichiometric vectors where appropriate. E.g. if for the previous model the reactions each involved two copies of $X$ (rather than a single one) we would use:
```@example programmatic_6
rxs = [
    Reaction(p, [], [X], [], [2]),
    Reaction(d, [X], [], [2], [])
]
```

### [Reactions with non-constant rates](@id programmatic_CRN_construction_reactions_options_rates)
Just like when creating models using the DSL, the rate can be any valid algebraic expressions (containing any combinations of parameters, species, and numeric values). E.g. to create a reaction that takes $X$ from its inactive form ($Xi$) to its activate form ($Xa$) catalysed by an enzyme ($E$) we would use:
```@example programmatic_6
@parameters k
@species Xi(t) Xa(t) E(t)
Reaction(k*E, [Xi], [Xa])
```

Here it is both possible to use known Julia functions:
```@example programmatic_6
Reaction(k*sqrt(E), [Xi], [Xa])
```
and ones declared by the user:
```@example programmatic_6
inhibition_function(x,p) = p/(p+x)
Reaction(inhibition_function(k,E), [Xi], [Xa])
```
Finally, Catalyst also pre-defined a few function commonly used in systems biology (listed [here](@ref ref)) which also can be used. E.g. to make the rate follow a [Michaelis-Menten function](https://en.wikipedia.org/wiki/Michaelis%E2%80%93Menten_kinetics) we use
```@example programmatic_6
@parameters v K
Reaction(mm(E,v,K), [Xi], [Xa])
```

### [Additional `Reaction` constructor arguments](@id programmatic_CRN_construction_reactions_options_)
The `Reaction` constructor accepts several additional optional arguments. E.g. to [disable to use of mass action to compute reaction propensities, and instead use the rate only], you can use the `only_use_rate = true` argument:
```@example programmatic_6
Reaction(k*E, [Xi], [Xa]; only_use_rate = true)
nothing # hide
```

To attach [metadata to a reaction](@ref ref), use the `metadata` argument, followed by a vector containing all metadata entries you wish to supply. Here, each metadata's name is designated using a `Symbol`, its value uses whichever form that metadata accepts, and these two are separated by a `=>`. E.g. to attach a description metadata to a reaction we use:
```@example programmatic_6
Reaction(k*E, [Xi], [Xa]; metadata = [:description => "The activation of species X"])
nothing # hide
```

## [Creation of `Reaction`s using the `@reaction`](@id programmatic_CRN_construction_reaction_macro)
Catalyst's DSL provides a simple way to declare reactions. By using the `@reaction` macro, it is possible to harness this during programmatic modelling. Here, `@reaction` accepts a single line of the same format as used with the DSL, and returns a single `Reaction` object (just as if it was created using the `Reaction` constructor).

E.g. here we create the reactions of our two-state model using the `@reaction` macro:
```@example programmatic_7
using Catalyst # hide
rxs = [
    @reaction k1, X1 --> X2,
    @reaction k2, X2 --> X1
]
```
We can now use these reactions as input to the `ReactionSystem` constructor to build our model:
```@example programmatic_7
t = default_t()
@named two_state_model = ReactionSystem(rxs, t)
nothing # hide
```

We note that when using `@reaction` it might not be necessary to declare the parameters and species using `@parameters` and `@species`. Instead, `@reaction` automatically infers and creates these from its input. Declaring these explicitly is still, however, possible, in which case they should be [interpolated into the `@reaction` macro](@ref programmatic_CRN_construction_symbolics_and_DSL_interpolation).

## [Working with symbolic variables and the DSL](@id programmatic_CRN_construction_symbolics_and_DSL)
The `@reaction` macro allowed us to use some of DSL-based modellings advantage (easy declaration of reactions) in programmatic modelling. Similarity, there are also ways to utilise concepts from programmatic modelling (like the declaration and use of symbolic variables) in DSL-based modelling. Below we briefly describe two of these.

### [Using `@unpack` to extract symbolic variables from `ReactionSystem`s](@id programmatic_CRN_construction_symbolics_and_DSL_unpack)
Let us consider a simple birth-death model created using the DSL:
```@example programmatic_8
using Catalyst # hide
bd_model = @reaction_network begin
    (p,d), 0 <--> X
end
nothing # hide
```
Since we have not explicitly declared `p`, `d`, and `X` using `@parameters` and `@species`, we cannot represent these symbolically (only using `Symbol`s). If we wish to do so, however, we can fetch these into our current scope using the `@unpack` macro:
```@example programmatic_8
@unpack p, d, X = bd_model
nothing # hide
```
This lists first the quantities we wish to fetch (does not need to be the model's full set of parameters and species), then `=`, followed by the model name. `p`, `d` and `X` are now symbolic variables in the current scope, just as if they had been declared using `@parameters` or `@species`. We can confirm this:
```@example programmatic_3
X
```
Next, we can now use these to e.g. designate initial conditions and parameter values for model simulations:
```@example programmatic_8
using OrdinaryDiffEq, Plots # hide
u0 = [X => 0.1]
tspan = (0.0, 10.0)
ps = [p => 1.0, d => 0.2]
oprob = ODEProblem(bd_model, u0, tspan, ps)
sol = solve(oprob)
plot(sol)
```

!!! warn
       Just like [when using `@parameters` and `@species`](@ref programmatic_CRN_construction_symbolic_variables), `@unpack` will overwrite any variables in the current scope which shares name with the imported quantities.

### [Interpolating variables into the DSL](@id programmatic_CRN_construction_symbolics_and_DSL_interpolation)

The DSL allows Julia variables to be interpolated for the network name, within
rate constant expressions, or for species/stoichiometry within reactions. Using
the lower-level symbolic interface we can then define symbolic variables and
parameters outside of the macro, which can then be used within expressions in
the DSL (see the [Programmatic Construction of Symbolic Reaction Systems](@ref programmatic_CRN_construction)
tutorial for details on the lower-level symbolic interface). For example,
```@example programmatic_9
using Catalyst
@parameters k α
@variables t
@species A(t)
spec = A
par = α
rate = k*A
name = :network
rn = @reaction_network $name begin
    $rate*B, 2*$spec + $par*B --> $spec + C
  end
```
As the parameters `k` and `α` were pre-defined and appeared via interpolation,
we did not need to declare them within the `@reaction_network` macro,
i.e. they are automatically detected as parameters:
```@example programmatic_9
parameters(rn)
```
as are the species coming from interpolated variables
```@example programmatic_9
species(rn)
```

!!! note
    When using interpolation, expressions like `2$spec` won't work; the
    multiplication symbol must be explicitly included like `2*$spec`.

## [Additional options for programmatic model creation](@id programmatic_CRN_construction_additional_options)
