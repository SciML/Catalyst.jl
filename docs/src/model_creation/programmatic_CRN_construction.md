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
While programmatically created models can also have their parameters and species designated using `Symbol`s, the reverse is not possible for DSL-created models. Here, the symbolic variables are never explicitly declared, and thus not available for designating their values. Symbolic designation can be enabled for DSL-created models [by using `@unpack`](@ref programmatic_CRN_construction_symbolic_representation_unpack).

Elsewhere, we also describe how e.g. `ODEProblem`s and [simulations solutions can be queried for the values of model quantities](@ref simulation_structure_interfacing). There, we use use `Symbol`s to represent model quantities, however, symbolic variables (when available) can again be used. E.g. we can use
```@example programmatic_2
sol[X1]
```
to retrieve $X$'s value across the simulation.

### [Using `@unpack` to extract symbolic variables from `ReactionSystem`s](@id programmatic_CRN_construction_symbolic_representation_unpack)
Let us consider a simple birth-death model created using the DSL:
```@example programmatic_3
using Catalyst # hide
bd_model = @reaction_network begin
    (p,d), 0 <--> X
end
nothing # hide
```
Since we have not explicitly declared `p`, `d`, and `X` using `@parameters` and `@species`, we cannot represent these symbolically (only using `Symbol`s). If we wish to do so, however, we can fetch these into our current scope using the `@unpack` macro:
```@example programmatic_3
@unpack p, d, X = bd_model
nothing # hide
```
This lists first the quantities we wish to fetch (does not need to be the model's full set of parameters and species), then `=`, followed by the model name. `p`, `d` and `X` are now symbolic variables in the current scope, just as if they had been declared using `@parameters` or `@species`. We can confirm this:
```@example programmatic_3
X
```
Next, we can now use these to e.g. designate initial conditions and parameter values for model simulations:
```@example programmatic_3
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

### [Designating parameter types](@id programmatic_CRN_construction_symbolic_variables_options_types)
Sometimes it is desired to designate that a parameter should have a specific [type](@ref ref). When supplying this parameters value to e.g. an `ODEProblem`, that parameter will then be restricted to that specific type. Designating a type is done by appending the parameter with `::` followed by its type. E.g. to specify that a parameter `n` should be an `Int64` we do:
```@example programmatic_4
@parameters n(t)::Int64
nothing # hide
```

If a parameter have a type, a metadata, and a default value, they are designated in the following order:
```@example programmatic_4
TODO: @parameters n(t)::Int64
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


## [Additional options for declaration of `Reaction`s](@id programmatic_CRN_construction_reactions_options)

## [Creation of `Reaction`s using the`@reaction` ](@id programmatic_CRN_construction_reaction_macro)

## [Additional options for programmatic model creation](@id programmatic_CRN_construction_additional_options)



While the DSL provides a simple interface for creating `ReactionSystem`s, it can
often be convenient to build or augment a [`ReactionSystem`](@ref)
programmatically. In this tutorial we show how to build the repressilator model
of the [Introduction to Catalyst](@ref introduction_to_catalyst) tutorial directly using symbolic variables, and
then summarize the basic API functionality for accessing information stored
within `ReactionSystem`s.

## Directly building the repressilator with `ReactionSystem`s
We first load Catalyst
```@example ex
using Catalyst
```
and then define symbolic variables for each parameter and species in the system
(the latter corresponding to a `variable` or `unknown` in ModelingToolkit
terminology)
```@example ex
t = default_t()
@parameters α K n δ γ β μ
@species m₁(t) m₂(t) m₃(t) P₁(t) P₂(t) P₃(t)
nothing    # hide
```
Note: each species is declared as a function of time. Here, we first import the *time independent variable*, and stores it in `t`, using `t = default_t()`, and then use it to declare out species.

!!! note
       For users familiar with ModelingToolkit, chemical species must be declared
       via the `@species` macro, and not the `ModelingToolkit.@variables` macro.
       `@species` wraps `@variables`, adding additional metadata to the symbolic
       variables that represent species which is used internally in Catalyst.

Next, we specify the chemical reactions that comprise the system using Catalyst
[`Reaction`](@ref)s
```@example ex
rxs = [Reaction(hillr(P₃,α,K,n), nothing, [m₁]),
       Reaction(hillr(P₁,α,K,n), nothing, [m₂]),
       Reaction(hillr(P₂,α,K,n), nothing, [m₃]),
       Reaction(δ, [m₁], nothing),
       Reaction(γ, nothing, [m₁]),
       Reaction(δ, [m₂], nothing),
       Reaction(γ, nothing, [m₂]),
       Reaction(δ, [m₃], nothing),
       Reaction(γ, nothing, [m₃]),
       Reaction(β, [m₁], [m₁,P₁]),
       Reaction(β, [m₂], [m₂,P₂]),
       Reaction(β, [m₃], [m₃,P₃]),
       Reaction(μ, [P₁], nothing),
       Reaction(μ, [P₂], nothing),
       Reaction(μ, [P₃], nothing)]
nothing    # hide
```
Here we use `nothing` where the DSL used ``\varnothing``. Finally, we are ready
to construct our [`ReactionSystem`](@ref) as
```@example ex
@named repressilator = ReactionSystem(rxs, t)
nothing     # hide
```
Notice, the model is named `repressilator`. A name must always be specified when
directly constructing a `ReactionSystem` (the DSL will auto-generate one if left
out). Using `@named` when constructing a `ReactionSystem` causes the name of the
system to be the same as the name of the variable storing the system.
Alternatively, one can use the `name = :repressilator` keyword argument to the
`ReactionSystem` constructor.

!!! warn
       All `ReactionSystem`s created via the symbolic interface (i.e. by calling `ReactionSystem` with some input, rather than using `@reaction_network`) are not marked as complete. To simulate them, they must first be marked as *complete*, indicating to Catalyst and ModelingToolkit that they represent finalized models. This can be done using the `complete` function, i.e. by calling `repressilator = complete(repressilator)`. An expanded description on *completeness* can be found [here](@ref completeness_note).

We can check that this is the same model as the one we defined via the DSL as
follows (this requires that we use the same names for rates, species and the
system)
```@example ex
repressilator2 = @reaction_network repressilator begin
    hillr(P₃,α,K,n), ∅ --> m₁
    hillr(P₁,α,K,n), ∅ --> m₂
    hillr(P₂,α,K,n), ∅ --> m₃
    (δ,γ), m₁ <--> ∅
    (δ,γ), m₂ <--> ∅
    (δ,γ), m₃ <--> ∅
    β, m₁ --> m₁ + P₁
    β, m₂ --> m₂ + P₂
    β, m₃ --> m₃ + P₃
    μ, P₁ --> ∅
    μ, P₂ --> ∅
    μ, P₃ --> ∅
end
repressilator == repressilator2
```

For more options in building `ReactionSystem`s, see the [`ReactionSystem`](@ref)
API docs. For a more extensive example of how to programmatically create a
`ReactionSystem`, see the [Smoluchowski Coagulation Equation example](@ref
smoluchowski_coagulation_equation).

## More general `Reaction`s
In the example above all the specified `Reaction`s were first or zero order. The
three-argument form of `Reaction` implicitly assumes all species have a
stoichiometric coefficient of one, i.e. for substrates `[S₁,...,Sₘ]` and
products `[P₁,...,Pₙ]` it has the possible forms
```julia
# rate, S₁ + ... + Sₘ --> P₁ + ... + Pₙ
Reaction(rate, [S₁,...,Sₘ], [P₁,...,Pₙ])

# rate, S₁ + ... + Sₘ --> ∅
Reaction(rate, [S₁,...,Sₘ], nothing)

# rate, ∅ --> P₁ + ... + Pₙ
Reaction(rate, nothing, [P₁,...,Pₙ])
```
To allow for other stoichiometric coefficients we also provide a five argument
form
```julia
# rate, α₁*S₁ + ... + αₘ*Sₘ --> β₁*P₁ + ... + βₙ*Pₙ
Reaction(rate, [S₁,...,Sₘ], [P₁,...,Pₙ], [α₁,...,αₘ], [β₁,...,βₙ])

# rate, α₁*S₁ + ... + αₘ*Sₘ --> ∅
Reaction(rate, [S₁,...,Sₘ], nothing, [α₁,...,αₘ], nothing)

# rate, ∅ --> β₁*P₁ + ... + βₙ*Pₙ
Reaction(rate, nothing, [P₁,...,Pₙ], nothing, [β₁,...,βₙ])
```
Finally, we note that the rate constant, `rate` above, does not need to be a
constant or fixed function, but can be a general symbolic expression:
```julia
t = default_t()
@parameters α, β
@species A(t), B(t)
rx = Reaction(α + β*t*A, [A], [B])
```
[See the FAQs](@ref user_functions) for info on using general user-specified
functions for the rate constant.

## The `@reaction` macro for constructing `Reaction`s
In some cases one wants to build reactions incrementally, as in the
repressilator example, but it would be nice to still have a short hand as in the
[`@reaction_network`](@ref) DSL. In this case one can construct individual
reactions using the [`@reaction`](@ref) macro.

For example, the repressilator reactions could also have been constructed like
```julia
t = default_t()
@species P₁(t) P₂(t) P₃(t)
rxs = [(@reaction hillr($P₃,α,K,n), ∅ --> m₁),
       (@reaction hillr($P₁,α,K,n), ∅ --> m₂),
       (@reaction hillr($P₂,α,K,n), ∅ --> m₃),
       (@reaction δ, m₁ --> ∅),
       (@reaction γ, ∅ --> m₁),
       (@reaction δ, m₂ --> ∅),
       (@reaction γ, ∅ --> m₂),
       (@reaction δ, m₃ --> ∅),
       (@reaction γ, ∅ --> m₃),
       (@reaction β, m₁ --> m₁ + P₁),
       (@reaction β, m₂ --> m₂ + P₂),
       (@reaction β, m₃ --> m₃ + P₃),
       (@reaction μ, P₁ --> ∅),
       (@reaction μ, P₂ --> ∅),
       (@reaction μ, P₃ --> ∅)]
@named repressilator = ReactionSystem(rxs, t)
```
Note, there are a few differences when using the `@reaction` macro to specify
one reaction versus using the full `@reaction_network` macro to create a
`ReactionSystem`. First, only one reaction (i.e. a single forward arrow type)
can be used, i.e. reversible arrows like `<-->` will not work (since they define
more than one reaction). Second, the `@reaction` macro does not have an option for designating what should be considered a species or parameter, and instead assumes that any symbol that appears as either a substrate or a product is a species, and everything else (including stoichiometric coefficients) are parameters. As such, the following are equivalent
```julia
rx = @reaction hillr(P,α,K,n), A --> B
```
is equivalent to
```julia
t = default_t()
@parameters P α K n
@variables A(t) B(t)
rx = Reaction(hillr(P,α,K,n), [A], [B])
```
Here `(P,α,K,n)` are parameters and `(A,B)` are species.

This behavior is the reason that in the repressilator example above we
pre-declared `(P₁(t),P₂(t),P₃(t))` as variables, and then used them via
interpolating their values into the rate law expressions using `$` in the macro.
This ensured they were properly treated as species and not parameters. See the
[`@reaction`](@ref) macro docstring for more information.

## Basic querying of `ReactionSystems`

The [Catalyst.jl API](@ref) provides a large variety of functionality for
querying properties of a reaction network. Here we go over a few of the most
useful basic functions. Given the `repressillator` defined above we have that
```@example ex
species(repressilator)
```
```@example ex
parameters(repressilator)
```
```@example ex
reactions(repressilator)
```

We can test if a `Reaction` is mass action, i.e. the rate does not depend on `t`
or other species, as
```@example ex
# Catalyst.hillr(P₃(t), α, K, n), ∅ --> m₁
rx1 = reactions(repressilator)[1]
ismassaction(rx1,repressilator)
```
while
```@example ex
# δ, m₁ --> ∅
rx2 = reactions(repressilator)[4]
ismassaction(rx2,repressilator)
```
Similarly, we can determine which species a reaction's rate law will depend on
like
```@example ex
rn = @reaction_network begin
       k*W, 2X + 3Y --> 5Z + W
     end
dependents(reactions(rn)[1], rn)
```
Basic stoichiometry matrices can be obtained from a `ReactionSystem` as
```@example ex
substoichmat(repressilator)
```
```@example ex
prodstoichmat(repressilator)
```
```@example ex
netstoichmat(repressilator)
```
Here the ``(i,j)`` entry gives the corresponding stoichiometric coefficient
of species ``i`` for reaction ``j``.

Finally, we can directly access fields of individual reactions like
```@example ex
rx1.rate
```
```@example ex
rx1.substrates
```
```@example ex
rx1.products
```
```@example ex
rx1.substoich
```
```@example ex
rx1.prodstoich
```
```@example ex
rx1.netstoich
```

See the [Catalyst.jl API](@ref) for much more detail on the various querying and
analysis functions provided by Catalyst.
