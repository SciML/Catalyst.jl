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
Previously, we have described how to [simulate Catalyst models declared via the DSL](@ref ref). Here, we used `Symbol`s (e.g. `:k1` and `:X1`) to designate parameters and species' values. However, in programmatic model creation we have explicitly [declared symbolic variables](@ref programmatic_CRN_construction_symbolic_variables) corresponding to our parameters and species. Here, it is instead possible to use these to designate our simulation conditions. Below we utilise this to declare initial conditions and parameters, and then using this to 

## [Additional options for declaration of symbolic variables](@id programmatic_CRN_construction_symbolic_variables_options)

## [Additional options for declaration of `Reaction`s](@id programmatic_CRN_construction_reactions_options)

## [Creation of `Reaction`s using the`@reaction` ](@id programmatic_CRN_construction_reaction_macro)

## [Additional options for programmatic model creation](@id programmatic_CRN_construction_additional_options)


!!! note
       While `Reaction`s are not explicitly created when models are created via the DSL, these structures are still declared and stored internally. When calling the 

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
