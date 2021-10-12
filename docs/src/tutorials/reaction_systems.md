# Symbolic Reaction Systems
While the DSL provides a simple interface for creating `ReactionSystem`s, it can
often be convenient to build or augment a [`ReactionSystem`](@ref)
programmatically. In this tutorial we show how to build the repressilator model
of the [Using Catalyst](@ref) tutorial directly using symbolic variables, and
then summarize the basic API functionality for accessing information stored
within `ReactionSystem`s.

## Directly Building the Repressilator with `ReactionSystem`s
We first load Catalyst
```julia
using Catalyst
```
and then define symbolic variables for each parameter and species in the system
(the latter corresponding to a `variable` or `state` in ModelingToolkit
terminology)
```julia
@parameters α, K, n, δ, γ, β, μ
@variables t, m₁(t), m₂(t), m₃(t), P₁(t), P₂(t), P₃(t)
```
*Note, each species is declared as a variable that is a function of time!*

Next we specify the chemical reactions that comprise the system using Catalyst
[`Reaction`](@ref)s
```julia
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
```
Here we use `nothing` where the DSL used ``\varnothing``. Finally, we are ready
to construct our [`ReactionSystem`](@ref) as
```julia
@named repressilator = ReactionSystem(rxs, t)
```
giving
```
Model repressilator with 15 equations
States (6):
  m₁(t)
  m₂(t)
  m₃(t)
  P₁(t)
  P₂(t)
  P₃(t)
Parameters (7):
  α
  K
  n
  δ
  γ
  β
  μ
```
Notice, the model is named `repressilator`. A name must always be specified when
directly constructing a `ReactionSystem` (the DSL will auto-generate one if
left out). Using `@named` when constructing a `ReactionSystem` causes the name
of the system to be the same as the name of the variable storing the system.
Alternatively, one can use the `name=:repressilator` keyword argument to the
`ReactionSystem` constructor.

We can check that this is the same model as the one we defined via the DSL as
follows (this requires that we use the same names for rates, species and the
system)
```julia
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
end α K n δ γ β μ
repressilator == repressilator2
```
which gives 
```
true
```

For more options in building `ReactionSystem`s, see the [`ReactionSystem`](@ref) API docs.

## More General `Reaction`s
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
@parameters α, β
@variables t, A(t), B(t)
rx = Reaction(α+β*t*A, [A], [B])
```
[See the FAQs](@ref user_functions) for info on using general user-specified
functions for the rate constant.

## [`ReactionSystem`](@ref) Accessors

A [`ReactionSystem`](@ref) is an instance of a
`ModelingToolkit.AbstractTimeDependentSystem`, and has a number of fields that
can be accessed using the [Catalyst.jl API](@ref) or the [ModelingToolkit.jl
Abstract System Interface](https://mtk.sciml.ai/dev/basics/AbstractSystem/).
Below we overview these components.

There are three basic sets of convenience accessors that will return information
either from a top-level system, the top-level system and all sub-systems that
are also `ReactionSystem`s (i.e. the full reaction-network), or the top-level
system and all sub-systems (i.e. the full model). To retrieve info from just a
base [`ReactionSystem`](@ref) `rn`, ignoring sub-systems of `rn`, one can use
the ModelingToolkit accessors:

* `get_states(rn)` is a vector that collects all the species defined within
  `rn`. 
* `get_ps(rn)` is a vector that collects all the parameters defined within
  reactions in `rn`. 
* `get_eqs(rn)` is a vector that collects all the [`Reaction`](@ref)s defined
  within `rn`.
* `get_iv(rn)` is the independent variable used in the system (usually `t` to
  represent time).
* `get_systems(rn)` is a vector of all sub-systems of `rn`.
* `get_defaults(rn)` is a dictionary of all the default values for parameters
  and species in `rn`.

These accessors do not allocate, directly accessing internal fields of the
`ReactionSystem`.

To retrieve the corresponding information from the full reaction network
represented by a system `rn`, which corresponds to information within both `rn`
and all sub-systems of type `ReactionSystem`, one can call:

* [`species(rn)`](@ref) is a vector collecting all the chemical species within
  the system and any sub-systems that are also `ReactionSystems`. 
* [`reactionparams(rn)`](@ref) is a vector of all the parameters
  within the system and any sub-systems that are also `ReactionSystem`s.
* [`reactions(rn)`](@ref) is a vector of all the `Reaction`s within the system
  and any sub-systems that are also `ReactionSystem`s.

These accessors will allocate unless there are no subsystems. In the latter case
they are equivalent to the corresponding `get_*` functions.

Finally, as some sub-systems may be other system types, for example specifying
algebraic constraints with a `NonlinearSystem`, it can also be convenient to
collect all state variables (e.g. species and algebraic variables) and such. The
following ModelingToolkit functions provide this information

* `ModelingToolkit.states(rn)` returns all species *and variables*
  across the system and *all sub-systems*.
* `ModelingToolkit.parameters(rn)` returns all parameters across the
  system and *all sub-systems*.
* `ModelingToolkit.equations(rn)` returns all [`Reaction`](@ref)s and
  all `Equations` defined across the system and *all sub-systems*.

`states` and `parameters` should be assumed to always allocate, while
`equations` will allocate unless there are no subsystems. In the latter case
`equations` is equivalent to `get_eqs`.

Empty `ReactionSystem`s can be generated via [`make_empty_network`](@ref) or
[`@reaction_network`](@ref) with no arguments (giving one argument to the latter
will specify a system name). `ReactionSystem`s can be programmatically extended
using [`addspecies!`](@ref), [`addparam!`](@ref), [`addreaction!`](@ref),
[`@add_reactions`](@ref), or composed using [`ModelingToolkit.extend`](@ref) and
[`ModelingToolkit.compose`](@ref).

### [`Reaction`](@ref) fields

Each `Reaction` within `reactions(rn)` has a number of subfields. For `rx` a
`Reaction` we have:
* `rx.substrates`, a vector of ModelingToolkit expressions storing each
  substrate variable.
* `rx.products`, a vector of ModelingToolkit expressions storing each product
  variable.
* `rx.substoich`, a vector storing the corresponding integer stoichiometry of
  each substrate species in `rx.substrates`.
* `rx.prodstoich`, a vector storing the corresponding integer stoichiometry of
  each product species in `rx.products`.
* `rx.rate`, a `Number`, `ModelingToolkit.Sym`, or ModelingToolkit expression
  representing the reaction rate. E.g., for a reaction like `k*X, Y --> X+Y`,
  we'd have `rate = k*X`.
* `rx.netstoich`, a vector of pairs mapping the ModelingToolkit expression for
  each species that changes numbers by the reaction to how much it changes. E.g.,
  for `k, X + 2Y --> X + W`, we'd have `rx.netstoich = [Y(t) => -2, W(t) => 1]`.
* `rx.only_use_rate`, a boolean that is `true` if the reaction was made with
  non-filled arrows and should ignore mass action kinetics. `false` by default.

