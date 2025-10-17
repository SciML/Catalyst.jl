# [Programmatic Construction of Symbolic Reaction Systems](@id programmatic_CRN_construction)

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
repressilator = complete(repressilator)
nothing     # hide
```

Notice, the model is named `repressilator`. A name must always be specified when
directly constructing a `ReactionSystem` (the DSL will auto-generate one if left
out). Using `@named` when constructing a `ReactionSystem` causes the name of the
system to be the same as the name of the variable storing the system.
Alternatively, one can use the `name = :repressilator` keyword argument to the
`ReactionSystem` constructor.

!!! warning
       All `ReactionSystem`s created via the symbolic interface (i.e. by calling `ReactionSystem` with some input, rather than using `@reaction_network`) are not marked as complete. To simulate them, they must first be marked as *complete*, indicating to Catalyst and ModelingToolkit that they represent finalized models. This can be done using the `complete` function, as above. An expanded description on *completeness* can be found [here](@ref completeness_note).

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

The [Catalyst.jl API](@ref api) provides a large variety of functionality for
querying properties of a reaction network. Here we go over a few of the most
useful basic functions. Given the `repressilator` defined above we have that

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

See the [Catalyst.jl API](@ref api) for much more detail on the various querying and
analysis functions provided by Catalyst.
