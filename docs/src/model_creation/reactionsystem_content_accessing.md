# [Accessing model properties](@id model_accessing)
Catalyst is based around the creation analysis, and simulation of chemical reaction network models. Throughout the package, these models are stored in [`ReactionSystem](@ref) structures. Here, we will describe some basic functions for accessing their content. This includes e.g. retrieving lists of species, parameters, and reactions, that a model consists of. An extensive lists of relevant functions for working with `ReactionSystem` models can be found in the [API](@ref ref).

!!! warn
    Generally, the field of Julia structures can be accessed through `struct.field`. E.g. a simulation's time vector can be retrieved using `simulation.t`. While Catalyst `ReactionSystem`s are structured, one should *never* access their fields using this approach, but rather using the accessor functions described below and in the [API](@ref ref) (attempting to do so can yield unexpected behaviours). E.g. to retrieve the species of a `ReactionsSystem` `rs`, use `species(rs), *not* `rs.species`.

## [Direct accessing of symbolic model parameter and species](@id model_accessing_symbolic_variables)
Previously we have described how the parameters and species (and [variables](@ref ref)) that Catalyst models consist of are represented using so-called [*symbolic variables*](@ref ref) (and how these enable the forming of [*symbolic expressions*](@ref ref)). We have described how, during [programmatic modelling](@ref ref), the user have [direct access to these](@ref ref) and can use it to [their advantage](@ref ref). We have also described how, during [DSL-based modelling](@ref ref), the need for symbolic representation can be circumvented by [using `@unpack`](@ref ref) or [creating an observable](@ref ref). However, sometimes, it is easier to *directly access a symbolic variable through the model itself*.

Let us consider the following [two-state model](@ref ref)
```@example model_accessing_symbolic_variables
using Catalyst
rs = @reaction_network begin
    (k1,k2), X1 <--> X2
end
```
If we wish to access onf of symbolic variable stored in it (here `X1`, `X2`, `k1`, and `k2`), we simply write 
```@example model_accessing_symbolic_variables
rs.X1
```
to access e.g. `X1`. This symbolic variable can be used just like those [declared using `@parameters` and `@species`](@ref ref):
```@example model_accessing_symbolic_variables
u0 = [rs.X1 => 1.0, rs.X2 => 2.0]
ps = [rs.k1 => 2.0, rs.k2 => 4.0]
oprob = ODEProblem(rs, u0, (0.0, 10.0), ps)
sol = solve(oprob)
nothing # hide
```
We can also use them to form symbolic expressions:
```@example model_accessing_symbolic_variables
Xtot = rs.X1 + rs.X2
```
wish can be used when we e.g. [plot our simulation](@ref ref):
```@example model_accessing_symbolic_variables
plot(sol; idxs = [rs.X1, rs.X2, Xtot])
```

Next we create our two-state model programmatically:
```@example model_accessing_symbolic_variables
t = default_t()
@species X1(t) X2(t)
@parameters k1 k2
rxs = [
    Reaction(k1, [X1], [X2]),
    Reaction(k2, [X2], [X1])
]
@named rs_prog = ReactionSystem(rxs, t)
rs_prog = compelte(rs_prog)
nothing # hide
```
here, we can confirm that the symbolic variables we access through out model are identical to those we used to create it:
```@example model_accessing_symbolic_variables
isequal(rs.k1, k1)
```

!!! warn
    When accessing model symbolic variables through the model (using e.g. `rs.X1`), it is important to first ensure that the [*model is complete](@ref ref).

## [Accessing basic model properties](@id model_accessing_basics)

### [Accessing model parameter and species](@id model_accessing_basics_parameters_n_species)
Previously we showed how to access individual parameters or species of a `ReactionSystem` model. Next, the `parameters` and `species` functions allow us to retrieve *all* parameters/species as a vector:
```@example model_accessing_basics
using Catalyst #
sir = @reaction_network begin
    α, S + I --> 2I
    β, I --> R
end
parameters(sir)
```
```@example model_accessing_basics
species(sir)
```
These vectors contain the exact same symbolic variables that we would access through the system:
```@example model_accessing_basics
issetequal([sir.S, sir.I, sir.R], species(sir))
```

If we wish to count the number of parameters or species in a system, we can do this directly through the `numparams` and `numspecies` functions:
```@example model_accessing_basics
numparams(sir)
```
```@example model_accessing_basics
numspecies(sir)
```

### [Accessing model reactions](@id model_accessing_basics_reactions)
A vector containing all a model's [reactions] can be retrieved using the `reactions` function:
```@example model_accessing_basics
reactions(sir)
```
We can count the number of reactions in a model using the `numreactions` function:
```@example model_accessing_basics
numreactions(sir)
```
Finally, a vector with all the reactions' rates can be retrieved using `reactionrates`:
```@example model_accessing_basics
reactionrates(sir)
```

### [Accessing content of models coupled to equations](@id model_accessing_basics_reactions)

Previously, we have shown how to [couple equations to a chemical reaction network model](@ref ref), creating models containing [non-species unknowns (variables)](@ref ref). Here we create a birth-death model where some nutrient supply (modelled through the variable $N$) is depleted in the presence of $X$:
```@example model_accessing_basics
using Catalyst #
coupled_crn = @reaction_network begin
    @equation D(N) ~ - N * X
    (p/(1+N),d), 0 <--> X
end
```
Here, the `unknown` function return all unknowns (species and variables):
```@example model_accessing_basics
unknowns(coupled_crn)
```
Meanwhile, `species` returns the species only, while `nonspecies` returns the variables only:
```@example model_accessing_basics
species(coupled_crn)
```
```@example model_accessing_basics
nonspecies(coupled_crn)
```

Similarly, the `equations` function returns a vector with all reactions and equations of the model (ordered so that reactions occur first and equations there after):
```@example model_accessing_basics
equations(coupled_crn)
```
Meanwhile, `reactions` returns the reactions only, while `nonreactions` returns any algebraic or differential equations:
```@example model_accessing_basics
reactions(coupled_crn)
```
```@example model_accessing_basics
nonreactions(coupled_crn)
```

### [Accessing other model properties](@id model_accessing_basics_others)
There exists several other functions for accessing model properties. 

The `observed`, `continuous_events`, `discrete_events` functions can be used to access a model's [observables](@ref ref), [continuous events]('ref ref), and [discrete events](@ref ref), respectively.

The `ModelingToolkit.get_iv` function can be used to retrieve a [model's independent variable](@ref ref):
```@example model_accessing_basics
ModelingToolkit.get_iv(sir)
```

## [Accessing properties of hierarchical models](@id model_accessing_hierarchical)
Previously, we have described how [compositional modelling can be used to create hierarchical models](@ref ref). There are some special considerations when accessing content of hierarchical models, which will be described below.

