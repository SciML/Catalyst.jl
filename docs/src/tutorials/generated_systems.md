# The generated [`ReactionSystem`](@ref) and [`Reaction`](@ref)s

### [`ReactionSystem`](@ref) Accessors

The `@reaction_network` macro generates a [`ReactionSystem`](@ref), an instance
of a `ModelingToolkit.AbstractTimeDependentSystem`. It has a number of fields
that can be accessed using the [Catalyst.jl API](@ref) or the
[ModelingToolkit.jl Abstract System
Interface](https://mtk.sciml.ai/dev/basics/AbstractSystem/). Below we overview
these components.

There are three basic sets of convenience accessors that will return information
either from a top-level system, the top-level system and all sub-systems that
are also `ReactionSystem`s (i.e. the full reaction-network), or the top-level
system and all sub-systems (i.e. the full model). To retrieve info from just a
base [`ReactionSystem`](@ref) `rn`, ignoring sub-systems of `rn`, one can use the ModelingToolkit
accessors:

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

* [`states(rn)`](@ref) returns all species *and variables* across the system and
  *all sub-systems*.
* [`parameters(rn)`](@ref) returns all parameters across the system and *all
  sub-systems*.
* [`equations(rn)`](@ref) returns all [`Reaction`](@ref)s and all
  [`Equations`](@ref) defined across the system and *all sub-systems*.

These accessors will allocate unless there are no subsystems. In the latter case
they are equivalent to the corresponding `get_*` functions.

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

