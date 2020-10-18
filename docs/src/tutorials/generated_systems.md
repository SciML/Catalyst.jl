# The generated [`ReactionSystem`](@ref) and [`Reaction`](@ref)s
The `@reaction_network` macro generates a [`ReactionSystem`](@ref) object, which
has a number of fields that can be accessed directly or via the [Catalyst.jl
API](@ref) (the recommended route). Below we list these components, with the recommended
API method listed first:

* [`species(rn)`](@ref), `states(rn)` and `rn.states` is a vector of all the chemical
  species within the system, each represented as a `ModelingToolkit.Variable`.
* [`params(rn)`](@ref), `parameters(rn)` and `rn.ps` is a vector of all the parameters
  within the system, each represented as a `ModelingToolkit.Variable`.
* [`reactions(rn)`](@ref), `equations(rn)` and `rn.eqs` is a vector of all the
  `Reaction`s within the system.
* `independent_variable(rn)` and `rn.iv` are the independent variable of the
  system, usually `t` for time, represented as a `ModelingToolkit.Variable`.

Each `Reaction` within `reactions(rn)` has a number of subfields. For `rx` a
`Reaction` we have:
* `rx.substrates`, a vector of `ModelingToolkit.Operation`s storing each
  substrate variable.
* `rx.products`, a vector of `ModelingToolkit.Operation`s storing each product
  variable.
* `rx.substoich`, a vector storing the corresponding integer stoichiometry of
  each substrate species in `rx.substrates`.
* `rx.prodstoich`, a vector storing the corresponding integer stoichiometry of
  each product species in `rx.products`.
* `rx.rate`, a `ModelingToolkit.Operation` representing the reaction rate. E.g.,
  for a reaction like `k*X, Y --> X+Y`, we'd have `rate = k*X`.
* `rx.netstoich`, a vector of pairs mapping the `ModelingToolkit.Variable` for
  each species that changes numbers by the reaction to how much it changes. E.g.,
  for `k, X + 2Y --> X + W`, we'd have `rx.netstoich = [Y => -2, W => 1]`.
* `rx.only_use_rate`, a boolean that is `true` if the reaction was made with
  non-filled arrows and should ignore mass action kinetics. `false` by default.

Empty `ReactionSystem`s can be generated via [`make_empty_network`](@ref) or
[`@reaction_network`](@ref) with no arguments. `ReactionSystem`s can be
programmatically extended using [`addspecies!`](@ref), [`addparam!`](@ref),
[`addreaction!`](@ref), [`@add_reactions`](@ref), or composed using `merge` and
`merge!`.
