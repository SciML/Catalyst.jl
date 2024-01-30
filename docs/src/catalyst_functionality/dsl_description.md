# [The Catalyst DSL - Introduction](@id dsl_description)
In the [introduction to Catalyst](@ref introduction_to_catalyst) we described how the `@reaction_network` macro can be used to create chemical reaction network models. This macro enables a so-called [domain-specific language](https://en.wikipedia.org/wiki/Domain-specific_language) (DSL) for creating chemical reaction network (CRN) models.  This tutorial will give a basic introduction on how to create Catalyst models using this macro (from now onwards called the "*Catalyst DSL*"). A [follow-up tutorial](@ref) will describe some of the DSLs more advanced features.

The Catalyst DSL generates a [`ReactionSystem`](@ref) (the [julia structure](https://docs.julialang.org/en/v1/manual/types/#Composite-Types) Catalyst uses to represent CRN models). These can be created through alternative methods (e.g. [programmatically](@ref programmatic_CRN_construction)). A descriptions of teh various ways to create `ReactionSystems`s can be found [here](@ref ref).

Before we  begin the tutorial, we will first load the `Catalyst` package (which is required to run the code).
```@example dsl_1
using Catalyst
```

### QUick-start summary


## [Basic syntax](@id dsl_description_basic_syntax)
The basic syntax of the DSL is
```@example dsl_1
rn = @reaction_network begin
  2.0, X --> Y
  1.0, Y --> X
end
```
Here, you start with `@reaction_network begin`, next list all of the model's reaction, and finish with an `end`. Each reaction consists of
- A rate
- A (potentially empty) set of substrates.
- A (potentially empty) set of products.

Each reaction line declares, in order, the rate, the substrate(s), and the products. The rate is separated from the substrate(s) by a `,`, ad the substrate(s) from the production by a `-->` (other arrows, however, are [also possible](@ref ref)). In the above examples our model consists of two reactions. In the first one, `X` (the single substrate) becomes `Y` (the single product) at rate `2.0`. 

Finally, `rn = ` is used to store the model in the variable `rn` (a normal Julia variable, which does not need to be called `rn`).

## [Defining parameters and species in the DSL](@id dsl_description_parameters_basics)
Typically, the rates are not constants, but rather parameters (which values can be set e.g. at [the beginning of each simulation](@ref ref)). To set parametric rates, simply use whichever symbol you which to represent your parameter with. E.g. to set the above rates to `a` and `b`, we use:
```@example dsl_1
rn1 = @reaction_network begin
  a, X --> Y
  b, Y --> X
end
```

Here we have used single-character symbols to designate all species and parameters. Multi-character symbols, however, are also permitted. E.g. we could call the rates `kX` and `kY`:
```@example dsl_1
rn1 = @reaction_network begin
  kX, X --> Y
  kY, Y --> X
end
nothing # hide
```
Generally, anything that is a [permitted Julia variable name](@id https://docs.julialang.org/en/v1/manual/variables/#man-allowed-variable-names) can be used to designate a species or parameter in Catalyst.

## [Different types of reactions](@id dsl_description_reactions)

### Reactions with multiple substrate(s) or product(s)
Previously, our reactions have had a single substrate and a single product. However, reactions with multiple substrates and/or products are possible. Here, all the substrates (or products) are listed and separated by a `+`. E.g. to create a model where `X` and `Y` binds (at rate `kB`) to form `XY` (which then can dissociate at, rate `kD`, to form `XY`) we use:
```@example dsl_1
rn2 = @reaction_network begin
  kB, X + Y --> XY
  kD, XY --> X + Y
end
```
Reactions can have any number of substrates and products, and their names does not need to have any relationship to each other, as demonstrated by the following mock-model:
```@example dsl_1
rn3 = @reaction_network begin
  k, X + Y + Z --> A + B + C + D
end
```

### [Reactions with degradation or production](@id dsl_description_reactions_degradation_and_production)
Some reactions has no products, in which case the substrate(s) are degraded (ie.e removed from the system). To denote this, set the reaction's right-hand side to `0`. Similarly, some reactions has no substrate(s), in which case the product(s) are produced (i.e. added to the system). This is denoted by setting the left-hand side to `0`. E.g. to create a model where a single species $X$ is both created and degraded, we use:
```@example dsl_1
rn4 = @reaction_network begin
  p, 0 --> X
  d, X --> 0
end
```

### Reactions with non-unitary stoichiometries
Reactions may include multiple copies of the same reactant (i.e. a substrate or a product). To specify this, the reactant is preceded by a number indicating the number of copies (also called the reactant's *stoichiometry*). E.g. to create a model where two copies of `X` dimerise to form `X2` (which may dissociate back to two `X` copies) we use:
```@example dsl_1
rn5 = @reaction_network begin
  kB, 2X --> X2
  kD, X2 --> 2X
end
```
Reactant which stoichiometry is not defined is assumed to have stoichiometry `1`. Any integer number can be used, furthermore, [decimal numbers and parameters can also be used as stoichiometries](@ref ref). A discussion of non-unitary (i.e. not equal to `1`) stoichiometries affects the created model can be found [here](@ref ref).

Stoichiometries can be combined with `( )` to define them for multiple reactants. Here, the following (mock) model declares the same reaction, both with and without this notation:
```@example dsl_1
rn6 = @reaction_network begin
  k, 2X + 3(Y + 2Z) --> 5(V + W)    
  k, 2X + 3Y + 6Z --> 5V + 5W    
end
nothing # hide
```

## [Bundling of similar reactions](@id dsl_description_reaction_bundling)

### Bi-directional (or reversible) reactions
As is the case for the following two-state model:
```@example dsl_1
rn7 = @reaction_network begin
  k1, X1 --> X2
  k2, X2 --> X1
end
```
it is common that reactions occurs in both direction. Here, it is possible ot bundle the reactions into a single line by using the `<-->` arrow. When we do this, we the rate term includes bother the rates (enclosed by a `()` and separated by a `,`). I.e. the two-state model can be declared using:
```@example dsl_1
rn7 = @reaction_network begin
  (k1,k2), X1 <--> X2
end
```
Here, the first rate (`k1`) denotes the *forward rate* and the second rate (`k2`) the *backwards rate*.

Catalyst also permits writing backwards reactions only. This uses the identical syntax to when forward reactions are created, but using the `<--` arrow:
```@example dsl_1
rn8 = @reaction_network begin
  k, X <-- Y
end
```
Here, the substrate(s) are on the right-hand side and the product(s) on the left-hand side. Hence, the above model can be written identically using:
```@example dsl_1
rn8 = @reaction_network begin
  k, Y --> X
end
```
Generally, using forward reactions is clearer than backwards ones, with the letter generally not being used.

### Bundling similar reactions on a singe line
There exists several other situation where models contain similar reactions (e.g. the systems where all system components degrade at the same rate). Reactions which share either rates, substrates, or products can be bundled into a single line. Here, the parts which are different for the reactions are written using `(,)` (containing one separate expression for each reaction). E.g., let us consider the following mode, where species $X$ and $Y$ both degrades at the rate $d$:
```@example dsl_1
rn8 = @reaction_network begin
  d, X --> 0
  d, Y --> 0
end
```
These share both the rates and products, however, the substrates are different. Hence, the reactions can be bundled into a single line using the common rate and product expression. However, we have to provide two separate substrate expressions:
```@example dsl_1
rn8 = @reaction_network begin
  d, (X,Y) --> 0
end
```
This declaration of the model is identical to the one previously.Reactions can share any subset of the rate, substrate, and product expression (the cases where they share all, or none, however, does not make sense to use). I.e. if the tw reactions had different degradation rates:
```@example dsl_1
rn9 = @reaction_network begin
  dX, X --> 0
  dY, Y --> 0
end
```
we could represent this using:
```@example dsl_1
rn9 = @reaction_network begin
  (dX,dY), (X,Y) --> 0
end
```

It is possible to use bundling for any number of reactions. E.g. in the following model we bundle the conversion of a species $X$ between its various forms (where all reactions uses the same rate $k$):
```@example dsl_1
rn10 = @reaction_network begin
  k, (X0,X1,X2,X3) --> (X1,X2,X3,X4)
end
```

It is possible to combine bundling with bi-directional reactions. In this case, the rate is first split into the forward and backwards rate(s). These may then (or may not) indicate several rates. We exemplify this using the two following (identical) networks, created with and without bundling.
```@example dsl_1
rn11 = @reaction_network begin
  kf, S --> P1
  kf, S --> P2
  kb_1, P1 --> S
  kb_2, P2 --> S
end
```
```@example dsl_1
rn11 = @reaction_network begin
  (kf, (kb_1, kb_2)), S <--> (P1,P2)
end
```

Like when we designated stoichiometries, reaction bundling can be applied very generally to create some truly complicated reactions:
```@example dsl_1
rn12 = @reaction_network begin
    ((pX, pY, pZ),d), (0, Y0, Z0) <--> (X, Y, Z1+Z2)
end
```
However, as for the above model, bundling reactions too zealously can reduce (rather improve) the model's readability.   

## [Non-constant reaction rates](@id dsl_description_nonconstant_rates)
So far we have assumed that all reaction rates are constant (being either a number of a parameter). Non-constant rates that depend on one (or several) species are also possible. More generally, the rate can be any valid expression of parameters and species.

Let us consider a model with an activator ($A$, which degraded at a constant rate) and a protein ($P$). The production rate of $P$ depends both on $A$ and parameter ($kp$). We model this through:
```@example dsl_1
rn_13 = @reaction_network begin
    d, A --> 0
    kp*A, 0 --> P
end
```
Here, $P$'s production rate will decay as $A$ is slowly removed from the system. We can [check the ODE this model produce by using `Latexify`]:
```@example dsl_1
using Latexify
latexify(rn_13; form=:ode)
```

In this case, we can generate an equivalent model by instead adding $A$ as both a substrate and a product to $A$'s production reaction:
```@example dsl_1
rn_13_alt = @reaction_network begin
    d, A --> 0
    kp, A --> A + P
end
```
We can confirm that this generate the same ODE:
```@example dsl_1
latexify(rn_13_alt; form=:ode)
```
Here, while the model will generate the same ODE, SDE, and jump simulations, the chemical reaction network models themselves are not equivalent. Generally, as pointed out in the two notes below, using the second form is preferable.
!!! warn
    While `rn_13` and `rn_13_alt` will generate equivalent simulations, for [jump simulations](@ref ref), the first model will have [reduced performance](@ref ref) (which generally are more performant when rates are constant).

!!! danger
    Catalyst automatically infers whether variables appearing in the DSL are species or parameters (as described [here](@ref ref)). Generally, anything that does not appear as a reactant is inferred to be a parameter. This mean that if you want to model a reaction activated by a species (e.g. `kp*A, 0 --> P`), but that species does not occur as a reactant, it will be interpreted as a parameter. This can be handled by [manually declaring the system species](@ref ref). A full example of how to do this for this example can be found [here](@ref ref).


### Time-dependant rates

### Using functions in rates

### Using pre-declared Michaelis-Menten or Hill function rate


## [Using special symbols](@id dsl_description_symbols)
Julia permits any unicode characters to be used in variable names, thus Catalyst are able to use these as well. Below we describe some case where these to alter the appearance of models. No functionality is, however, tied to this.

### [Using ∅ to denote degradation/production reactions](@id dsl_description_symbols_empty_set)
Previously, we described how `0` could be used to [denote degradation or production reactions](@ref dsl_description_reactions_degradation_and_production). Catalyst permits the user to instead use the `∅`. E.g. the production/degradation system alternatively can be written as:
```@example dsl_1
rn4 = @reaction_network begin
  p, ∅ --> X
  d, X --> ∅
end
```

### Using special arrow symbols
Catalyst uses `-->`, `<-->`, and `<--` to denote forward, reversible, and backwards reactions, respectively. Several unicode representation of these arrows are available. Here,
- `>`, `→`, `↣`, `↦`, `⇾`, `⟶`, `⟼`, `⥟`, `⥟`, `⇀`, and `⇁` can be used to represent forward reactions.
- `↔`, `⟷`, `⇄`, `⇆`, `⇌`, `⇋`, , and `⇔` can be used to represent reversible reactions.
- `<`, `←`, `↢`, `↤`, `⇽`, `⟵`, `⟻`, `⥚`, `⥞`, `↼`, , and `↽` can be used to represent backwards reactions. 

E.g. the production/degradation system alternatively can be written as:
```@example dsl_1
rn4 = @reaction_network begin
  p, ∅ → X
  d, X → ∅
end
```

### Using special symbols to denote species or parameters
A range of possible characters are available which can be incorporated into species and parameter names. This includes, but is not limited to:
- Greek letters (e.g `α`, `σ`, `τ`, and `Ω`).
- Superscript and subscript characters (to create e.g `k₁`, `k₂`, `Xₐ`, and `Xᴾ`).
- Non-latin, non-greek, letters (e.g. `ä`, `Д`, `س`, and `א`).
- Other symbols (e.g. `£`, `ℂ`, `▲`, and `♠`).

An example of how this can be used to create neat-looking model can be found in [Schwall et al. (2021)](https://www.embopress.org/doi/full/10.15252/msb.20209832) where it was used it the creation a model of the sigma factor V circuit in the bacteria *Bacillus subtilis*:
```@example dsl_1
σᵛ_model = @reaction_network begin
  v₀ + hill(σᵛ,v,K,n), ∅ → σᵛ + A
  kdeg, (σᵛ, A, Aσᵛ) → ∅
  (kB,kD), A + σᵛ ↔ Aσᵛ
  L, Aσᵛ → σᵛ
end
nothing # hide
```

This functionality can also be used to create less-serious models:
rn_13 = @reaction_network begin
    🍦, 😢 --> 😃
end

It should be noted that the following symbols are *not permitted* to be used as species or parameter names:
- `ℯ` (used in Julia to denote [Euler's constant](https://en.wikipedia.org/wiki/Euler%27s_constant))
- `pi` and `π` (used in Julia to denote [`3.1415926535897...`](https://en.wikipedia.org/wiki/Pi)).
- `t` (used to denote the [time variable](@ref ref)).
- `∅` ([used for production/degradation reactions](@ref dsl_description_symbols_empty_set))
- `im` (used in Julia to represent [complex numbers](https://docs.julialang.org/en/v1/manual/complex-and-rational-numbers/#Complex-Numbers)).
- `nothing` ([used in Julia](https://docs.julialang.org/en/v1/base/constants/#Core.nothing)).
- `Γ` (used by Catalyst to represent [conserved quantities](@ref ref)).






## [Basic syntax](@id basic_examples)

The `@reaction_network` macro allows the (symbolic) specification of reaction
networks with a simple format. Its input is a set of chemical reactions, and
from them it generates a symbolic [`ReactionSystem`](@ref) reaction network
object. The `ReactionSystem` can be used as input to ModelingToolkit
`ODEProblem`, `NonlinearProblem`, `SteadyStateProblem`, `SDEProblem`,
`JumpProblem`, and more. `ReactionSystem`s can also be incrementally extended as
needed, allowing for programmatic construction of networks and network
composition.

The basic syntax is:

```@example dsl_1
rn = @reaction_network begin
  2.0, X + Y --> XY
  1.0, XY --> Z1 + Z2
end
```
where each line of the [`@reaction_network`](@ref) macro corresponds to a
chemical reaction. Each reaction consists of a reaction rate (the expression on
the left-hand side of  `,`), a set of substrates (the expression in-between `,`
and `-->`), and a set of products (the expression on the right-hand side of
`-->`). The substrates and the products may contain one or more reactants,
separated by `+`. The naming convention for these is the same as for normal
variables in Julia.

The chemical reaction model is generated by the `@reaction_network` macro and
stored in the `rn` variable (a normal Julia variable, which does not need to be
called `rn`). It corresponds to a [`ReactionSystem`](@ref), a symbolic
representation of the chemical network.  The generated `ReactionSystem` can be
converted to a symbolic differential equation model via
```@example dsl_1
osys  = convert(ODESystem, rn)
```

We can then convert the symbolic ODE model into a compiled, optimized
representation for use in the SciML ODE solvers by constructing an `ODEProblem`.
Creating an `ODEProblem` also requires our specifying the initial conditions for
the model. We do this by creating a mapping from each symbolic variable
representing a chemical species to its initial value
```@example dsl_1
# define the symbolic variables
@variables t
@species X(t) Y(t) Z(t) XY(t) Z1(t) Z2(t)

# create the mapping
u0 = [X => 1.0, Y => 1.0, XY => 1.0, Z1 => 1.0, Z2 => 1.0]
```
Alternatively, we can create a mapping using Julia `Symbol`s for each variable,
and then convert them to a mapping involving symbolic variables like
```@example dsl_1
u0 = symmap_to_varmap(rn, [:X => 1.0, :Y => 1.0, :XY => 1.0, :Z1 => 1.0, :Z2 => 1.0])
```
Given the mapping, we can then create an `ODEProblem` from our symbolic `ODESystem`
```@example dsl_1
tspan = (0.0, 1.0)  # the time interval to solve on
oprob = ODEProblem(osys, u0, tspan, [])
```

Catalyst provides a shortcut to avoid having to explicitly `convert` to an
`ODESystem` and/or use `symmap_to_varmap`, allowing direct construction
of the `ODEProblem` like
```@example dsl_1
u0 = [:X => 1.0, :Y => 1.0, :XY => 1.0, :Z1 => 1.0, :Z2 => 1.0]
oprob = ODEProblem(rn, u0, tspan, [])
```

For more detailed examples, see the [Basic Chemical Reaction Network
Examples](@ref basic_CRN_examples).

## Defining parameters and species
Numeric parameter values do not need to be set when the model is created, i.e.
Catalyst supports symbolic parameters too:
```@example dsl_1
rn = @reaction_network begin
  k1, X --> Y
  k2, Y --> X
end
```
All symbols that do not appear as a substrate or product in a reaction are
designated by Catalyst as a parameter (i.e. all symbols appearing only within
rate expressions and/or as [stoichiometric coefficients](@ref parametric_stoichiometry)). In this example `X` and `Y`
appear as a substrates and products, but neither `k1` nor `k2`. Hence `k1` and `k2` are
designated as parameters. Later in this tutorial, we will describe how to manually specify what should be
considered a species or parameter.

## Production, Destruction, and Stoichiometry
Sometimes reactants are produced/destroyed from/to nothing. This can be
designated using either `0` or `∅`:
```@example dsl_1
rn = @reaction_network begin
  2.0, 0 --> X
  1.0, X --> 0
end
```
If several molecules of the same reactant are involved in a reaction, the
stoichiometry of a reactant in a reaction can be set using a number. Here, two
molecules of species `X` form the dimer `X2`:
```@example dsl_1
rn = @reaction_network begin
  1.0, 2X --> Y
end
```
this corresponds to the differential equation:
```@example dsl_1
convert(ODESystem, rn)
```
Other numbers than 2 can be used, and parenthesis can be used to reuse the same
stoichiometry for several reactants:
```@example dsl_1
rn = @reaction_network begin
  1.0, X + 2(Y + Z) --> W
end
```
Note, one can explicitly multiply by integer coefficients too, i.e.
```@example dsl_1
rn = @reaction_network begin
  1.0, X + 2*(Y + Z) --> W
end
```

## Arrow variants
A variety of Unicode arrows are accepted by the DSL in addition to `-->`. All of
these work:  `>`, `→` `↣`, `↦`, `⇾`, `⟶`, `⟼`, `⥟`, `⥟`, `⇀`, `⇁`. Backwards
arrows can also be used to write the reaction in the opposite direction. For
example, these reactions are equivalent:
```@example dsl_1
rn = @reaction_network begin
  1.0, X + Y --> XY
  1.0, X + Y → XY
  1.0, XY ← X + Y
  1.0, XY <-- X + Y
end
```

## Bi-directional arrows for reversible reactions
Bi-directional arrows, including bidirectional Unicode arrows like ↔, can be
used to designate a reversible reaction. For example, these two models are
equivalent:
```@example dsl_1
rn = @reaction_network begin
  2.0, X + Y --> XY
  2.0, X + Y <-- XY
end
```
```@example dsl_1
rn2 = @reaction_network begin
  (2.0,2.0), X + Y <--> XY
end
```

If the reaction rates in the backward and forward directions are different, they
can be designated in the following way:
```@example dsl_1
rn = @reaction_network begin
  (2.0,1.0), X + Y <--> XY
end
```
which is identical to
```@example dsl_1
rn = @reaction_network begin
  2.0, X + Y --> XY
  1.0, X + Y <-- XY
end
```

Finally, Catalyst also 

## Combining several reactions in one line
Several similar reactions can be combined in one line by providing a tuple of
reaction rates and/or substrates and/or products. If several tuples are provided,
they must all be of identical length. These pairs of reaction networks are all
identical.

Pair 1:
```@example dsl_1
rn1 = @reaction_network begin
  1.0, S --> (P1,P2)
end
```
```@example dsl_1
rn2 = @reaction_network begin
  1.0, S --> P1
  1.0, S --> P2
end
```
Pair 2:
```@example dsl_1
rn1 = @reaction_network begin
  (1.0,2.0), (S1,S2) --> P
end
```
```@example dsl_1
rn2 = @reaction_network begin
  1.0, S1 --> P
  2.0, S2 --> P
end
```
Pair 3:
```@example dsl_1
rn1 = @reaction_network begin
  (1.0,2.0,3.0), (S1,S2,S3) --> (P1,P2,P3)
end
```
```@example dsl_1
rn2 = @reaction_network begin
  1.0, S1 --> P1
  2.0, S2 --> P2
  3.0, S3 --> P3
end
```
This can also be combined with bi-directional arrows, in which case separate
tuples can be provided for the backward and forward reaction rates.
These reaction networks are identical
```@example dsl_1
rn1 = @reaction_network begin
 (1.0,(1.0,2.0)), S <--> (P1,P2)
end
```
```@example dsl_1
rn2 = @reaction_network begin
  1.0, S --> P1
  1.0, S --> P2
  1.0, P1 --> S
  2.0, P2 --> S
end
```

## Variable reaction rates
Reaction rates do not need to be a single parameter or a number, but can also be
expressions depending on time or the current amounts of system species (when, for
example, one species can activate the production of another). For instance, this
is a valid notation:
```@example dsl_1
rn = @reaction_network begin
  1.0, X --> ∅
  k*X, Y --> ∅
end
```
corresponding to the ODE model
```@example dsl_1
convert(ODESystem,rn)
```

With respect to the corresponding mass action ODE model, this is actually
equivalent to the reaction system
```@example dsl_1
rn = @reaction_network begin
  1.0, X --> ∅
  k, X + Y --> X
end
```
```@example dsl_1
convert(ODESystem,rn)
```
!!! note
    While the ODE models corresponding to the preceding two reaction systems are
    identical, in the latter example the `Reaction` stored in `rn` will be classified as
    [`ismassaction`](@ref) while in the former it will not, which can impact optimizations
    used in generating `JumpSystem`s. For this reason, it is recommended to use the
    latter representation when possible.

Most expressions and functions are valid reaction rates, e.g.:
```@example dsl_1
using SpecialFunctions
rn = @reaction_network begin
  2.0*X^2, 0 --> X + Y
  t*gamma(Y), X --> ∅
  pi*X/Y, Y --> ∅
end
```
where here `t` always denotes Catalyst's time variable. Please note that many
user-defined functions can be called directly, but others will require
registration with Symbolics.jl ([see the faq](@ref user_functions)).

## Explicit specification of network species and parameters
Recall that the `@reaction_network` macro automatically designates symbols used
in the macro as either parameters or species, with symbols that appear as a
substrate or product being species, and all other symbols becoming parameters
(i.e. those that only appear within a rate expression and/or as [stoichiometric coefficients](@ref parametric_stoichiometry)). Sometimes, one might want to manually override this default
behavior for a given symbol. E.g one might want something to be considered as a
species, even if it only appears within a rate expression. In the following
network
```@example dsl_1
rn = @reaction_network begin
  k*X, Y --> 0
end
```
`X` (as well as `k`) will be considered a parameter.

By using the `@species` and `@parameters` options within the `@reaction_network`
macro, one can manually declare that specified symbols should be considered a
species or parameter. E.g in:
```@example dsl_1
rn = @reaction_network begin
  @species X(t) Y(t)
  k*X, Y --> 0
end
```
`X` and `Y` are set as species. Please note that when declaring species using
the `@species` option, their dependant variable (almost always `t`) also needs
to be designated. Similarly in
```@example dsl_1
rn = @reaction_network begin
  @parameters k
  k*X, Y --> 0
end
```
both `X` and `k` will be considered as parameters. It is also possible to use
both options simultaneously, allowing users to fully specify which symbols are
species and/or parameters:
```@example dsl_1
rn = @reaction_network begin
  @species X(t) Y(t)
  @parameters k
  k*X, Y --> 0
end
```
Here, `X` and `Y` are designated as species and `k` as a parameter.

The lists provided to the `@species` and `@parameters` options do not need to be extensive. Any symbol that appears in neither list will use the default option as determined by the macro. E.g. in the previous example, where we only want to change the default designation of `X` (making it a species rather than a parameter), we can simply write:
```@example dsl_1
rn = @reaction_network begin
  @species X(t)
  k*X, Y --> 0
end
```

Finally, note that the `@species` and `@parameters` options can also be used in
`begin ... end` block form, allowing more formatted lists of species/parameters:
```@example dsl_1
rn = @reaction_network begin
  @parameters begin
      d1
      d2
  end
  @species begin
      X1(t)
      X2(t)
  end
  d2, X2 --> 0
  d1, X1 --> 0
end
```
This can be especially useful when declaring default values for clarity of model
specification (see the next section).

## [Setting default values for initial conditions and parameters](@id dsl_description_defaults)
When using the `@species` and ` @parameters` macros to declare species and/or
parameters, one can also provide default initial conditions for each species and
values for each parameter:
```@example dsl_1
rn = @reaction_network begin
  @species X(t)=1.0
  @parameters p=1.0 d=0.1
  p, 0 --> X
  d, X --> ∅
end
```
This system can now be simulated without providing initial condition or
parameter vectors to the DifferentialEquations.jl solvers:
```@example dsl_1
using DifferentialEquations, Plots
u0 = []
tspan = (0.0, 10.0)
p = []
oprob = ODEProblem(rn, u0, tspan, p)
sol = solve(oprob)
plot(sol)
```

When providing default values, it is possible to do so for only a subset of the
species or parameters, in which case the rest can be specified when constructing
the problem type to solve:
```@example dsl_1
rn = @reaction_network begin
  @species X(t)
  @parameters p=1.0 d
  p, 0 --> X
  d, X --> 0
end

u0 = [:X => 1.0]
tspan = (0.0, 10.0)
p = [:d => .1]
oprob = ODEProblem(rn, u0, tspan, p)
sol = solve(oprob)
plot(sol)
```

Finally, default values can be overridden by passing mapping vectors to the
DifferentialEquations.jl problem being constructed. Only those initial conditions
or parameters for which we want to change their value from the default will need to be passed
```@example dsl_1
u0 = [:X => 1.0]
tspan = (0.0, 10.0)
p = [:p => 2.0, :d => .1]   # we change p to 2.0
oprob = ODEProblem(rn, u0, tspan, p)
sol = solve(oprob)
plot(sol)
```

## [Setting initial conditions that depend on parameters](@id dsl_description_parametric_initial_conditions)
It is possible to set the initial condition of one (or several) species so that they depend on some system parameter. This is done in a similar way as default initial conditions, but giving the parameter instead of a value. When doing this, we also need to ensure that the initial condition parameter is a variable of the system:
```@example dsl_1
rn = @reaction_network begin
  @parameters X0
  @species X(t)=X0
  p, 0 --> X
  d, X --> ∅
end
```
We can now simulate the network without providing any initial conditions:
```@example dsl_1
u0 = []
tspan = (0.0, 10.0)
p = [:p => 2.0, :d => .1, :X0 => 1.0]
oprob = ODEProblem(rn, u0, tspan, p)
sol = solve(oprob)
plot(sol)
```

## Naming the generated `ReactionSystem`
ModelingToolkit uses system names to allow for compositional and hierarchical
models. To specify a name for the generated `ReactionSystem` via the
[`@reaction_network`](@ref) macro, just place the name before `begin`:
```@example dsl_1
rn = @reaction_network production_degradation begin
  p, ∅ --> X
  d, X --> ∅
end
ModelingToolkit.nameof(rn) == :production_degradation
```

## Pre-defined functions
Hill functions and a Michaelis-Menten function are pre-defined and can be used
as rate laws. Below, the pair of reactions within `rn1` are equivalent, as are
the pair of reactions within `rn2`:
```@example dsl_1
rn1 = @reaction_network begin
  hill(X,v,K,n), ∅ --> X
  v*X^n/(X^n+K^n), ∅ --> X
end
```
```@example dsl_1
rn2 = @reaction_network begin
  mm(X,v,K), ∅ --> X
  v*X/(X+K), ∅ --> X
end
```
Repressor Hill (`hillr`) and Michaelis-Menten (`mmr`) functions are also
provided:
```@example dsl_1
rn1 = @reaction_network begin
  hillr(X,v,K,n), ∅ --> X
  v*K^n/(X^n+K^n), ∅ --> X
end
```
```@example dsl_1
rn2 = @reaction_network begin
  mmr(X,v,K), ∅ --> X
  v*K/(X+K), ∅ --> X
end
```

Please see the API [Rate Laws](@ref api_rate_laws) section for more details.

## Including non-species variables
Non-species state variables can be specified in the DSL using the `@variables`
macro. These are declared similarly to species. For example,
```@example dsl_1
rn_with_volume = @reaction_network begin
  @variables V(t)
  k*V, 0 --> A
end
```
creates a network with one species
```@example dsl_1
species(rn_with_volume)
```
and one non-species
```@example dsl_1
nonspecies(rn_with_volume)
```
giving two state variables, always internally ordered by species and then
nonspecies:
```@example dsl_1
states(rn_with_volume)
```

`rn_with_volume` could then be extended with constraint equations for how `V(t)`
evolves in time, see the [associated tutorial](@ref constraint_equations).

## Specifying alternative time variables and/or extra independent variables
While the DSL defaults to allowing `t` as the time variable, one can use the
`@ivs` macro to specify an alternative independent variable. For example, to
make `s` the default time variable one can say
```@example dsl_1
rn_with_s = @reaction_network begin
    @ivs s
    @variables V(s)
    @species B(s)
    k, A + V*B --> C
end
show(stdout, MIME"text/plain"(), rn_with_s)  # hide
```
where we see all states are now functions of `s`.

Similarly, if one wants states to be functions of more than one independent
variable, for example to encode a spatial problem, one can list more than one
variable, i.e. `@ivs t x y`. Here the first listed independent variable is
always chosen to represent time. For example,
```@example dsl_1
rn_with_many_ivs = @reaction_network begin
    @ivs s x
    @variables V1(s) V2(s,x)
    @species A(s) B(s,x)
    k, V1*A --> V2*B + C
end
show(stdout, MIME"text/plain"(), rn_with_many_ivs)  # hide
```
Here again `s` will be the time variable, and any inferred species, `C` in this
case, are made functions of both variables, i.e. `C(s, x)`.

## [Interpolation of Julia variables](@id dsl_description_interpolation_of_variables)
The DSL allows Julia variables to be interpolated for the network name, within
rate constant expressions, or for species/stoichiometry within reactions. Using
the lower-level symbolic interface we can then define symbolic variables and
parameters outside of the macro, which can then be used within expressions in
the DSL (see the [Programmatic Construction of Symbolic Reaction Systems](@ref programmatic_CRN_construction)
tutorial for details on the lower-level symbolic interface). For example,
```@example dsl_1
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
```@example dsl_1
parameters(rn)
```
as are the species coming from interpolated variables
```@example dsl_1
species(rn)
```

!!! note
    When using interpolation, expressions like `2$spec` won't work; the
    multiplication symbol must be explicitly included like `2*$spec`.
