# [The Catalyst DSL - Introduction](@id dsl_description)
In the [introduction to Catalyst](@ref introduction_to_catalyst) we described how the `@reaction_network` macro can be used to create chemical reaction network models. This macro enables a so-called [domain-specific language](https://en.wikipedia.org/wiki/Domain-specific_language) (DSL) for creating chemical reaction network (CRN) models.  This tutorial will give a basic introduction on how to create Catalyst models using this macro (from now onwards called the "*Catalyst DSL*"). A [follow-up tutorial](@ref) will describe some of the DSLs more advanced features.

The Catalyst DSL generates a [`ReactionSystem`](@ref) (the [julia structure](https://docs.julialang.org/en/v1/manual/types/#Composite-Types) Catalyst uses to represent CRN models). These can be created through alternative methods (e.g. [programmatically](@ref programmatic_CRN_construction)). A descriptions of the various ways to create `ReactionSystems`s can be found [here](@ref ref). [Previous](@ref ref) and [following](@ref ref) tutorials describe how to simulate models once they have been created using the DSL. This tutorial wil solely focus on options for model creation.

Before we  begin the tutorial, we will first load the `Catalyst` package (which is required to run the code).
```@example dsl_1
using Catalyst
```

### Quick-start summary


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
Here, $P$'s production rate will decay as $A$ is slowly removed from the system. We can [check the ODE this model produce by using `Latexify`](@ref ref):
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

Above we used a simple example where the rate was the product of a species and a parameter. However, any valid Julia expression of parameters, species, and values can be used. E.g the following is a valid model:
```@example dsl_1
rn_14 = @reaction_network begin
  2.0 + X^2, 0 --> X + Y
  k1+k2^k3, X --> ∅
  pi*X/(sqrt(2)+Y), Y → ∅
end
```

### Time-dependant rates
Previously we have assumed that that the rates are independent on the [time variable, $t$](@ref ref). However, time-dependent reactions are also possible. Here, simply use `t` to represent the time variable. E.g., to create a production/degradation model where the production rate decays as time progress, we can use:
```@example dsl_1
rn_14 = @reaction_network begin
    kp/t, 0 --> P
    d, P --> 0
end
```

Like previously, `t` can be part of any valid expression. E.g. to create a reaction with a cyclic rate (e.g. to representing a [circadian system](https://en.wikipedia.org/wiki/Circadian_rhythm)) we can use:
```@example dsl_1
rn_15 = @reaction_network begin
    A*(sin(2π*f*t - ϕ)+1)/2, 0 --> P
    d, P --> 0
end
```

### Using functions in rates
It is possible for the rate to use Julia function. These can either be functions from Julia's standard library:  
```@example dsl_1
rn_16 = @reaction_network begin
    d, A --> 0
    kp*sqrt(A), 0 --> P
end
```
or ones defined by the user:
```@example dsl_1
custom_function(p1,p2,X) = (p1+E)/(p2+E)
rn_17 = @reaction_network begin
    d, A --> 0
    custom_function(k1,k2,E), 0 --> P
end
```

### Using pre-declared Michaelis-Menten or Hill function rate
Two function frequently used within systems biology are the [*Michaelis-Menten*](https://en.wikipedia.org/wiki/Michaelis%E2%80%93Menten_kinetics) and [*Hill*](https://en.wikipedia.org/wiki/Hill_equation_(biochemistry)) functions. These are pre-defined in Catalyst and can be called using `mm(X,v,K)` and `hill(X,v,K,n)`. E.g. a self-activation loop where $X$ activates its own production through a hill function can be created using:
```@example dsl_1
custom_function(p1,p2,X) = (p1+E)/(p2+E)
rn_18 = @reaction_network begin
    hill(X,v,K,n), 0 --> P
    d, X --> 0
end
```

Catalyst comes with the following predefined functions:
- The Michaelis-Menten function: $mm(X,v,K) = v*X/(X + K)$.
- The repressive Michaelis-Menten function: $mmr(X,v,K) = v*K/(X + K)$.
- The Hill function: $hill(X,v,K,n) = v*(X^n)/(X^n + K^n)$.
- The repressive Hill function: $hillr(X,v,K,n) = v*(K^n)/(X^n + K^n)$.
- The activating/repressive Hill function: $hillar(X,Y,v,K,n) = v*(X^n)/(X^n + Y^n + K^n)$.


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







