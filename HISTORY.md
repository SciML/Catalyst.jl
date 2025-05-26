# Breaking updates and feature summaries across releases

## Catalyst unreleased (master branch)

## Catalyst 15.0
- The Catalyst release process is changing; certain core dependencies of
  Catalyst will now be capped to ensure Catalyst releases are only installed
  with versions of dependencies for which Catalyst CI and doc build tests pass
  (at the time the release is made). If you need a dependency version increased,
  please open an issue and we can update it and make a new Catalyst release once
  testing against the newer dependency version is complete.
- **BREAKING:** New formula for inferring variables from equations (declared
  using the `@equations` options) in the DSL. The order of inference of
  species/variables/parameters is now:
    1. Every symbol explicitly declared using `@species`, `@variables`, and
       `@parameters` are assigned to the correct category.
    2. Every symbol used as a reaction reactant is inferred as a species.
    3. Every symbol not declared in (1) or (2) that occurs in an expression
       provided after `@equations` is inferred as a variable.
    4. Every symbol not declared in (1), (2), or (3) that occurs either as a
       reaction rate or stoichiometric coefficient is inferred to be a
       parameter. E.g. in
       ```julia
       @reaction_network begin
           @equations V1 + S ~ V2^2
           (p + S + V1), S --> 0
       end
       ```
       `S` is inferred as a species, `V1` and `V2` as variables, and `p` as a
       parameter. The previous special cases for the `@observables`, `@compounds`,
       and `@differentials` options still hold. Finally, the `@require_declaration`
       options (described in more detail below) can now be used to require everything
       to be explicitly declared.
- **BREAKING:** New formula for determining whether the default differentials
  have been used within an `@equations` option. Now, if any expression `D(...)`
  is encountered (where `...` can be anything), this is inferred as usage of the
  default differential D. E.g. in the following equations `D` is inferred as a
  differential with respect to the default independent variable:
  ```julia
  @reaction_network begin
    @equations D(V) + V ~ 1
  end
  @reaction_network begin
    @equations D(D(V)) ~ 1
  end
  ```
  Please note that this cannot be used at the same time as `D` is used to
  represent a species, variable, or parameter (including if these are implicitly
  designated as such by e.g. appearing as a reaction reactant). 
- **BREAKING:** Array symbolics support is more consistent with ModelingToolkit
  v9. Parameter arrays are no longer scalarized by Catalyst, while species and
  variables arrays still are (as in ModelingToolkit). As such, parameter arrays
  should now be specified as arrays in value mappings, i.e.
  ```julia
  @parameters k[1:4]
  pmap = [k => rand(4)]
  ```
  While one can still manually scalarize a parameter array, it is recommended
  *not* to do this as it has signifcant performance costs with ModelingToolkit
  v9. Note, scalarized parameter arrays passed to the two-argument
  `ReactionSystem` constructor may become unscalarized.
- **BREAKING:** We have introduced a restriction on bundling of reactions in the
  DSL. Now, bundling is not permitted if multiple rates are provided but only
  one set each of substrates/products. E.g. this model:
  ```julia
  @reaction_network begin
    (k1,k2), X --> Y
  end
  ```
  will now throw an error. The reason that users attempting to write
  bi-directional reactions but typing `-->` instead of `<-->` would get a wrong
  model. We decided that this kind of bundling was unlikely to be used, and
  throwing errors for people who made the typo was more important. If you use
  this type of bundling and it indeed is useful to you, please raise and issue
  and we will see if we can sort something out.
- **BREAKING:** Catalyst's network visualization capability has shifted from
  using Graphviz to [GraphMakie.jl](https://graph.makie.org/stable/). To use
  this functionality, load the GraphMakie extension by installing `Catalyst` and
  `GraphMakie`, along with a Makie backend like `CairoMakie` or `GLMakie`. There
  are two new methods for visualizing graphs: `plot_network` and
  `plot_complexes`, which respectively display the species-reaction graph and
  complex graph.
  ```julia
  using Catalyst, GraphMakie, GLMakie
  brusselator = @reaction_network begin
     A, ∅ --> X
     1, 2X + Y --> 3X
     B, X --> Y
     1, X --> ∅
  end
  plot_network(brusselator)
  ```
- **BREAKING:** The letter Ø (used in Danish/Norwegian alphabet) is now
  considered the same as ∅ (empty set). It can no longer be used as a
  species/parameter.
- **BREAKING:** When converting a Catalyst `ReactionSystem` to a ModelingToolkit
  system, for example an `ODESystem`, Catalyst defined functions like
  `hill(A,B,C,D)` are now replaced with the explicit rational function they
  represent in the equations of the generated system. For example `mm(X,v,K)`
  will be replaced with `v*X / (X + K)`. This can be disabled by passing the
  keyword argument `expand_catalyst_funs = false`. e.g.
  ```julia
  using Catalyst
  rn = @reaction_network begin
    hill(X,v,K,n), A --> 0
  end
  osys = convert(ODESystem, rn)
  ```
  generates an ODE system with `D(A) ~ -((v*A(t)*(X^n)) / (K^n + X^n))`, while
  ```julia
  osys = convert(ODESystem, rn; expand_catalyst_funs = false)
  ```
  generates an ODE system with `D(A) ~ -A(t)*hill(X, v, K, n)`. This keyword
  argument can also be passed to problems defined over `ReactionSystem`s, i.e.
  when calling `ODEProblem(rn, u0, tspan, p; expand_catalyst_funs = false)`.
- It is no longer recommended to install and use the full OrdinaryDiffEq library
  to access specific ODE solvers. Instead, only install the specific
  OrdinaryDiffEq sub-libraries that contain the desired solver. This
  significantly reduces installation and package loading times. I.e. to use the
  default solver that auto-switches between explicit and implicit methods,
  install `OrdinaryDiffEqDefault`. To use `Tsit5` install `OrdinaryDiffEqTsit5`,
  etc. The possible sub-libraries, each containing different solvers, can be
  viewed [here](https://github.com/SciML/OrdinaryDiffEq.jl/tree/master/lib).
- It should now be safe to use `remake` on problems which have had conservation
  laws removed with the exception of `NonlinearProblem`s or `NonlinearSystem`s.
  For NonlinearProblems it is safe to use `remake` if only updating `u0` values,
  but it is not safe to update the value of the conserved constant, `Γ`. See
  [the FAQ](https://docs.sciml.ai/Catalyst/stable/faqs/#faq_remake_nonlinprob)
  for details.
- Functional (e.g. time-dependent) parameters can now be used in Catalyst
  models. These can e.g. be used to incorporate arbitrary time-dependent
  functions (as a parameter) in a model. For more details on how to use these,
  please read:
  https://docs.sciml.ai/Catalyst/stable/model_creation/functional_parameters/.
- Scoped species/variables/parameters are now treated similar to the latest MTK
  releases (≥ 9.49).
- A tutorial on making interactive plot displays using Makie has been added.
- The BifurcationKit extension has been updated to v.4.
- There is a new DSL option `@require_declaration` that will turn off automatic
  inferring for species, parameters, and variables in the DSL. For example, the
  following will now error:
  ```julia
  rn = @reaction_network begin
    @require_declaration
    (k1, k2), A <--> B
  end
  ```
  When this flag is set, all symbolics must be explicitly declared. 
  ```julia
  rn = @reaction_network begin
    @species A(t) B(t)
    @parameters k1 k2
    @require_declaration
    (k1, k2), A <--> B
  end
  ```

## Catalyst 14.4.1
- Support for user-defined functions on the RHS when providing coupled equations 
  for CRNs using the @equations macro. For example, the following now works: 
  ```julia
  using Catalyst
  f(A, t) = 2*A*t
  rn = @reaction_network begin
    @equations D(A) ~ f(A,t)
  end
  ```
  Note that user-defined functions will not work on the LHS of equations. 

## Catalyst 14.4
- Symbolics 6 support.



## Catalyst 14.3
- Support for simulating stochastic chemical kinetics models with explicitly
  time-dependent propensities (i.e. where the resulting `JumpSystem` contains
  `VariableRateJump`s). As such `JumpProblem`s need to be defined over
  `ODEProblem`s or `SDEProblem`s instead of `DiscreteProblem`s we have
  introduced a new input struct, `JumpInputs`, that handles selecting via
  analysis of the generated `JumpSystem`, i.e. one can now say
  ```julia
  using Catalyst, OrdinaryDiffEq, JumpProcesses, Plots
  rn = @reaction_network begin
      k*(1 + sin(t)), 0 --> A
  end
  jinput = JumpInputs(rn, [:A => 0], (0.0, 10.0), [:k => .5])
  # note that jinput.prob isa ODEProblem
  jprob = JumpProblem(jinput)
  sol = solve(jprob, Tsit5())
  plot(sol, idxs = :A)

  rn = @reaction_network begin
      k, 0 --> A
  end
  jinput = JumpInputs(rn, [:A => 0], (0.0, 10.0), [:k => .5])
  # note that jinput.prob isa DiscreteProblem
  jprob = JumpProblem(jinput)
  sol = solve(jprob)
  plot(sol, idxs = :A)
  ```
  When calling solve for problems with explicit time-dependent propensities,
  i.e. where `jinput.prob isa ODEProblem`, note that one must currently
  explicitly select an ODE solver to handle time-stepping and integrating the
  time-dependent propensities.
- Note that solutions to jump problems with explicit time-dependent
  propensities, i.e. a `JumpProblem` over an `ODEProblem`, require manual
  selection of the variables to plot. That is, currently `plot(sol)` will error
  in this case due to limitations in the SciMLBase plot recipe.

## Catalyst 14.2
- Support for auto-algorithm selection in `JumpProblem`s. For systems with only
  propensities that do not have an explicit time-dependence (i.e. that are not
  `VariableRateJump`s in JumpProcesses), one can now run model simulations via
  ```julia
  using Catalyst, JumpProcesses
  model = @reaction_network begin
    kB, S + E --> SE
    kD, SE --> S + E
    kP, SE --> P + E
  end
  u0 = [:S => 50, :E => 10, :SE => 0, :P => 0]
  tspan = (0., 200.)
  ps = [:kB => 0.01, :kD => 0.1, :kP => 0.1]
  dprob = DiscreteProblem(model, u0, tspan, ps)
  jprob = JumpProblem(model, dprob)
  sol = solve(jprob)
  ```
  For small systems this will just use Gillespie's `Direct` method, transitioning to using `RSSA` and `RSSACR` as system size increase. Once can still manually select a given SSA, but no longer needs to specify `SSAStepper` when calling `solve`, i.e.
  ```julia
  # use the SortingDirect method instead
  jprob = JumpProblem(model, dprob, SortingDirect())
  sol = solve(jprob)
  ```
- Latexify recipe improvements including display fixes for array symbolics.
- Deficiency one and concentration robustness checks.

## Catalyst 14.1.1
The expansion of `ReactionSystem` models to spatial lattices has been enabled. Here follows a
simple example where a Brusselator model is expanded to a 20x20 grid of compartments, with diffusion
for species X, and then simulated using ODEs. Finally, an animation of the simulation is created.
```julia
using Catalyst, CairoMakie, OrdinaryDiffEq

# Create `LatticeReactionSystem` model.
brusselator = @reaction_network begin
    A, ∅ --> X
    1, 2X + Y --> 3X
    B, X --> Y
    1, X --> ∅
end
diffusion_rx = @transport_reaction D X
lattice = CartesianGrid((20,20))
lrs = LatticeReactionSystem(brusselator, [diffusion_rx], lattice)

# Create a spatial `ODEProblem`.
u0 = [:X => rand(20, 20), :Y => 10.0]
tspan = (0.0, 40.0)
ps = [:A => 1.0, :B => 4.0, :D => 0.2]
oprob = ODEProblem(lrs, u0, tspan, ps)

# Simulate the ODE and plot the results.
sol = solve(oprob, FBDF())
lattice_animation(sol, :X, lrs, "brusselator.mp4")
```
The addition of spatial modelling in Catalyst contains a large number of new features, all of which are
described in the [corresponding documentation](https://docs.sciml.ai/Catalyst/stable/spatial_modelling/lattice_reaction_systems/).

## Catalyst 14.0.1
Bug fix to address that independent variables, like time, should now be `@parameters`
according to MTKv9. Converted internal time variables to consistently use `default_t()`
to hopefully avoid such issues going forward.

## Catalyst 14.0

#### Breaking changes
Catalyst v14 was prompted by the (breaking) release of ModelingToolkit v9, which
introduced several breaking changes to Catalyst. A summary of these (and how to
handle them) can be found
[here](https://docs.sciml.ai/Catalyst/stable/v14_migration_guide/). These are
briefly summarised in the following bullet points:
- `ReactionSystem`s must now be marked *complete* before they are exposed to
  most forms of simulation and analysis. With the exception of `ReactionSystem`s
  created through the `@reaction_network` macro, all `ReactionSystem`s are *not*
  marked complete upon construction. The `complete` function can be used to mark
  `ReactionSystem`s as complete. To construct a `ReactionSystem` that is not
  marked complete via the DSL the new `@network_component` macro can be used.
- The `states` function has been replaced with `unknowns`. The `get_states`
  function has been replaced with `get_unknowns`.
- Support for most units (with the exception of `s`, `m`, `kg`, `A`, `K`, `mol`,
  and `cd`) has currently been dropped by ModelingToolkit, and hence they are
  unavailable via Catalyst too. Its is expected that eventually support for
  relevant chemical units such as molar will return to ModelingToolkit (and
  should then immediately work in Catalyst too).
- Problem parameter values are now accessed through `prob.ps[p]` (rather than
  `prob[p]`).
-  ModelingToolkit currently does not support the safe application of the
   `remake` function, or safe direct mutation, for problems for which
   `remove_conserved = true` was used when updating the values of initial
   conditions. Instead, the values of each conserved constant must be directly
   specified.
- The `reactionparams`, `numreactionparams`, and `reactionparamsmap` functions
  have been deprecated and removed.
- To be more consistent with ModelingToolkit's immutability requirement for
  systems, we have removed API functions that mutate `ReactionSystem`s such as
  `addparam!`, `addreaction!`, `addspecies`, `@add_reactions`, and `merge!`.
  Please use `ModelingToolkit.extend` and `ModelingToolkit.compose` to generate
  new merged and/or composed `ReactionSystem`s from multiple component systems.

#### General changes
- `default_t()` and `default_time_deriv()` functions should be used for creating
  the default time independent variable and its differential. i.e.
  ```julia
  # do
  t = default_t()
  @species A(t)

  # avoid
  @variables t
  @species A(t)
- It is now possible to add metadata to individual reactions, e.g. using:
  ```julia
  rn = @reaction_network begin
      @parameters η
      k, 2X --> X2, [description="Dimerisation"]
  end
  getdescription(rn)
  ```
  a more detailed description can be found [here](https://docs.sciml.ai/Catalyst/dev/model_creation/dsl_advanced/#dsl_advanced_options_reaction_metadata).
- `SDEProblem` no longer takes the `noise_scaling` argument. Noise scaling is
  now handled through the `noise_scaling` metadata (described in more detail
  [here](https://docs.sciml.ai/Catalyst/stable/model_simulation/simulation_introduction/#simulation_intro_SDEs_noise_saling))
- Fields of the internal `Reaction` structure have been changed.
  `ReactionSystems`s saved using `serialize` on previous Catalyst versions
  cannot be loaded using this (or later) versions.
- A new function, `save_reactionsystem`, which permits the writing of
  `ReactionSystem` models to files, has been created. A thorough description of
  this function can be found
  [here](https://docs.sciml.ai/Catalyst/stable/model_creation/model_file_loading_and_export/#Saving-Catalyst-models-to,-and-loading-them-from,-Julia-files)
- Updated how compounds are created. E.g. use
  ```julia
  @variables t C(t) O(t)
  @compound CO2 ~ C + 2O
  ```
  to create a compound species `CO2` that consists of `C` and two `O`.
- Added documentation for chemistry-related functionality (compound creation and
  reaction balancing).
- Added function `isautonomous` to check if a `ReactionSystem` is autonomous.
- Added function `steady_state_stability` to compute stability for steady
  states. Example:
  ```julia
  # Creates model.
  rn = @reaction_network begin
      (p,d), 0 <--> X
  end
  p = [:p => 1.0, :d => 0.5]

  # Finds (the trivial) steady state, and computes stability.
  steady_state = [2.0]
  steady_state_stability(steady_state, rn, p)
  ```
  Here, `steady_state_stability` takes an optional keyword argument `tol =
  10*sqrt(eps())`, which is used to check that the real part of all eigenvalues
  are at least `tol` away from zero. Eigenvalues within `tol` of zero indicate
  that stability may not be reliably calculated.
- Added a DSL option, `@combinatoric_ratelaws`, which can be used to toggle
  whether to use combinatorial rate laws within the DSL (this feature was
  already supported for programmatic modelling). Example:
  ```julia
  # Creates model.
  rn = @reaction_network begin
      @combinatoric_ratelaws false
      (kB,kD), 2X <--> X2
  end
  ```
- Added a DSL option, `@observables` for [creating
  observables](https://docs.sciml.ai/Catalyst/stable/model_creation/dsl_advanced/#dsl_advanced_options_observables)
  (this feature was already supported for programmatic modelling).
- Added DSL options `@continuous_events` and `@discrete_events` to add events to
  a model as part of its creation (this feature was already supported for
  programmatic modelling). Example:
  ```julia
  rn = @reaction_network begin
      @continuous_events begin
          [X ~ 1.0] => [X ~ X + 1.0]
      end
      d, X --> 0
  end
  ```
- Added DSL option `@equations` to add (algebraic or differential) equations to
  a model as part of its creation (this feature was already supported for
  programmatic modelling). Example:
  ```julia
  rn = @reaction_network begin
      @equations begin
          D(V) ~ 1 - V
      end
      (p/V,d/V), 0 <--> X
  end
  ```
  couples the ODE $dV/dt = 1 - V$ to the reaction system.
- Coupled reaction networks and differential equation (or algebraic differential
  equation) systems can now be converted to `SDESystem`s and `NonlinearSystem`s.

#### Structural identifiability extension
- Added CatalystStructuralIdentifiabilityExtension, which permits
  StructuralIdentifiability.jl to be applied directly to Catalyst systems. E.g.
  use
  ```julia
  using Catalyst, StructuralIdentifiability
  goodwind_oscillator = @reaction_network begin
      (mmr(P,pₘ,1), dₘ), 0 <--> M
      (pₑ*M,dₑ), 0 <--> E
      (pₚ*E,dₚ), 0 <--> P
  end
  assess_identifiability(goodwind_oscillator; measured_quantities=[:M])
  ```
  to assess (global) structural identifiability for all parameters and variables
  of the `goodwind_oscillator` model (under the presumption that we can measure
  `M` only).
- Automatically handles conservation laws for structural identifiability
  problems (eliminates these internally to speed up computations).
- A more detailed of how this extension works can be found
  [here](https://docs.sciml.ai/Catalyst/stable/inverse_problems/structural_identifiability/).

#### Bifurcation analysis extension
- Add a CatalystBifurcationKitExtension, permitting BifurcationKit's
  `BifurcationProblem`s to be created from Catalyst reaction networks. Example
  usage:
  ```julia
  using Catalyst
  wilhelm_2009_model = @reaction_network begin
      k1, Y --> 2X
      k2, 2X --> X + Y
      k3, X + Y --> Y
      k4, X --> 0
      k5, 0 --> X
  end

  using BifurcationKit
  bif_par = :k1
  u_guess = [:X => 5.0, :Y => 2.0]
  p_start = [:k1 => 4.0, :k2 => 1.0, :k3 => 1.0, :k4 => 1.5, :k5 => 1.25]
  plot_var = :X
  bprob = BifurcationProblem(wilhelm_2009_model, u_guess, p_start, bif_par; plot_var = plot_var)

  p_span = (2.0, 20.0)
  opts_br = ContinuationPar(p_min = p_span[1], p_max = p_span[2], max_steps = 1000)

  bif_dia = bifurcationdiagram(bprob, PALC(), 2, (args...) -> opts_br; bothside = true)

  using Plots
  plot(bif_dia; xguide = "k1", guide = "X")
  ```
- Automatically handles elimination of conservation laws for computing
  bifurcation diagrams.
- Updated Bifurcation documentation with respect to this new feature.

## Catalyst 13.5
- Added a CatalystHomotopyContinuationExtension extension, which exports the `hc_steady_state` function if HomotopyContinuation is exported. `hc_steady_state` finds the steady states of a reaction system using the homotopy continuation method. This feature is only available for julia versions 1.9+. Example:
```julia
wilhelm_2009_model = @reaction_network begin
    k1, Y --> 2X
    k2, 2X --> X + Y
    k3, X + Y --> Y
    k4, X --> 0
end
ps = [:k1 => 8.0, :k2 => 2.0, :k3 => 1.0, :k4 => 1.5]
hc_steady_states(wilhelm_2009_model, ps)
```

## Catalyst 13.4
- Added the ability to create species that represent chemical compounds and know
  their constituents. For example, water can be created and queried as
  ```julia
  @variables t
  @species H(t) O(t)
  @compound H2O(t) 2*H O
  iscompound(H2O) == true
  isspecies(H2O) == true   # compounds are also species, so can be used in reactions
  isequal(components(H2O), [H, O])
  coefficients(H2O) == [2, 1]
  ```
- Added reaction balancing via the `balance_reaction` command, which returns a
  vector of balanced reaction versions, i.e.
  ```julia
  @variables t
  @species H(t) O(t) C(t)
  @compound CH4(t) C 4H
  @compound O2(t) 2O
  @compound CO2(t) C 2O
  @compound H2O(t) 2H O

  # unbalanced reaction to balance
  rx = Reaction(1.0, [CH4, O2], [CO2, H2O])

  # calculate a balanced version, this returns a vector
  # storing a single balanced version of the reaction in this case
  brxs = balance_reaction(rx)

  # what one would calculate by hand
  balanced_rx = Reaction(1.0, [CH4, O2], [CO2, H2O], [1, 2], [1, 2])

  # testing equality
  @test isequal(balanced_rx, first(brxs))
  ```
- Note that balancing works via calculating the nullspace of an associated
  integer matrix that stores in entry `(i,j)` a signed integer representing the
  number of times the `i`'th atom appears within the `j`th compound. The entry
  is positive for a substrate and negative for a product. One cannot balance a
  reaction involving compounds of compounds currently. A non-empty solution
  vector is returned if the reaction can be balanced in exactly one way with
  minimal coefficients while preserving the set of substrates and products, i.e.
  if the dimension of the nullspace is one. If the dimension is greater than one
  we return a `Reaction` for each nullspace basis vector, but note that they may
  currently interchange substrates and products (i.e. we do not solve for if
  there is a linear combination of them that preserves the set of substrates and
  products). An empty `Reaction` vector indicates it is not possible to balance
  the reaction.

## Catalyst 13.2
- Array parameters, species, and variables can be use in the DSL if explicitly
  declared with `@parameters`, `@species`, or `@variables` respectively, i.e.
  ```julia
  rn = @reaction_network begin
      @parameters k[1:2] a
      @variables (V(t))[1:2] W(t)
      @species (X(t))[1:2] Y(t)
      k[1]*a+k[2], X[1] + V[1]*X[2] --> V[2]*W*Y + B*C
  end
  ```

## Catalyst 13.1
- Non-species states can be declared in the DSL using `@variables`, and custom
  independent variables (instead of just `t`) using `@ivs`. For the latter, the
  first independent variable is always interpreted as the time variable, and all
  *discovered* species are created to be functions of all the `ivs`. For example
  in
  ```julia
  rn = @reaction_network begin
      @ivs s x
      @variables A(s) B(x) C(s,x)
      @species D(s) E(x) F(s,x)
      k*C, A*D + B*E --> F + H
  end
  ```
  `s` will be the time variable, `H = H(s,x)` will be made a function of `s` and
  `x`, and `A(s)`, `B(x)`, and `C(s,x)` will be non-species state variables.
- `Catalyst.isequal_ignore_names` has been deprecated for `isequivalent(rn1,
  rn2)` to test equality of two networks and ignore their name. To include names
  in the equality check continue to use `rn1 == rn2` or use `isequivalent(rn1,
  rn2; ignorenames = false)`.

## Catalyst 13.0
- **BREAKING:** Parameters should no longer be listed at the end of the DSL
  macro, but are instead inferred from their position in the reaction statements
  or via explicit declarations in the DSL macro. By default, any symbol that appears
  as a substrate or product is a species, while any other is a parameter. That
  is, parameters are those that only appear within a rate expression and/or as a
  stoichiometric coefficient. E.g. what previously was
  ```julia
  using Catalyst
  rn = @reaction_network begin
    p, 0 --> X
    d, X --> 0
  end p d
  ```
  is now
  ```julia
  using Catalyst
  rn = @reaction_network begin
    p, 0 --> X
    d, X --> 0
  end
  ```
  More generally, in the reaction system
  ```julia
  rn = @reaction_network begin
      k*k1*A, A --> B
      k2, k1 + k*A --> B
  end
  ```
  `k` and `k2` are inferred as parameters by the preceding convention, while
  `A`, `B` and `k1` are species.

- Explicit control over which symbols are treated as parameters vs. species is
  available through the new DSL macros, `@species` and `@parameters`. These can
  be used to designate when something should be a species or parameter,
  overriding the default DSL assignments. This allows setting that a symbol
  which would by default be interpreted as a parameter should actually be a
  species (or vice-versa). E.g. in:
  ```julia
  using Catalyst
  rn = @reaction_network begin
    @species X(t)
    k*X, 0 --> Y
  end
  ```
  `X` and `Y` will be considered species, while `k` will be considered a
  parameter. These options take the same arguments as standalone the `@species`
  (i.e. `ModelingToolkit.@variables`) and `ModelingToolkit.@parameters` macros,
  and support default values and setting metadata. E.g you can set default
  values using:
  ```julia
  using Catalyst
  rn = @reaction_network begin
    @species X(t)=1.0
    @parameters p=1.0 d=0.1
    p, 0 --> X
    d, X --> 0
  end
  ```
  or designate a parameter as representing a constant species using metadata:
  ```julia
  using Catalyst
  rn = @reaction_network begin
    @parameters Y [isconstantspecies=true]
    k, X + Y --> 0
  end
  ```
- **BREAKING:** A standalone `@species` macro was added and should be used in
  place of `@variables` when declaring symbolic chemical species, i.e.
  ```julia
  @parameters k
  @variables t
  @species A(t) B(t)
  rx = Reaction(k, [A], [B])
  @named rs = ReactionSystem([rx], t)
  ```
  This will no longer work as substrates and products must be species
  ```julia
  @parameters k
  @variables t A(t) B(t)
  rx = Reaction(k, [A], [B]) # errors as neither A or B are species
  rx = Reaction(k, [A], nothing) # errors as A is not a species
  rx = Reaction(k, nothing, [B]) # errors as B is not a species

  # this works as the rate or stoichiometry can be non-species
  @species C(t) D(t)
  rx = Reaction(k*A, [C], [D], [2], [B])
  @named rs = ReactionSystem([rx], t)
  ```
  `@variables` is now reserved for non-chemical species state variables (for
  example, arising from constraint equations). Internally, species are normal
  symbolic variables, but with added metadata to indicate they represent
  chemical species.
- To check if a symbolic variable is a species one can use `isspecies`:
  ```julia
  @variables t
  @species A(t)
  @variables B(t)
  isspecies(A) == true
  isspecies(B) == false
  ```
- **BREAKING:** Constraint subsystems and the associated keyword argument to
  `ReactionSystem` have been removed. Instead, one can simply add ODE or
  algebraic equations into the list of `Reaction`s passed to a `ReactionSystem`.
  i.e. this should now work
  ```julia
  @parameters k α
  @variables t V(t)
  @species A(t)
  rx = Reaction(k*V, nothing, [A])
  D = Differential(t)
  eq = D(V) ~ α
  @named rs = ReactionSystem([rx, eq], t)
  osys = convert(ODESystem, rs)
  ```
  which gives the ODE model
  ```
  julia> equations(osys)
  2-element Vector{Equation}:
    Differential(t)(A(t)) ~ k*V(t)
    Differential(t)(V(t)) ~ α
  ```
  Mixing ODEs and algebraic equations is allowed and should work when converting
  to an `ODESystem` or `NonlinearSystem` (if only algebraic equations are
  included), but is not currently supported when converting to `JumpSystem`s or
  `SDESystem`s.
- API functions applied to a `ReactionSystem`, `rs`, now have:

    - `species(rs)` give the chemical species of a system.
    - `states(rs)` give all the variables, both chemical species and
      non-chemical species of a system.

  Catalyst now orders species before non-species in `states(rs)` such that
  `states(rs)[1:length(species(rs))]` and `species(rs)` should be the same.
  Similarly:

    - `equations(rs)` gives the set of `Reaction`s and `Equation`s of a
      system.
    - `reactions(rs)` gives the `Reaction`s of a system.

  As with species, `Reaction`s are always ordered before `Equation`s so that
  `equations(rs)[1:length(reactions(rs))]` should be the same ordered list of
  `Reaction`s as given by `reactions(rs)`.
- Catalyst has been updated for Symbolics v5, and requires Symbolics v5.0.3 or
  greater and ModelingToolkit v8.47.0 or greater.
- The accessors for a given system, `rs`, that return the internal
  arrays at the top-level (i.e. ignoring sub-systems) now have
    - `ModelingToolkit.get_states(rs)` to get the list of all species and
      non-species variables.
    - `Catalyst.get_species(rs)` to get the list of all species variables. Note
      that `get_states(rs)[1:length(get_species(rs))]` should be the same
      ordered list of species as `get_species(rs)`.
    - `ModelingToolkit.get_eqs(rs)` gives the list of all `Reaction`s and then
      `Equation`s in the system.
    - `Catalyst.get_rxs(rs)` gives the list of all `Reaction`s, such that
      `get_eqs(rs)[1:length(get_rx(rs))]` is the same ordered list of
      `Reaction`s as returned by `get_rxs(rs)`.

- **BREAKING:** Chemical species specified or inferred via the DSL are now
  created via the same mechanism as `@species`, and therefore have the
  associated metadata that is missing from a normal symbolic variable.
- Deprecated functions `params` and `merge` have been removed.
- **BREAKING:** The old notation for the constants representing conserved
  quantities, `_Conlaw`, has been replaced with uppercase unicode gamma, "Γ".
  This can be entered in notebooks, the REPL, or many editors by typing the
  corresponding Latex command, "\Gamma", and hitting tab. This leads to much
  cleaner equations when Latexifying systems where conservation laws have been
  applied. The underlying symbol can also be accessed via
  `Catalyst.CONSERVED_CONSTANT_SYMBOL`.
- Modelingtoolkit symbolic continuous and discrete events are now supported when
  creating `ReactionSystem`s via the `continuous_events` and `discrete_events`
  keyword arguments. As in ModelingToolkit, species, states, and parameters that
  appear only within events are not detected automatically, and hence the
  four-argument `ReactionSystem` constructor, where states and parameters are
  explicitly passed, must be used unless every variable, state, or parameter in
  the events appears within a `Reaction` or `Equation` too. See the
  [ModelingToolkit
  docs](https://docs.sciml.ai/ModelingToolkit/stable/basics/Events/) for more
  information on using events. Note that `JumpSystem`s only support discrete
  events at this time.


## Catalyst 12.3.2
- Support for states/species that are functions of multiple variables. This
  enables (symbolically) building PDEs to solve with
  [MethodOfLines](https://github.com/SciML/MethodOfLines.jl/). To use multiple
  independent variables one can say:
  ```julia
  using Catalyst
  using ModelingToolkit: scalarize
  @parameters k[1:7]
  @variables t x y U(x,y,t) V(x,y,t) W(x,y,t)
  rxs = [Reaction(k[1], [U, W], [V, W]),
        Reaction(k[2], [V], [W], [2], [1]),
        Reaction(k[3], [W], [V], [1], [2]),
        Reaction(k[4], [U], nothing),
        Reaction(k[5], nothing, [U]),
        Reaction(k[6], [V], nothing),
        Reaction(k[7], nothing, [V])]
  pars = scalarize(k)
  @named rn = ReactionSystem(rxs, t, [U, V, W], pars; spatial_ivs = [x, y])
  ```
  The `spatial_ivs` keyword lets Catalyst know which independent variables
  correspond to spatial variables. Note that rate expressions can depend on `x`
  and `y` too, i.e. `k[1] * x + y*t` would be valid. See the [work in progress
  PDE
  tutorial](https://github.com/SciML/Catalyst.jl/blob/master/docs/src/tutorials/pdes.md)
  to solve the resulting system and add spatial transport.

## Catalyst 12.3
- API functions to generate substrate, product, and net stoichiometry matrices
  should now work with floating point stoichiometric coefficients. Note,
  symbolic coefficients are still not supported by such functions.

## Catalyst 12.0
- **BREAKING:** Modified how constant and boundary condition species (in the
  SBML sense) work. Constant species should now be specified as ModelingToolkit
  `@parameters` with the `isconstantspecies=true` metadata, while non-constant
  boundary condition species should be specified as ModelingToolkit `@variables`
  with the `isbcspecies=true` metadata. As before, boundary condition species
  are treated as constant with respect to reactions, but since they are
  considered variables their dynamics should be defined in a constraint system.
  Moreover, it is required that BC species appear in a balanced manner (i.e. in
  each reaction for which a BC species is a reactant it must appear as a
  substrate and a product with the same stoichiometry).  Right now only
  conversion of `ReactionSystem`s to an `ODESystem` with a constraint
  `ODESystem` or `NonlinearSystem`, or conversion to a `NonlinearSystem` with a
  constraint `NonlinearSystem`, are supported. Constraints are not supported in
  `SDESystem` or `JumpSystem` conversion, and so boundary condition species are
  effectively constant when converting to those model types (but still left as
  states instead of parameters). Defining constant and boundary condition
  species is done by
  ```julia
  @parameters k A [isconstantspecies=true]
  @variables t  B(t) [isbcspecies=true] C(t)
  rx = Reaction(k, [A,B], [B,C], [1,2], [1,1])
  ```
  Here `A` is a constant species, `B` is a non-constant boundary condition
  species, and `C` is a normal species. Constant and boundary condition species
  can be used in creating `Reaction`s like normal species as either substrates
  or products. Note that network API functions such as `netstoichmat`,
  `conservationlaws`, or `reactioncomplexes` ignore constant species. i.e. for
  `A` a constant species the reaction `2A + B --> C` is treated as equivalent to
  ``B --> C`` with a modified rate constant, while `B --> A` would be identical
  to `B --> 0`. Boundary condition species are checked to be balanced by default
  when `ReactionSystem`s are constructed, i.e.
  ```julia
  rx = Reaction(k, [A,B], [C], [1,2], [1])
  @named rs = ReactionSystem(rs, t)
  ```
  would error since `B` only appears as a substrate. This check can be disabled
  with
  ```julia
  @named rs = ReactionSystem(rs, t; balanced_bc_check=false)
  ```
  Note that network analysis functions assume BC species appear in a balanced
  manner, so may not work correctly if one appears in an unbalanced fashion.
  (Conversion to other system types should still work just fine.)

## Catalyst 11.0
- **BREAKING:** Added the ability to eliminate conserved species when generating
  ODEs, nonlinear problems, SDEs, and steady state problems via the
  `remove_conserved=true` keyword that can be passed to `convert` or to
  `ODEProblem`, `NonlinearProblem`, `SDEProblem`, or `SteadyStateProblem` when
  called with a `ReactionSystem`. For example,
  ```julia
  rn = @reaction_network begin
     k, A + B --> C
     k2, C --> A + B
     end k k2
  osys = convert(ODESystem, rn; remove_conserved=true)
  equations(osys)
  ```
  gives
  ```
  Differential(t)(A(t)) ~ k2*(_ConLaw[2] - A(t)) - k*(A(t) + _ConLaw[1])*A(t)
  ```
  Initial conditions should still be specified for all the species in `rn`, and
  the conserved constants will then be calculated automatically. Eliminated
  species are stored as observables in `osys` and still accessible via solution
  objects. Breaking as this required modifications to the `ReactionSystem` type
  signature.
- **BREAKING:** Added an internal cache in `ReactionSystem`s for network
  properties, and revamped many of the network analysis functions to use this
  cache (so just a `ReactionSystem` can be passed in). Most of these functions
  will now only calculate the chosen property the first time they are called,
  and in subsequent calls will simply returned that cached value. Call
  `reset_networkproperties!` to clear the cache and allow properties to be
  recalculated. The new signatures for `rn` a `ReactionSystem` are
  ```julia
  reactioncomplexmap(rn)
  reactioncomplexes(rn)
  complexstoichmat(rn)
  complexoutgoingmat(rn)
  incidencemat(rn)
  incidencematgraph(rn)
  linkageclasses(rn)
  deficiency(rn)
  sns = subnetworks(rn)
  linkagedeficiencies(rn)
  isreversible(rn)
  isweaklyreversible(rn, sns)
  ```
  Breaking as this required modifications to the `ReactionSystem` type
  signature.
- **BREAKING** `ReactionSystem`s now store a default value for
  `combinatoric_ratelaws=true`. This default value can be set in the
  `ReactionSystem` constructor call as a keyword argument. Passing
  `combinatoric_ratelaws` as a keyword to `convert` or problem calls involving a
  `ReactionSystem` is still allowed, and will override the `ReactionSystem`'s
  default.
- Fixed a bug where `ODESystem` constraint systems did not propagate
  `continuous_events` during calls to `convert(ODESystem, rn::ReactionSystem)`.
- Added constant and boundary condition species (in the SBML sense). During
  conversion constant species are converted to parameters, while boundary
  condition species are kept as state variables. Note that boundary condition
  species are treated as constant with respect to reactions, so their dynamics
  must be defined in a constraint system. Right now only conversion of
  `ReactionSystem`s to an `ODESystem` with a constraint `ODESystem` or
  `NonlinearSystem`, or conversion to a `NonlinearSystem` with a constraint
  `NonlinearSystem`, are supported. Constraints are not supported in `SDESystem`
  or `JumpSystem` conversion, and so boundary condition species are effectively
  constant when converting to those model types (but still left as states
  instead of parameters). Defining constant and boundary condition species is
  done by
  ```julia
  @variables t A(t) [isconstant=true] B(t) [isbc=true] C(t)
  ```
  Here `A` is a constant species, `B` is a boundary condition species, and `C`
  is a normal species. Note that network API functions do not make use of these
  labels, and treat all species as normal -- these properties are only made use
  of when converting to other system types.

## Catalyst 10.8
- Added the ability to use symbolic stoichiometry expressions via the DSL. This should now work
  ```julia
  rn = @reaction_network rs begin
    t*k, (α+k+B)*A --> B
    1.0, α*A + 2*B --> k*C + α*D
  end k α
  ```
  Here Catalyst will try to preserve the order of symbols within an expression,
  taking the rightmost as the species and everything multiplying that species as
  stoichiometry. For example, we can interpret the above reaction as `S1 A -->
  S2 b` where `S1 = (α+k+B)` is the stoichiometry of the reactant `A` and `1` is
  the stoichiometry of the reactant `B`. For
  ```julia
  rn = @reaction_network rs begin
    1.0, 2X*(Y + Z) --> XYZ
  end
  ```
  all of `X`, `Y` and `Z` will be registered as species, with substrates `(Y,Z)` having associated stoichiometries of
  `(2X,2X)`. As for rate expressions, any symbols that appear and are not defined as parameters will be declared to be species.

  In contrast, when declaring reactions
  ```julia
  rx = @reaction t*k, (k+α)*A --> B
  ```
  will work, with every symbol declared a parameter except the leftmost symbol in the reaction line. So
  ```julia
  rx = @reaction 1.0, 2X*(Y + Z) --> XYZ
  ```
  will make `X` a parameter and `Y`, `Z` and `XYZ` species.
- Symbolic stoichiometry supports interpolation of expressions in
  `@reaction_network` and `@reaction`.

## Catalyst 10.7
- Added the ability to use symbolic variables, parameters and expressions for
  stoichiometric coefficients. See the new tutorial on [Parametric
  Stoichiometry](https://docs.sciml.ai/Catalyst/stable/tutorials/symbolic_stoich/) for
  details, and note the caveat about ModelingToolkit converting integer
  parameters to floating point types that must be worked around to avoid calls
  to `factorial` that involve `float`s.

## Catalyst 10.6
- Added the ability to use floating point stoichiometry (currently only tested
  for generating ODE models). This should now work
  ```julia
  rn = @reaction_network begin
    k, 2.5*A --> 3*B
  end k
  ```
  or directly
  ```julia
  @parameters k b
  @variables t A(t) B(t) C(t) D(t)
  rx1 = Reaction(k,[B,C],[B,D], [2.5,1],[3.5, 2.5])
  rx2 = Reaction(2*k, [B], [D], [1], [2.5])
  rx3 = Reaction(2*k, [B], [D], [2.5], [2])
  @named mixedsys = ReactionSystem([rx1,rx2,rx3],t,[A,B,C,D],[k,b])
  osys = convert(ODESystem, mixedsys; combinatoric_ratelaws=false)
  ```
  Note, when using `convert(ODESystem, mixedsys; combinatoric_ratelaws=false)`
  the `combinatoric_ratelaws=false` parameter must be passed. This is also true
  when calling `ODEProblem(mixedsys,...; combinatoric_ratelaws=false)`. This
  disables Catalyst's standard rescaling of reaction rates when generating
  reaction rate laws, see the
  [docs](https://catalyst.sciml.ai/dev/tutorials/introduction_to_catalyst/#Reaction-rate-laws-used-in-simulations).
  Leaving this out for systems with floating point stoichiometry will give an
  error message.

## Catalyst 10.5
- Added `@reaction` macro
  ```julia
  rx = @reaction k*v, A + B --> C + D

  # is equivalent to
  @parameters k v
  @variables t A(t) B(t) C(t) D(t)
  rx == Reaction(k*v, [A,B], [C,D])
  ```
  Here `k` and `v` will be parameters and `A`, `B`, `C` and `D` will be
  variables. Interpolation of existing parameters/variables also works
  ```julia
  @parameters k b
  @variables t A(t)
  ex = k*A^2 + t
  rx = @reaction b*$ex*$A, $A --> C
  ```
  Any symbols arising in the rate expression that aren't interpolated are
  treated as parameters, while any in the reaction part (`A + B --> C + D`) are
  treated as species.

## Catalyst 10.4
- Added `symmap_to_varmap`, `setdefaults!`, and updated all `*Problem(rn,...)`
  calls to allow setting initial conditions and parameter values using symbol
  maps. See the [Catalyst API](https://docs.sciml.ai/Catalyst/stable/) for details.
  These allow using regular Julia `Symbols` to specify parameter values and
  initial conditions. i.e. to set defaults we can do
  ```julia
  rn = @reaction_network begin
      α, S + I --> 2I
      β, I --> R
  end α β
  setdefaults!(rn, [:S => 999.0, :I => 1.0, :R => 0.0, :α => 1e-4, :β => .01])
  op    = ODEProblem(rn, [], (0.0,250.0), [])
  sol   = solve(op, Tsit5())
  ```
  To explicitly pass initial conditions and parameters using symbols we can do
  ```julia
  rn = @reaction_network begin
      α, S + I --> 2I
      β, I --> R
  end α β
  u0 = [:S => 999.0, :I => 1.0, :R => 0.0]
  p  = (:α => 1e-4, :β => .01)
  op    = ODEProblem(rn, u0, (0.0,250.0), p)
  sol   = solve(op, Tsit5())
  ```
  In each case ModelingToolkit symbolic variables can be used instead of
  `Symbol`s, e.g.
  ```julia
  @parameters α β
  @variables t S(t) I(t) R(t)
  setdefaults!(rn, [S => 999.0, I => 1.0, R => 0.0, α => 1e-4, β => .01])
  ```

## Catalyst 10.3
- **BREAKING:** The order of the parameters in the `ReactionSystem`'s `.ps`
  field has been changed (only when created through the `@reaction_network`
  macro). Previously they were ordered according to the order with which they
  appeared in the macro. Now they are ordered according the to order with which
  they appeared after the `end` part. E.g. in
  ```julia
  rn = @reaction_network begin
    (p,d), 0 <--> X
  end d p
  ```
  previously the order was `[p,d]`, while now it is `[d, p]`.


## Catalyst 10.1
- Added support for `@unpack observable_variable = rn` and
  `rn.observable_variable`. This requires a new inner constructor definition for
  `ReactionSystem`s, but is not considered breaking as the inner constructor is
  considered private.
- Support added for ModelingToolkit 7 and Symbolics 4.

## Catalyst 10.0
- `ReactionSystem(rxs::Vector{Reaction}, t)` should now work and will infer the
  species and parameters.
- **BREAKING:** Any undeclared variables in the DSL are now inferred to be
  species. i.e. this no longer errors, and `B` is assumed to be a species
  ```julia
  rn = @reaction_network begin
    k*B, A --> C
  end k
  ```
- **BREAKING:** Internal changes mean the order of species or parameters in
  generated systems may have changed. Changes that induce different orders will
  not be considered breaking in the future.
- Added interpolation in the DSL for species, variables, and the network name.
  i.e. this is now valid
  ```julia
  @parameters k
  @variables t, A(t)
  spec = A
  rate = k*A
  name = :network
  rn = @reaction_network $name begin
    $rate*B, 2*$spec + B --> $spec + C
    end
  ```
- Added the ability to compose `ReactionSystem`s via subsystems, and include
  either `ODESystem`s or `NonlinearSystem`s as subsystems. Note, if using
  non-`ReactionSystem` subsystems it is not currently possible to convert to
  a `JumpSystem` or `SDESystem`. It is also not possible to include either
  `SDESystem`s or `JumpSystems` as subsystems.
- Added `extend(sys, reactionnetwork, name=nameof(sys))` to extend
  `ReactionSystem`s with constraint equations (algebraic equations or ODEs), or
  other `ReactionSystem`s. Algebraic or differential constraints are stored as a `NonlinearSystem` or
  `ODESystem` within the `ReactionSystem`, and accessible via
  `get_constraints(reactionnetwork)`.
- Added `Catalyst.flatten(rn)` to allow flattening of a `ReactionSystem` with
  sub-systems into one `ReactionSystem`. Non-`ReactionSystem` subsystems are
  merged into the constraints of the flattened `ReactionSystem`, and accessible
  via `get_constraints`.
- **BREAKING:** `ReactionSystem`s are now always flattened when calling
  `convert`. This should only affect models that use `subsystem`s.
- Added `incidencematgraph`, `linkageclasses`, `deficiency`, `subnetworks`,
  `linkagedeficiency`, `isreversible` and `isweaklyreversible` API functions.
- Deprecated `merge`, use `ModelingToolkit.extend` instead.
- Deprecated `params` and `numparams` (use `ModelingToolkit.parameters` to get
  all parameters of a system and all subsystems, or use `reactionparams` to get
  all parameters of a system and all `ReactionSystem` subsystems. The latter
  correspond to those parameters used within `Reaction`s.)
- **BREAKING:** Added a custom `hash` for `Reaction`s to ensure they work in
  `Dict`s and `Set`s properly, ensuring set-type comparisons between collections
  of `Reaction`s work.
- Updated the docs and added a new tutorial on using compositional tooling.

## Catalyst 9.0
*1.* **BREAKING:** `netstoichmat`, `prodstoichmat` and `substoichmat` are now
transposed to be number of species by number of reactions. This is more
consistent with the chemical reaction network literature for stoichiometry
matrices.

*2.* `reactioncomplexmap` added to provide a mapping from reaction complexes to
reactions they participate in.

*3.* Most API `*mat` functions now take an optional `sparse` keyword argument.
If passed `sparse=true` a sparse matrix representation is generated, otherwise
the default `sparse=false` value returns dense `Matrix` representations.

## Catalyst 8.3
*1.* Network representations for the reaction complexes of a system along with
associated graph functionality:
```julia
rn = @reaction_network begin
           k₁, 2A --> B
           k₂, A --> C
           k₃, C --> D
           k₄, B + D --> E
           k₅, B --> E
           k₆, D --> C
     end k₁ k₂ k₃ k₄ k₅ k₆
smap  = speciesmap(rn)
rcs,B = reactioncomplexes(rn; smap=smap)
Z     = complexstoichmat(rn; rcs=rcs)
Δ     = complexoutgoingmat(rn; B=B)
complexgraph(rn; complexdata=(rcs,B))
```
which gives

![rn_complexes](https://user-images.githubusercontent.com/9385167/130252763-4418ba5a-164f-47f7-b512-a768e4f73834.png)

*2.* Support for units via ModelingToolkit and
[Unitful.jl](https://github.com/PainterQubits/Unitful.jl) in directly constructed
`ReactionSystem`s:
```julia
# ]add Unitful
using Unitful
@parameters α [unit=u"μM/s"] β [unit=u"s"^(-1)] γ [unit=u"μM*s"^(-1)]
@variables t [unit=u"s"] A(t) [unit=u"μM"] B(t) [unit=u"μM"] C(t) [unit=u"μM"]
rxs = [Reaction(α, nothing, [A]),
       Reaction(β, [A], [B]),
       Reaction(γ, [A,B], [B], [1,1], [2])]
@named rs = ReactionSystem(rxs, t, [A,B,C], [α,β,γ])
```
By default, during construction of `rs` Catalyst will call
```julia
validate(rs)
```
which will print warnings and return `false` if either
1. The `species(rs)` do not all have the same units.
2. The implicit (ODE) rate laws for each reaction do not have units of (species
   units) / (time units), where the time units are the units of `t`.

(Note, at this time the `@reaction_network` macro does not support units.)

*3.* Calculation of conservation laws
```julia
rn = @reaction_network begin
  (k₊,k₋), A + B <--> C
  end k₊ k₋
clawmat = conservationlaws(netstoichmat(rn))
```
giving
```
 1  -1  0
 0   1  1
```
and
```julia
cquants = conservedquantities(species(rn), clawmat)
```
giving
```
 A(t) - B(t)
 B(t) + C(t)
```

See the [API docs](https://docs.sciml.ai/Catalyst/stable/api/catalyst_api/) for more
details about each of these new features.

## Catalyst 8.2
*1.* Basic unit validation has been added following its addition for all
ModelingToolkit systems.

## Catalyst 8.1
*1.* `reactioncomplexes`, `ReactionComplex`, `reactionrates`, `complexstoichmat`
and `complexoutgoingmat` are added to allow the calculation of reaction complex-based
network matrix representations.

## Catalyst 8.0
**BREAKING:** This is a breaking release, with all ModelingToolkit `ReactionSystem` and
`Reaction` functionality migrated to Catalyst.

## Catalyst 6.11
*1.* Plain text arrows "<--" and "<-->" for backward and reversible reactions are
   available if using Julia 1.6 or higher:
```julia
rn = @reaction_network begin
  (k1,k2), A + B <--> C
  k3, 0 <-- C
end k1 k2 k3
```
*2.* **BREAKING:** Reaction networks can be named
```julia
rn = @reaction_network Reversible_Reaction begin
  k1, A --> B
  k2, B --> A
  end k1 k2
ModelingToolkit.nameof(rn) == :Reversible_Reaction
```
Note, empty networks can no longer be created with parameters, i.e. only
```julia
rn = @reaction_network          # uses a randomly generated name
rn = @reaction_network MyName   # is named MyName
```
are allowed.

*3.* Compositional modeling with generated `ODESystem`s, see
[here](https://github.com/SciML/Catalyst.jl/blob/master/test/reactionsystem_components.jl)
for an example that composes three gene modules to make the repressilator.
