# Breaking updates and feature summaries across releases

## Catalyst unreleased (master branch)
- **BREAKING:** Added the ability to eliminate conserved species when generating
  ODEs, nonlinear problems, and steady-state problems via the
  `remove_conserved=true` keyword that can be passed to `convert` or to
  `ODEProblem`, `NonlinearProblem` or `SteadyStateProblem` when called with a
  `ReactionSystem`. For example,
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
- **BREAKING:** Added an internal cache in `ReactionSystem`s for network properties, and
  revamped many of the network analysis functions to use this cache (so just a
  `ReactionSystem` can be passed in). Most of these functions will now only
  calculate the chosen property the first time they are called, and in
  subsequent calls will simply returned that cached value. Call
  `reset_networkproperties!` to clear the cache and allow properties to be
  recalculated. The new signatures for `rn` a `ReactionSystem` are
  ```julia
  reactioncomplexmap(rn)
  reactioncomplexes(rn)
  complexstoichmat(rn)
  complexoutgoingmat(rn)
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
  instead of parameters). Defining constant and boundary condition species is done by
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
  Stoichiometry](https://catalyst.sciml.ai/dev/tutorials/symbolic_stoich/) for
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
  [docs](https://catalyst.sciml.ai/dev/tutorials/using_catalyst/#Reaction-rate-laws-used-in-simulations).
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
  maps. See the [Catalyst API](https://catalyst.sciml.ai/dev/) for details.
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
  they appeard after the `end` part. E.g. in
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
[Uniftul.jl](https://github.com/PainterQubits/Unitful.jl) in directly constructed
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

See the [API docs](https://catalyst.sciml.ai/dev/api/catalyst_api/) for more
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
nameof(rn) == :Reversible_Reaction
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
