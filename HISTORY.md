# Breaking updates and feature summaries across releases

## Catalyst unreleased (master branch) 
- **BREAKING:** The order of the parameters in the `ReactionSystem`'s `.ps` field has been changed (only when created through the `@reaction_network` macro). Previously they were ordered acording to the order with which they appeared in the macro. Now they are ordered acoridng the to order with which the appeard after the `end` part. E.g. in
  ```julia
  rn = @reaction_network begin
    (p,d), 0 <--> X
  end d p
  ```
  previousl the order was `[p,d]`, which now it is `[d, p]`.


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
