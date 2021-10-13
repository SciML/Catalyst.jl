# Compositional Modeling of Reaction Systems
Catalyst supports the construction of models in a compositional fashion, based
on ModelingToolkit's subsystem functionality. In this tutorial we'll see how we
can construct the earlier repressilator model by composing together three
identically repressed genes, and how to use compositional modeling to create
compartments.

## Compositional Modeling Tooling
Catalyst supports two ModelingToolkit interfaces for composing multiple
[`ReactionSystem`](@ref)s together into a full model. The first mechanism for
extending a system is the `extend` command
```@example ex1
using Catalyst
basern = @reaction_network rn1 begin
           k, A + B --> C
         end k
newrn = @reaction_network rn2 begin
        r, C --> A + B
      end r
@named rn = extend(newrn, basern)
show(stdout, MIME"text/plain"(), rn) # hide
```
with reactions
```julia
reactions(rn)
```
```@example ex1
show(stdout, MIME"text/plain"(), reactions(rn)) # hide
```
Here we extended `basern` with `newrn` giving a system with all the
reactions. Note, if a name is not specified via `@named` or the `name` keyword
then `rn` will have the same name as `newrn`.

The second main compositional modeling tool is the use of subsystems. Suppose we
now create `rn` with `basern` and `newrn` as subsystems of it, we get a
different result:
```@example ex1
newestrn = @reaction_network rn3 begin
            v, A + D --> 2D
           end v
@named rn = compose(basern, [newrn, newestrn])
show(stdout, MIME"text/plain"(), rn) # hide
```
with reactions
```julia
reactions(rn)
```
```@example ex1
show(stdout, MIME"text/plain"(), reactions(rn)) # hide
```
Here we have created a new `ReactionSystem` that adds `newrn` and `newestrn` as
subsystems of `basern`. The variables and parameters in the sub-systems are
considered distinct from those in other systems, and so are namespaced (i.e.
prefaced) by the name of the system they come from. 

We can see the subsystems of a given system like
```@example ex1
ModelingToolkit.get_systems(rn)
```

## Compositional Models for the Repressilator
