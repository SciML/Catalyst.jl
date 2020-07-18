# Catalyst.jl API
```@meta
CurrentModule = Catalyst
```

## Reaction Network Generation and Representation
Catalyst provides the [`@reaction_network`](@ref) macro for generating a
complete network, stored as a [`ModelingToolkit.ReactionSystem`](@ref), which in
turn is composed of [`ModelingToolkit.Reaction`](@ref)s. `ReactionSystem`s can
be converted to other `ModelingToolkit.AbstractSystem`s, including a
`ModelingToolkit.ODESystem`, `ModelingToolkit.SDESystem` or
`ModelingToolkit.JumpSystem`. 

An empty network can be generated using [`@reaction_network`](@ref) with no arguments or
the [`make_empty_network`](@ref) function. These can then be extended
programmatically using [`addspecies!`](@ref), [`addparam!`](@ref), and
[`addreaction!`](@ref). 

It is important to note for [`@reaction_network`](@ref) that species which are
used *within the macro* as part of a rate expression, but not as a substrate or
product of some reaction, may lead to undefined behavior. i.e. avoid
```julia
rn = @reaction_network begin
    k*X, Y --> W
end k
```
as here `X` is never defined as either a species or parameter. 

```@docs
@reaction_network
make_empty_network
@add_reactions
ModelingToolkit.Reaction
ModelingToolkit.ReactionSystem
```

## Basic properties
```@docs
species
speciesmap
params
paramsmap 
reactions
numspecies 
numparams
numreactions
```

## Reaction Properties
```@docs
ModelingToolkit.ismassaction
dependents
dependants
```

## Functions to extend a Network
```@docs
addspecies!
addparam!
addreaction!
merge!(network1::ReactionSystem, network2::ReactionSystem)
merge(network1::ReactionSystem, network2::ReactionSystem)
```

## Generated ModelingToolkit Operations
As the underlying [`ReactionSystem`](@ref) is comprised of
`ModelingToolkit.Operation`s and `ModelingToolkit.Variable`s, one can directly
access the generated rate laws, and using `ModelingToolkit` tooling generate
functions or Julia `Expr`s from them.
```@docs
ModelingToolkit.oderatelaw
ModelingToolkit.jumpratelaw
```

## Network Comparison Functions
```@docs
==(rn1::ReactionSystem, rn2::ReactionSystem)
==(rn1::Reaction, rn2::Reaction)
```
