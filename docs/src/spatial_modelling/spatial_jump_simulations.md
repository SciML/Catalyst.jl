# [Spatial Jump Simulations](@id spatial_dspace_jump_simulations)
Our [introduction to spatial space simulations](@ref spatial_dspace_modelling_intro) has already described how to simulate [`DiscreteSpaceReactionSystem`](@ref)s using ODEs. Jump simulations of [`DiscreteSpaceReactionSystem`](@ref) are carried out using an almost identical approach. Spatial jump simulations in Catalyst are built on top of JumpProcesses.jl's spatial jump simulators, more details on which can be found [here](https://docs.sciml.ai/JumpProcesses/stable/tutorials/spatial/).

Below we perform a spatial jump simulation of a simple [birth-death process](@ref basic_CRN_library_bd).
```@example spatial_jump
using Catalyst, JumpProcesses
bd_model = @reaction_network begin
    (p,d), 0 <--> X
end
diffusion_rx = @transport_reaction D X
space = CartesianGrid((10,10))
dsrs = DiscreteSpaceReactionSystem(bd_model, [diffusion_rx], space)

u0 = [:X => 0]
tspan = (0.0, 200.0)
ps = [:p => 10.0, :d => 1.0, :D => 0.1]
dprob = DiscreteProblem(dsrs, u0, tspan, ps)
jprob = JumpProblem(dsrs, u0, tspan, ps)
sol = solve(jprob, SSAStepper())
nothing # hide
```
We can now access the values of `sol` using the interfaces described [here](@ref dspace_simulation_structure_interaction_simulation_species), and plot it using the functions described [here](@ref dspace_simulation_plotting).

Aggregators can be designated through the `aggregator` key word argument, By default, `NSM()` is used. However, `DirectCRDirect()` is also available and expected to perform better for large networks.

!!! note
    Currently, spatial jump simulations are only supported when all reaction of the non-spatial `ReactionSystem` are `MassActionJump`s, i.e. have constant rates. This means that reactions with e.g. Michaelis-Menten rates are currently not supported.