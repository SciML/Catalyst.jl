# [Spatial jump simulations](@id spatial_lattice_jump_simulations)
Our [introduction to spatial lattice simulations](@ref spatial_lattice_modelling_intro) has already described how to simulate `LatticeReactionSystem`s using ODEs. Jump simulations of `LatticeReactionSystem` are carried out using an almost identical approach. However, just like for non-spatial models, we must first create a `DiscreteProblem`, which is then used as input to our `JumpProblem`. Furthermore, a spatial [jump aggregator](@ref simulation_intro_jumps_solver_designation) (like `NSM`) can be used. Spatial jump simulations in Catalyst are built on top of JumpProcesses.jl's spatial jump simulators, more details on which can be found [here](https://docs.sciml.ai/JumpProcesses/stable/tutorials/spatial/).

Below we perform a spatial jump simulation of a simple [birth-death process](@ref basic_CRN_library_bd). Note that we use our `LatticeReactionSystem` as input to both our `DiscreteProblem` and `JumpProblem`, and that we [provide] the spatial `NSM` jump aggregator to `JumpProblem`.
```@example spatial_jump
using Catalyst, JumpProcesses
bd_model = @reaction_network begin
 (p,d), 0 <--> X
end
diffusion_rx = @transport_reaction D X
lattice = CartesianGrid((10,10))
lrs = LatticeReactionSystem(bd_model, [diffusion_rx], lattice)

u0 = [:X => 0]
tspan = (0.0, 200.0)
ps = [:p => 10.0, :d => 1.0, :D => 0.1]
dprob = DiscreteProblem(lrs, u0, tspan, ps)
jprob = JumpProblem(lrs, dprob, NSM())
sol = solve(jprob, SSAStepper())
nothing # hide
```
We can now access the values of `sol` using the interfaces described [here](@ref lattice_simulation_structure_interaction_simulation_species), and plot it using the functions described [here](@ref lattice_simulation_plotting).

Currently, the only available spatial jump aggregators are `NSM` and `DirectCRDirect`, with `DirectCRDirect` expected to perform better for large networks. It is also possible to use non-spatial solvers (like `Direct` and `RSSA`) when creating spatial `JumpProblem`s.