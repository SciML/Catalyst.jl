# [Solving the chemical master equation using FiniteStateProjection.jl](@id finite-state_projection)

## [Finite state projection simulation of single-species model](@id state_projection_one_species)
For this example we will use a simple [birth-death model](@ref basic_CRN_library_bd), where a single species ($X$) is created and degraded at constant rates ($p$ and $d$, respectively).
```@example state_projection_one_species
using Catalyst
rs = @reaction_network begin
    (p,d), 0 <--> X
end
```
Next we perform jump simulations (using the [ensemble simulation interface](@ref ensemble_simulations_monte_carlo)) of the model. Here, we can see how the it develops from an initial condition and reaches a steady state distribution.
```@example state_projection_one_species
using JumpProcesses, Plots
u0 = [:X => 5]
tspan = (0.0, 10.0)
ps = [:p => 20.0, :d => 0.5]
jprob = JumpProblem(JumpInputs(rs, u0, tspan, ps))
eprob = EnsembleProblem(jprob)
esol = solve(eprob, SSAStepper(); trajectories = 10)
plot(esol; ylimit = (0.0, Inf))
```
Using chemical master equation simulations, we want to simulate how the *full probability distribution* of these jump simulations develop across the simulation time-frame. 

As a first step, we import the FiniteStateProjection package. Next, we convert our [`ReactionSystem`](@ref) to a `FSPSystem` (from which we later will generate an ODE).
```@example state_projection_one_species
using FiniteStateProjection
fsp_sys = FSPSystem(rs)
nothing # hide
```
Next, we set our initial condition. For normal simulations, $X$'s initial condition would be a single value. Here, however, we will simulate $X$'s probability distribution. Hence, its initial condition will also be a probability distribution. In FiniteStateProjection's interface, the initial condition is a vector, where the $i$'th index is the probability that $X$ have an initial value of $i-1$. The total sum of all probabilities across the vector should be normalised to $1$. Here we assume that $X$'s initial conditions is known to be $5$ (hence the corresponding probability is $1$, and the remaining ones are $0$):
```@example state_projection_one_species
u0 = zeros(75)
u0[6] = 1.0
bar(u0)
```
We also plot the full distribution using the `bar` function. Finally, the initial condition vector also define the finite space onto which we project the chemical master equation. I.e. we will assume that, throughout the entire simulation, the probability of $X$ reaching values outside this initial vector is negligible. 
!!! warn
    This last fact is important. Even if the probability seems to ve very small on the boundary provided by the initial condition, there is still a risk that probability will "leak". Here, it can be good to make simulations using different projections, ensuring that the results are consistent (especially for longer simulations).

Now, we can finally create an `ODEProblem` using out `FSPSystem`, our initial conditions, and the parameters declared previously. We can simulate this `ODEProblem` like any other ODE.
```@example state_projection_one_species
using OrdinaryDiffEqDefault
oprob = ODEProblem(fsp_sys, u0, tspan, ps)
osol = solve(oprob)
nothing # hide
```
Finally, we can plot the $X$'s probability distribution at various time point's of the simulation. Again, we will uset he `bar` function to plot teh distribution, and the interface described [here](@ref simulation_structure_interfacing_solutions) to acess the simulation at various timestamps.
```@example state_projection_one_species
bar(osol(1.0);  width = 1.0, linecolor = 1, alpha = 0.7, linealpha = 0.7, label = "t = 1.0")
bar!(osol(2.0); linecolor = 1, alpha = 0.7, linealpha = 0.7, label = "t = 2.0")
bar!(osol(5.0); linecolor = 1, alpha = 0.7, linealpha = 0.7, label = "t = 5.0")
bar!(osol(10.0); linecolor = 1, alpha = 0.7, linealpha = 0.7, label = "t = 10.0")
```

## [Finite state projection simulation of multi-species model](@id state_projection_multi_species)

```@example state_projection_multi_species
using Catalyst # hide
rs = @reaction_network begin
    (p,d), 0 <--> X
    (kB,kD), 2X <--> X2
end
```

```@example state_projection_multi_species
ps = [:p => 1.0, :d => 0.2, :kB => 2.0, :kD => 5.0]
u0 = zeros(25,25)
u0[1,1] = 1.0
```

```@example state_projection_multi_species
using OrdinaryDiffEqRosenbrock, Plots # hide
fsp_sys = FSPSystem(rs, [:X, :X2])
oprob = ODEProblem(fsp_sys, u0, 100.0, ps)
@time osol = solve(oprob, Rodas5P())
heatmap(osol[end]; xguide = "X2", yguide = "X")
```

## [Finite state projection steady state simulations](@id state_projection_steady_state_sim)
Previously we described how the [SteadyStateDiffEq.jl](https://github.com/SciML/SteadyStateDiffEq.jl) package can be used to [find an ODE's steady state through forward simulation](@ref steady_state_stability). The same interface can be used for ODEs generated through FiniteStateProjection.jl. Below we use this to find the steady state of the dimerisation example studied in the last example.
```@example state_projection_multi_species
using SteadyStateDiffEq, OrdinaryDiffEqRosenbrock
ssprob = SteadyStateProblem(fsp_sys, u0, ps)
sssol = solve(ssprob, DynamicSS(Rodas5P()))
heatmap(sssol; xguide = "X2", yguide = "X")
```