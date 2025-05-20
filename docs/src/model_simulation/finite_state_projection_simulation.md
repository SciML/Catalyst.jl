# [Solving the chemical master equation using FiniteStateProjection.jl](@id finite-state_projection)
```@raw html
<details><summary><strong>Quick-start example</strong></summary>
```
The following code provides a brief example of how to simulate the chemical master equation using the [FiniteStateProjection.jl](https://github.com/SciML/FiniteStateProjection.jl) package.
```julia
# Create reaction network model (a bistable switch).
using Catalyst
rs_bistable = @reaction_network begin
    hillr(Y, v, K, n), ∅ --> X
    hillr(X, v, K, n), ∅ --> Y
 d, (X,Y) --> 0
end

# Create a FSPSystem. The second argument denotes species order in u0 and sol.
using FiniteStateProjection
fsp_sys = FSPSystem(rs_bistable, [:X, :Y])

# Create ODEProblem. Initial condition is the system's initial distribution.
# Initial condition also designates projection space (outside of which solution should be very close to 0).
u0 = zeros(20,20)
u0[6,2] = 1.0
tspan = (0.0, 50.0)
ps = [:v => 1.0, :K => 3.0, :n => 3, :d => 0.2]
oprob = ODEProblem(fsp_sys, u0, tspan, ps)

# Simulate ODE (it can be quite large, so consider performance options).
# Plot solution as a heatmap at a specific time point.
using OrdinaryDiffEqRosenbrock, Plots
osol = solve(oprob, Rodas5P())
heatmap(osol(50.0); xguide = "Y", yguide = "X")
```
```@raw html
</details>
```
  \
  
Previously, we have shown how [*stochastic chemical kinetics*](@ref math_models_in_catalyst_sck_jumps) describe how chemical reaction network models can be [exactly simulated](@ref simulation_intro_jumps) (using e.g. [Gillespie's algorithm](https://en.wikipedia.org/wiki/Gillespie_algorithm)). We also described how the [SDE](@ref math_models_in_catalyst_cle_sdes) and [ODE](@ref math_models_in_catalyst_rre_odes) approaches were approximations of these jump simulations, and only valid for large copy numbers. To gain a good understanding of the system's time development, we typically have to carry out a large number of jump simulations. An alternative approach, however, is to instead simulate the *full probability distribution of the system*. This corresponds to the distribution from which these jump simulations are drawn.

[*The chemical master equation*](https://en.wikipedia.org/wiki/Master_equation) (CME) describes the time development of this probability distribution[^1], and is given by a (possibly infinite) coupled system of ODEs (with one ODE for each possible chemical state, i.e. number configuration, of the system). For a system with a single species $X$, the CME looks like
```math
\begin{aligned}
\frac{dp(x=0)}{dt} &= f_0(p(x=0), p(x=1), ...) \\
\frac{dp(x=1)}{dt} &= f_0(p(x=0), p(x=1), ...) \\
\frac{dp(x=2)}{dt} &= f_0(p(x=0), p(x=1), ...) \\
 &\vdots\\
\end{aligned}
```
We note that, since (for almost all chemical reaction networks) there is a non-zero probability that $X$ is at any specific integer value, the CME have an infinite number of variables (corresponding to the infinite number of potential states). Hence, it cannot be solved practically. However, for high enough species values, the probability of the system attaining such values becomes negligibly small. Here, a truncated version of the CME can be solved practically. An approach for this is the *finite state projection*[^2]. Below we describe how to use the [FiniteStateProjection.jl](https://github.com/SciML/FiniteStateProjection.jl) package to solve the truncated CME (with the package's [documentation](https://docs.sciml.ai/FiniteStateProjection/dev/) providing a more extensive description). While the CEM approach can be very powerful, we note that even for systems with few species, the truncated CME typically have too many states to be feasibly solved.

## [Finite state projection simulation of single-species model](@id state_projection_one_species)
For this example, we will use a simple [birth-death model](@ref basic_CRN_library_bd), where a single species ($X$) is created and degraded at constant rates ($p$ and $d$, respectively).
```@example state_projection_one_species
using Catalyst
rs = @reaction_network begin
    (p,d), 0 <--> X
end
```
Next, we perform jump simulations (using the [ensemble simulation interface](@ref ensemble_simulations_monte_carlo)) of the model. Here, we can see how it develops from an initial condition and reaches a steady state distribution.
```@example state_projection_one_species
using JumpProcesses, Plots
u0 = [:X => 5]
tspan = (0.0, 10.0)
ps = [:p => 20.0, :d => 0.5]
jprob = JumpProblem(JumpInputs(rs, u0, tspan, ps))
eprob = EnsembleProblem(jprob)
esol = solve(eprob, SSAStepper(); trajectories = 10)
plot(esol; ylimit = (0.0, Inf))
Catalyst.PNG(plot(esol; ylimit = (0.0, Inf), fmt = :png, dpi = 200)) # hide
```
Using chemical master equation simulations, we want to simulate how the *full probability distribution* of these jump simulations develops across the simulation time frame. 

As a first step, we import the FiniteStateProjection package. Next, we convert our [`ReactionSystem`](@ref) to a `FSPSystem` (from which we later will generate an ODE).
```@example state_projection_one_species
using FiniteStateProjection
fsp_sys = FSPSystem(rs)
nothing # hide
```
Next, we set our initial condition. For normal simulations, $X$'s initial condition would be a single value. Here, however, we will simulate $X$'s probability distribution. Hence, its initial condition will also be a probability distribution. In FiniteStateProjection's interface, the initial condition is an array, where the $i$'th index is the probability that $X$ have an initial value of $i-1$. The total sum of all probabilities across the vector should be normalised to $1.0$. Here we assume that $X$'s initial conditions is known to be $5$ (hence the corresponding probability is $1.0$, and the remaining ones are $0.0$):
```@example state_projection_one_species
u0 = zeros(75)
u0[6] = 1.0
bar(u0, label = "t = 0.0")
```
We also plot the full distribution using the `bar` function. Finally, the initial condition vector defines the finite space onto which we project the CME. I.e. we will assume that, throughout the entire simulation, the probability of $X$ reaching values outside this initial vector is negligible. 

!!! warning
    This last bit is important. Even if the probability seems to be very small on the boundary provided by the initial condition, there is still a risk that probability will "leak". Here, it can be good to make simulations using different projections, ensuring that the results are consistent (especially for longer simulations).

Now, we can finally create an `ODEProblem` using our `FSPSystem`, initial conditions, and the parameters declared previously. We can simulate this `ODEProblem` like any other ODE.
```@example state_projection_one_species
using OrdinaryDiffEqDefault
oprob = ODEProblem(fsp_sys, u0, tspan, ps)
osol = solve(oprob)
nothing # hide
```
Finally, we can plot $X$'s probability distribution at various simulation time points. Again, we will use the `bar` function to plot the distribution, and the interface described [here](@ref simulation_structure_interfacing_solutions) to access the simulation at specified time points.
```@example state_projection_one_species
bar(osol(1.0);  bar_width = 1.0, linewidth = 0, alpha = 0.7, label = "t = 1.0")
bar!(osol(2.0); bar_width = 1.0, linewidth = 0, alpha = 0.7, label = "t = 2.0")
bar!(osol(5.0); bar_width = 1.0, linewidth = 0, alpha = 0.7, label = "t = 5.0")
bar!(osol(10.0); bar_width = 1.0, linewidth = 0, alpha = 0.7, label = "t = 10.0",
  xguide = "X (copy numbers)", yguide = "Probability density")
```

## [Finite state projection simulation of multi-species model](@id state_projection_multi_species)
Next, we will consider a system with more than one species. The workflow will be identical, however, we will have to make an additional consideration regarding our initial conditions, simulation performance, and plotting approach.

For this example, we will consider a simple dimerisation model. In it, $X$ gets produced and degraded at constant rates, and can also dimerise to form $X₂$.
```@example state_projection_multi_species
using Catalyst # hide
rs = @reaction_network begin
    (p,d), 0 <--> X
    (kB,kD), 2X <--> X₂
end
```
Next, we will declare our parameter values and initial condition. In this case, the initial condition is a matrix where element $(i,j)$ denotes the initial probability that $(X(0),X₂(0)) = (i-1,j-1)$. In this case, we will use an initial condition where we know that $(X(0),X₂(0)) = (5,0)$.
```@example state_projection_multi_species
ps = [:p => 1.0, :d => 0.2, :kB => 2.0, :kD => 5.0]
u0 = zeros(25,25)
u0[6,1] = 1.0
nothing # hide
```
In the next step, however, we have to make an additional consideration. Since we have more than one species, we have to define which dimension of the initial condition (and hence also the output solution) corresponds to which species. Here we provide a second argument to `FSPSystem`, which is a vector listing all species in the order they occur in the `u0` array.
```@example state_projection_multi_species
using FiniteStateProjection # hide
fsp_sys = FSPSystem(rs, [:X, :X₂])
nothing # hide
```
Finally, we can simulate the model just like in the 1-dimensional case. In this case, however, we are simulating an ODE with $25⋅25 = 625$ states, which means we need to make some considerations regarding performance. In this case, we will simply specify the `Rodas5P()` ODE solver (more extensive advice on performance can be found [here](@ref ode_simulation_performance)).
```@example state_projection_multi_species
using Plots # hide
using OrdinaryDiffEqRosenbrock
oprob = ODEProblem(fsp_sys, u0, 200.0, ps)
osol = solve(oprob, Rodas5P())
heatmap(osol[end]; xguide = "X₂", yguide = "X")
```
Here we perform a simulation with a long time span ($t = 200.0$) aiming to find the system's steady state distribution. Next, we plot it using the `heatmap` function.

## [Finite state projection steady state simulations](@id state_projection_steady_state_sim)
Previously, we have shown how the [SteadyStateDiffEq.jl](https://github.com/SciML/SteadyStateDiffEq.jl) package can be used to [find an ODE's steady state through forward simulation](@ref steady_state_stability). The same interface can be used for ODEs generated through FiniteStateProjection. Below, we use this to find the steady state of the dimerisation example studied in the last example.
```@example state_projection_multi_species
using SteadyStateDiffEq, OrdinaryDiffEqRosenbrock
ssprob = SteadyStateProblem(fsp_sys, u0, ps)
sssol = solve(ssprob, DynamicSS(Rodas5P()))
heatmap(sssol; xguide = "X₂", yguide = "X")
```


---
## References
[^1]: [Daniel T. Gillespie, *A rigorous derivation of the chemical master equation*, Physica A: Statistical Mechanics and its Applications (1992).](https://www.sciencedirect.com/science/article/abs/pii/037843719290283V)
[^2]: [Brian Munsky, Mustafa Khammash, *The finite state projection algorithm for the solution of the chemical master equation*, Journal of Chemical Physics (2006).](https://pubs.aip.org/aip/jcp/article-abstract/124/4/044104/561868/The-finite-state-projection-algorithm-for-the?redirectedFrom=fulltext)