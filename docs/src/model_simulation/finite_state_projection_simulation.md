# [Solving the chemical master equation using FiniteStateProjection.jl](@id finite-state_projection)
Previously we have shown how *stochastic chemical kinetics* describe how chemical reaction network models can be [exactly simulated](@ref simulation_intro_jumps) (using e.g. Gillespie's algorithm). We also described how the [SDE](@ref simulation_intro_SDEs) and [ODE](@ref simulation_intro_ODEs) approaches were approximation of these jump simulations, and only valid for large copy numbers. To gain a good understanding of the system's time-development, we typically have to carry out a large number of jump simulations. Here, an alternative approach is instead to simulate the *full probability distribution of the system*. This corresponds to the distribution from which these jump simulations are drawn.

[*The chemical master equation*](https://en.wikipedia.org/wiki/Master_equation) (CME) describes the time development of this distribution[^1]. In fact, this equation is at the core of chemical reaction network kinetics, with all other approaches (such as ODE, SDE, and Jump simulations) being derived as various approximations of it. The CME is a system of ODEs, with one variable *for each possible state of the system*. Each of these variables describes the probability of the system being in that state, over time. For a system with a single species $X$, the CME looks like
```math
\begin{aligned}
\frac{dp(x=0)}{dt} &= f_0(p(x=0), p(x=1), ...) \\
\frac{dp(x=1)}{dt} &= f_0(p(x=0), p(x=1), ...) \\
\frac{dp(x=2)}{dt} &= f_0(p(x=0), p(x=1), ...) \\
                &\vdots\\
\end{aligned}
```
Here, we note that, since (for almost all chemical reaction networks) there is a non-zero probability that $X$ is any specific integer value, the CME have an infinite number of variables. Hence, it cannot be solved practically. However, one notes that for high enough species values, the probability of the system attaining such values becomes negligibly small. Here, a truncated version of the CME can be solved practically. An approach for this is the *finite state projection*[^2]. Below we describe how to use the [FiniteStateProjection.jl](https://github.com/SciML/FiniteStateProjection.jl) package to solve the truncated CME. While this approach can be very powerful, we note that for system's with many species, even the truncated CME typically have too many states to be feasible to solve.

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
bar(osol(1.0);  bar_width = 1.0, linewidth = 0, alpha = 0.7, label = "t = 1.0")
bar!(osol(2.0); bar_width = 1.0, linewidth = 0, alpha = 0.7, label = "t = 2.0")
bar!(osol(5.0); bar_width = 1.0, linewidth = 0, alpha = 0.7, label = "t = 5.0")
bar!(osol(10.0); bar_width = 1.0, linewidth = 0, alpha = 0.7, label = "t = 10.0")
```

## [Finite state projection simulation of multi-species model](@id state_projection_multi_species)
Next, we will consider a system with more than one species. The workflow will be identical, however, we will have to make an additional consideration regarding our initial conditions. We will also need to use a different plotting approach.

For this example we will consider a simple dimerisation model. In it, $X$ gets produced and degraded at constant rates, and can also dimerise to form $X₂$.
```@example state_projection_multi_species
using Catalyst # hide
rs = @reaction_network begin
    (p,d), 0 <--> X
    (kB,kD), 2X <--> X₂
end
```

Next, we will declare our parameter values and initial condition. In this case, the initial condition is a matrix where element $(i,j)$ denotes the initial probability that $(X(0),X₂(0)) = (i-1,j-1)$. In this case, we will use an initial condition where we know that $(X(0),X₂(0)) = (0,0)$.
```@example state_projection_multi_species
ps = [:p => 1.0, :d => 0.2, :kB => 2.0, :kD => 5.0]
u0 = zeros(25,25)
u0[1,1] = 1.0
```
In the next step, however, we have to make an additional consideration. Since we have more than one species, we have to define which dimension of the initial condition (and hence also the output solution) correspond to which species. Here we provide a second argument to `FSPSystem`, which is a vector listing all species in the order they occur in the `u0` array.
```@example state_projection_multi_species
using FiniteStateProjection # hide
fsp_sys = FSPSystem(rs, [:X, :X₂])
nothing # hide
```
Finally, we can simulate the model just like in the 1-dimensional case. Here we will specifically use the `Rodas5P()` ODE solver (as it performs noticeably better than teh default choice).
```@example state_projection_multi_species
using Plots # hide
using OrdinaryDiffEqRosenbrock
oprob = ODEProblem(fsp_sys, u0, 100.0, ps)
@time osol = solve(oprob, Rodas5P())
heatmap(osol[end]; xguide = "X₂", yguide = "X")
```
Here we perform a simulation with a long time span ($t = 100$) aiming to find the system's steady state distribution. Next, we plot it using the `heatmap` function.

## [Finite state projection steady state simulations](@id state_projection_steady_state_sim)
Previously we described how the [SteadyStateDiffEq.jl](https://github.com/SciML/SteadyStateDiffEq.jl) package can be used to [find an ODE's steady state through forward simulation](@ref steady_state_stability). The same interface can be used for ODEs generated through FiniteStateProjection.jl. Below we use this to find the steady state of the dimerisation example studied in the last example.
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