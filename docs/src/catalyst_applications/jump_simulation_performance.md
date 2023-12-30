# [Advice for performant jump simulations](@id jump_simulation_performance)
We have previously described how to perform simulations of stochastic chemical kinetics *chemical reaction network* (CRN) jump process models using e.g. Gillespie's algorithm. These simulations can, however, be highly computationally intensive. Fortunately, there are several ways to increase their performance (thus reducing runtime). Here, we describe various considerations for performant stochastic chemical kinetics simulations, which we will subsequently refer to as jump process simulations. All jump process simulations arising from stochastic chemical kinetics representations of Catalyst models are performed using stochastic simulation algorithms (SSAs) from JumpProcesses.jl. Please see the [JumpProcesses documentation](https://github.com/SciML/JumpProcesses.jl) for a more extensive introduction to the package and the available solvers.

#### Brief (and optional) introduction to jump simulations
Jump processes are continuous-time, discrete-space stochastic processes. Exact realizations of these processes can be generated using Stochastic Simulation Algorithms (SSAs), of which Gillespie's Direct method is one popular choice. In the chemical reaction modeling context the discrete-state variables typically correspond to the integer-valued number of each chemical species at each time. A system's state changes at discrete time points, the jump times, when the amount of one or more species are increased by integer amounts (for example, the creation of a new protein due to translation, or the removal of one protein due to degradation). For CRNs, these jumps correspond to the occurrence of individual reactions. Typically, the frequency of each reaction depends on its *propensity* (which in turn depends on its *rate* and *substrates*). The propensity of a reaction represents its rate law, i.e. probability per time that it occurs (also known as the associated jump process' intensity function). For example, the reaction `k, A + B --> C + D` has a propensity of $k*A(t)*B(t)$ at time $t$. See [Reaction rate laws used in simulations](@ref) for more details of what propensity function Catalyst generates for a given stochastic chemical kinetics reaction.

The time of the next occurrence of some reaction, and the type of reaction that occurs, are sampled from distributions that depend on the current values of the propensity. As the latter depend on the state of the system, they must be recomputed whenever the system's state changes (for example due to a reaction occurring). Hence, jump simulations' run-times are heavily dependent on how frequently these propensities must be recomputed, and how many must be recomputed when a reaction occurs.

Typically, propensities are recomputed only whenever a jump occurs. This means that jump simulations' runtimes, very roughly, scale with the number of jumps. Runtimes typically become prohibitively expensive for:
- Simulations of large models (i.e. with many different species, where some species occur in large copy numbers, or where many reactions are present in a system).
- Simulations over long time spans.
- Simulations that are performed a large number of times.

A more thorough overview of simulation methods for Stochastic Chemical Kinetics jump process models and their computational efficiency is given in [^1].

## Managing of solution saving
By default, `solve` saves the value of the solution at the time of every jump. For simulations where a large number of jumps occur, this can cause memory to quickly fill up. Typically, for simulations with a large number of jumps, we want to [disable this feature](https://docs.sciml.ai/JumpProcesses/dev/tutorials/discrete_stochastic_example/#save_positions_docs) and instead set the save frequency manually. To exemplify this, let us consider a simple production/degradation model:
```@example jump_simulation_performance_1
using Catalyst, JumpProcesses, Plots

rn = @reaction_network begin
    (p,d), 0 <--> X
end
u0 = [:X => 10]
tspan = (0.0, 1000.0)
ps = [:p => 1000.0, :d => 100.0]

dprob = DiscreteProblem(rn, u0, tspan, ps)
nothing # hide
```
Let us simulate this model using the default options, and plot the results. Furthermore, we use [BenchmarkTools.jl's](https://github.com/JuliaCI/BenchmarkTools.jl) `@btime` macro to measure the time it takes to plot the output solution:
```@example jump_simulation_performance_1
using BenchmarkTools
jprob = JumpProblem(rn, dprob, Direct())
sol = solve(jprob, SSAStepper())
@btime plot(sol)
```
This simulation generates a very large number of jumps, with the solution saved and plotted at each time such a jump occurs. If the number of jumps is high enough, the memory limits may be exceeded (causing Julia to crash). Here, we do not reach that limit, however, the performance of the `plot(sol)` is affected by the number of points it must render (rendering too much in a single plot is another potential cause of crashed due to memory strain).

Next, we provide the `save_positions = (false, false)` option to `JumpProblem`. This turns off the saving of the solution both before and after each jump. 
```@example jump_simulation_performance_1
jprob = JumpProblem(rn, dprob, Direct(); save_positions = (false, false))
nothing # hide
```
However, if we were to simulate our model now, the solution would only actually be saved at its initial and final times. To remedy this, we provide the `saveat = 1.0` argument to the `solve` command, causing the solution to be saved at every `1.0`th time unit (at `t` = `0.0`, `1.0`, `2.0`, ... `999.0`, `1000.0`).
```@example jump_simulation_performance_1
sol = solve(jprob, SSAStepper(), saveat = 1.0)
nothing # hide
```
we can now plot the new solution:
```@example jump_simulation_performance_1
@btime plot(sol)
```
Here, we note that the time to plot the simulation was reduced (for this example, admittedly from an already low level). Furthermore, the plots look different since the trajectory is sampled at much sparser time points.

## Types of jumps
Each reaction in a chemical reaction network model corresponds to a possible jump of the jump simulation. These jumps can be divided into 3 categories:
- `MassActionJump`s: These correspond to reactions which rates are constant throughout the simulation. They are typically generated when the rate contains parameters only. 
- `ConstantRateJump`s: These correspond to reactions which rates remain constant between individual jumps, but may change in response to a jump occurring. They are typically generated when the rate contains variables and parameters only. 
- `VariableRateJump`s: These correspond to reactions which rates may change at any time during the simulation. They are typically generated when the rate depends on the time variable ($t$).

Here are some example reactions for the different types of jumps:
```@example jump_simulation_performance_1
# `MassActionJump`s
@reaction 1.0, X --> Y
@reaction k, 2X + Y --> X2Y
@reaction k1*k2+k3, 0 --> X

# `ConstantRateJump`s
@reaction k*log(X), Y --> X
@reaction mm(X,v,K), 0 --> X

# `VariableRateJump`s
@reaction k*(1+sin(t)), 0 --> X
@reaction X/t, X + Y --> XY
nothing # hide
```
Throughout a simulation, `VariableRateJump`s' rates require updating more frequently than `ConstantRateJump`s' rates, which in turn requires updating more frequently than `MassActionJump`s' rates. Since the performance of jump simulations is heavily affected by how frequently jump rates are computed, keeping the number of `ConstantRateJump`s and `VariableRateJump`s as small as possible is advantageous. Primarily, there exist two common cases where models are written in a way so that sub-optimal jump types are generated, both of which we describe below.

#### Unnecessarily putting species in rates rather than in reactions
Sometimes, a rate has been made dependent on a species, where that species instead could have been made part of the actual reaction. Consider a reaction where an enzyme ($E$) catalyses the phosphorylation of a protein ($X$) to phosphorylated form ($Xᵖ$). It can be written in two different forms:
```@example jump_simulation_performance_1
r1 = @reaction k*E, X --> Xᵖ 
r2 = @reaction k, X + E --> Xᵖ + E
nothing # hide
```
These two reactions will generate identical simulations (this holds for ODE, SDE, and Jump simulations). However, while `r1` will generate a `ConstantRateJump`, `r2` will generate a `MassActionJump`. Hence, if the `r2` form is used, jump simulation performance is, at no cost, improved. Since the two forms are otherwise identical, it is always preferable to, whenever possible, put species in the reactions, rather than the rates.

#### Using piecewise constant, time-dependant, function instead of callbacks.
Let us consider a reaction with some input that depends on time. We assume that the input is initially $0$, but after some critical time $t=5.0$ it is increased to $1$. This can be implemented through a step function:
```@example jump_simulation_performance_1
input(t) = t > 5.0
r = @reaction input(t), 0 --> X
nothing # hide
```
Here, the production of species $X$ is switched on at the time $t=5.0$. This reaction (which rate depends on `t`) will generate a `VariableRateJump`, which is bad for jump simulation performance. However, since the rate is piecewise constant, it can instead be implemented by setting it to a constant parameter $i$, and then use a [*callback*](@ref advanced_simulations_callbacks) to update it at the critical times. This will again generate an equivalent model, but with the reaction encoded as a `MassActionJump` (rather than a `VariableRateJump`).

## Jump solver selection
When creating a `JumpProblem`, a specific solver is designated using its third argument.
```@example jump_simulation_performance_2
using Catalyst, JumpProcesses, Plots

rn = @reaction_network begin
    (p,d), 0 <--> X
end
u0 = [:X => 10]
tspan = (0.0, 1000.0)
ps = [:p => 1000.0, :d => 100.0]

dprob = DiscreteProblem(rn, u0, tspan, ps)
jprob = JumpProblem(rn, dprob, Direct())
nothing # hide
```
Here (as throughout most of Catalyst's documentation) we have used the `Direct()` solver (which corresponds to Gillespie's original direct method [^2][^3], also called the *stochastic simulation algorithm*). This method was originally published in 1976, and since then, many additional methods for performing jump simulations of CRN models have been developed. 

Gillespie's direct method will, after a jump has been performed, recompute the rates of *all* possible jumps in the system. This is typically not required. E.g. consider the following system:
```@example jump_simulation_performance_2
@reaction_network begin
    k1, X1 --> X2
    k2, X2 --> X3
    k3, X3 --> 0
end
```
Here, the rate of the `k1, X1 --> X2` and `k2, X2 --> X3` reactions does not depend on the amount of $X3$ in the system. Hence, their rates are unaffected by the occurrence of the `k3, X3 --> 0` reaction. Performant jump simulation methods have clever ways to determine which rates require recomputing after the occurrence of each reaction, which improves their performance. Many of these depend on so-called dependency graphs (which track which reactions' rates are affected by the occurrence of which reactions). Catalyst automatically builds such dependency graphs, which means that most jump simulators can be used without any additional input.

A full list of jump simulation method implemented by JumpProcesses can be found [here](https://docs.sciml.ai/JumpProcesses/stable/jump_types/#Jump-Aggregators-for-Exact-Simulation). Generally, `RSSA()` (the rejection SSA method [^4][^5]) is recommended for small models, with `RSSACR()` (the rejection SSA with composition-rejection method [^6]) typically being more performant for larger models. For models that are simulated a large number of times, it can be worthwhile to try a few different jump simulation methods to determine which one is most performant in each given case.

## Hybrid simulations
For some models, copy numbers may vary greatly between different species. E.g. consider a genetic promoter which can either be in an inactive form ($Pᵢ$) or an active form ($Pₐ$). The active promoter produces a molecule ($M$):
```@example jump_simulation_performance_3
using Catalyst # hide
rn = @reaction_network begin
    (kA,kI), Pᵢ <--> Pₐ
    p, Pₐ --> Pₐ + M
    d, M --> ∅
end
```
Let us simulate this model and consider the copy numbers of each individual component:
```@example jump_simulation_performance_2
using JumpProcesses, Plots # hide
u0 = [:Pᵢ => 1, :Pₐ => 0, :M => 10000]
tspan = (0.0, 5000.0)
p = [:kA => 0.05, :kI => .01, :p => 500.0, :d => 0.0005]

dprob = DiscreteProblem(rn, u0, tspan, p)
jprob = JumpProblem(rn, dprob, RSSA())
sol = solve(jprob, SSAStepper())

plt_1 = plot(sol; idxs=2)
plt_2 = plot(sol; idxs=3)
plot(plt_1, plt_2, size=(1200,400))
```
We note that the copy numbers $Pₐ$ are highly stochastic. Meanwhile, $M$ exists in so large copy numbers that its trajectory can be considered approximately deterministic. For models like this one, jump simulations are required to capture the stochastic behaviour of the promoter. However, the large number of jumps generated by $M$ makes such simulations excessively expensive. Here, *hybrid simulations* can be used. Hybrid simulations employ different simulation strategies for their different components. In our example, the state of the promoter could be simulated using a jump simulation, while $M$ could be simulated using an ODE.

Hybrid simulations for Catalyst models are currently not supported, though they are supported by the lower-level solvers in JumpProcesses and DifferentialEquations. It is work in progress to interface Catalyst to such lower-level solvers, and to provide better optimized hybrid methods. If you require such simulations, please [raise an issue](https://github.com/SciML/Catalyst.jl/issues) and we can notify you of the current state of our implementation, or give suggestions on how one can directly interface with the lower-level hybrid solver interface.


---
## References
[^1]: [L. Marchetti, C. Priami, V. H. Thanh, *Simulation Algorithms for Computational Systems Biology*, Springer (2017).](https://link.springer.com/book/10.1007/978-3-319-63113-4)
[^2]: [D. T. Gillespie, *A general method for numerically simulating the stochastic time evolution of coupled chemical reactions*, Journal of Computational Physics (1976).](https://www.sciencedirect.com/science/article/abs/pii/0021999176900413)
[^3]: [D. T. Gillespie, *Exact Stochastic Simulation of Coupled Chemical Reactions*, The Journal of Physical Chemistry (1977).](https://pubs.acs.org/doi/10.1021/j100540a008)
[^4]: [V. H. Thanh, C. Priami and R. Zunino, *Efficient rejection-based simulation of biochemical reactions with stochastic noise and delays*, Journal of Chemical Physics (2014).](https://pubmed.ncbi.nlm.nih.gov/25296793/)
[^5]: [V. H. Thanh, R. Zunino and C. Priami, *On the rejection-based algorithm for simulation and analysis of large-scale reaction networks*, Journal of Chemical Physics (2015).](https://pubmed.ncbi.nlm.nih.gov/26133409/)
[^6]: [V. H. Thanh, R. Zunino, and C. Priami, *Efficient constant-time complexity algorithm for stochastic simulation of large reaction networks*, IEEE/ACM Transactions on Computational Biology and Bioinformatics (2017).](https://pubmed.ncbi.nlm.nih.gov/26890923/)
