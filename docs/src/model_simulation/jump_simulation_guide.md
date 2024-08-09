# [Advice for Stochastic Chemical Kinetics Jump Simulations](@id jump_simulation_guide)
We have previously introduced how Catalyst *chemical reaction network* (CRN)
models can be converted to *stochastic chemical kinetics* jump process models,
which can then be (exactly) simulated (i.e. sampled) using Stochastic Simulation
Algorithms (SSAs) such as Gillespie's Direct method.

In this tutorial we will review basic information and configuration options for
such simulations, and discuss several ways to increase their performance (thus
reducing run time). All jump simulations arising from stochastic chemical
kinetics representations of Catalyst models are performed using SSAs from
JumpProcesses.jl. Please see the [JumpProcesses
documentation](https://github.com/SciML/JumpProcesses.jl) for a more extensive
introduction to the package and its available solvers.

### [Brief (and optional) introduction to stochastic chemical kinetics jump simulations](@id jump_simulation_guide_intro)
Jump processes are continuous-time, discrete-space, stochastic processes. Exact
realizations of these processes can be generated using SSAs (of which
Gillespie's Direct method is the most well-known choice). In the chemical
reaction modelling context, the discrete-state variables typically correspond to
the integer-valued number of each chemical species at each time. A system's
state changes at discrete time points corresponding to when reactions occur (the
jump times). A these these times the amount of one (or more) species are changed
by integer amount(s) (for example, the creation of a new protein due to
translation, or the removal of one protein due to degradation).

The frequency of each reaction's occurrence depends on its *propensity* function
(which in turn depends on its *rate constant* and *substrate* amounts). The
propensity of a reaction is analogous to the reaction's *rate law* in the ODE
context, and represents the probability per unit time the reaction can occur
given the current state of the system. In probability and statistics
propensities are also often called intensity functions or transition rate
functions. For example, the reaction $A + A \overset{k}{\to} B$ has a rate
constant of $k$ and a propensity of $k A (A-1) / 2$, while the reaction $A + B
\overset{\gamma}{\to} C + D$ has a rate constant of $\gamma$ and a propensity of
$k A B$. See [Reaction rate laws used in simulations](@ref
introduction_to_catalyst_ratelaws) for more details of what propensity functions
Catalyst generates for a given reaction, and [Mathematical Models Catalyst can
Generate](@ref math_models_in_catalyst) for details on both propensity functions
and the mathematical jump process models Catalyst can generate from a CRN (called the time-change representation in the stochastic chemical kinetics literature).

During a typical simulation of a jump process model, the simulation algorithm
samples both the time to the next reaction and which reaction will occur at this
time point. The sampled time determines the timestep within a simulation, while
the sampled reaction determines how the state is updated after stepping to the
new time. Both are computed by sampling probability distributions that depend on
the current values for the reaction propensities. More precisely, in the
notation of the [Mathematical Models Catalyst can Generate](@ref
math_models_in_catalyst) guide, when the current time is $t$, the current system
state is $X(t)$ (i.e. the vector of the number of each chemical species), and
the current propensity values are $\{a_k(X(t))\}_{k=1}^K$, a typical timestep of
an SSA
1. Samples both a *random* time $\Delta t$ until the next reaction occurs and
   which reaction, $k \in \{1,\dots,K\}$, occurs at this time. This sampling
   process uses probability distributions that are built from the current
   propensity values, i.e. $\{a_k(X(t))\}_{k=1}^K$.
2. Executes the $k$'th reaction calculating $X(t + \Delta t)$ from $X(t)$.
3. Recalculates one or more propensity functions, $a_k(X(t + \Delta
   t))$, based on the new state.
4. Updates the time $t \to t + \Delta t$.

Hence, the average work per one time step of a jump simulation is heavily
dependent on how many propensities must typically be recomputed when a reaction
occurs. Similarly, the total work across a full jump simulation will also depend
on typical values for the timestep. Propensity function values generally
increase as populations increase (i.e. for $A + B \overset{k}{\to} C$ the
propensity is $k A B$ which increases in the amount of $A$ or $B$). Larger
propensities mean that reactions occur more frequently, and hence typical values
for $\Delta t$ will decrease and more timesteps will be taken when simulating to
a fixed (physical) time, $T$.

A variety of factors can therefore impact the run time of different SSAs. Such
simulations can become increasingly expensive for:
- Simulations of large models (i.e. with many different species, where some
  species occur in large copy numbers, or where many reactions are present in a
  system).
- Simulations over long time spans or with large numbers of reactions occurring
  over the desired time span.
- Simulations that are performed a large number of times for statistical sampling.
- Simulations that are saving large amounts of data (such as saving the values
  of all species each time a reaction occurs).

A more thorough overview of simulation methods for Stochastic Chemical Kinetics
jump process models and their computational efficiency is given in [^1].

## [Basic simulation behavior]


## [Managing of solution saving](@id jump_simulation_guide_solution_saving)
By default, `solve` saves the value of the solution at the time of every jump. For simulations with a large number of jump events, this can cause memory to quickly fill up. Typically, for simulations with a large number of jumps, we want to [disable this feature](https://docs.sciml.ai/JumpProcesses/dev/tutorials/discrete_stochastic_example/#save_positions_docs) and instead set the save frequency manually. Let us consider a simple [birth-death model](@ref basic_CRN_library_bd):
```@example jump_simulation_guide_1
using Catalyst, JumpProcesses, Plots

bd_model = @reaction_network begin
    (p,d), 0 <--> X
end
u0 = [:X => 10]
tspan = (0.0, 1000.0)
ps = [:p => 1000.0, :d => 100.0]

dprob = DiscreteProblem(bd_model, u0, tspan, ps)
nothing # hide
```
Let us simulate it using the default options, and plot the results. Furthermore, we use [BenchmarkTools.jl's](https://github.com/JuliaCI/BenchmarkTools.jl) `@btime` macro to measure the time it takes to plot the output solution:
```@example jump_simulation_guide_1
using BenchmarkTools
jprob = JumpProblem(bd_model, dprob, Direct())
sol = solve(jprob, SSAStepper())
sol = solve(jprob, SSAStepper(); seed = 1234) # hide
@btime plot(sol)
```
This simulation generates a very large number of jumps, with the solution saved and plotted at each time such a jump occurs. If the number of jumps is high enough, the memory limits may be exceeded (causing Julia to crash). Here, we do not reach this limit, however, the run time of the `plot(sol)` is affected by the number of points it must render (rendering too much in a single plot is another potential cause of crashed due to memory strain).

Next, we provide the `save_positions = (false, false)` option to `JumpProblem`. This turns off the saving of the solution both before and after each jump.
```@example jump_simulation_guide_1
jprob = JumpProblem(bd_model, dprob, Direct(); save_positions = (false, false))
nothing # hide
```
However, if we were to simulate our model now, the solution would only actually be saved at its initial and final times. To remedy this, we provide the `saveat = 1.0` argument to the `solve` command, ensuring that the solution to be saved at every `1.0`th time unit (that is, at `t` = `0.0`, `1.0`, `2.0`, ... `999.0`, `1000.0`).
```@example jump_simulation_guide_1
sol = solve(jprob, SSAStepper(); saveat = 1.0)
sol = solve(jprob, SSAStepper(); saveat = 1.0, seed = 1234) # hide
nothing # hide
```
we can now plot the new solution:
```@example jump_simulation_guide_1
@btime plot(sol)
```
Here, we note that the time to plot the simulation was reduced (for this example, admittedly from an already low level). Furthermore, the plots look different since the trajectory is sampled at much sparser time points.

!!! note
    With the default saving behaviour, [evaluating `sol(t)`](@ref simulation_structure_interfacing_solutions) at any time `t` within the `tspan` gives the *exact* value of the solution at that time (i.e. the vector of the exact number of each species). When using `saveat` and `save_positions = (false,false)` the solution is only saved at the selected time points, and as such `sol(t)` should not be evaluated except at times at which the solution was saved. While evaluating the solution at other times will return values they will not be the exact value of the jump process!

## [Types of jumps](@id jump_simulation_guide_jump_types)
Each reaction in a chemical reaction network model corresponds to a possible jump of the jump simulation. These jumps can be divided into 3 categories:
- `MassActionJump`s: These correspond to reactions which rates remain constant throughout the simulation.They are typically generated by reactions' which rates contain time-independent parameters only.
- `ConstantRateJump`s: These correspond to reactions which rates remain constant between individual jumps, but may change in response to a jump occurring. They are typically generated by reactions' which rates contain species and time-independent parameters only.
- `VariableRateJump`s: These correspond to reactions which rates may change at any time during the simulation. They are typically generated by reactions' which rates contain the time variable ($t$).

Here are some example reactions for the different types of jumps:
```@example jump_simulation_guide_1
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
Updating `MassActionJump`s' propensities is more computationally efficient (due to their constrained form) than updating them for `ConstantRateJump`s. Thus simulations are more performant if a larger fraction of reactions are of the `MassActionJump` form (rather than the `ConstantRateJump` form). Furthermore, simulations containing `VariableRateJump` requires additional routines to find the times between jump events (as unlike non-`VariableRateJump` simulations, jump propensities change *between* jump events, making this computation more difficult). Hence, the existence of `VariableRateJump` can significantly slow down a simulation Primarily, there exist two common situation where models are written in a way so that sub-optimal jump types are generated, both of which we describe below.

### [Unnecessarily putting species in rates rather than in reactions](@id jump_simulation_guide_jump_types_unnecessary_constantratejumps)
Sometimes, a rate has been made dependent on a species, where that species instead could have been made part of the actual reaction. Consider a reaction where an enzyme ($E$) catalyses the phosphorylation of a protein ($X$) to phosphorylated form ($Xᵖ$). It can be written in two different forms:
```@example jump_simulation_guide_1
r1 = @reaction k*E, X --> Xᵖ
r2 = @reaction k, X + E --> Xᵖ + E
nothing # hide
```
These two reactions will generate identical simulations (this holds for ODE, SDE, and Jump simulations). However, while `r1` will generate a `ConstantRateJump`, `r2` will generate a `MassActionJump`. Hence, if the `r2` form is used, jump simulation performance is (at no cost) improved. Since the two forms are otherwise identical, it is always preferable to, whenever possible, put species in the reactions, rather than the rates.

### [Using piecewise constant, time-dependant, function instead of events](@id jump_simulation_guide_jump_types_unnecessary_variableratejumps)
Let us consider a reaction with some input that depends on time. We assume that the input is initially $0$, but after some critical time $t = 5.0$ it is increased to $1$. This can be implemented through a step function:
```@example jump_simulation_guide_1
input(t) = t > 5.0
r = @reaction input(t), 0 --> X
nothing # hide
```
Here, the production of species $X$ is switched on at the time $t = 5.0$. This reaction (for which the rate depends on `t`) will generate a `VariableRateJump`, which is the least performant jump type in simulations. However, since the rate is piecewise constant, it can alternatively be implemented by setting it to a constant parameter $i$, and then using a [*discrete event*](@ref constraint_equations_events) to update it at the switching times. This will again generate an equivalent model, but with the reaction encoded as a `MassActionJump` (rather than a `VariableRateJump`). Generally such explicit time-discontinuities should be encoded via discrete callbacks instead of as `VariableRateJump`s if possible (as simulation methods for the latter typically assume the system's propensities evolve continuously in-between jumps).

## [Jump solver selection](@id jump_simulation_guide_solver_selection)
When creating a `JumpProblem`, a specific solver is designated using its third argument.
```@example jump_simulation_guide_2
using Catalyst, JumpProcesses, Plots

bd_model = @reaction_network begin
    (p,d), 0 <--> X
end
u0 = [:X => 10]
tspan = (0.0, 1000.0)
ps = [:p => 1000.0, :d => 100.0]

dprob = DiscreteProblem(bd_model, u0, tspan, ps)
jprob = JumpProblem(bd_model, dprob, Direct())
nothing # hide
```
Here (as throughout most of Catalyst's documentation) we have used the `Direct()` SSA solver (which corresponds to Gillespie's original direct method [^2][^3]). This method was originally published in 1976, and since then, many additional methods for simulating stochastic chemical kinetics models have been developed.

Gillespie's direct method will, after a jump has been performed, recompute the propensities of *all* possible jumps in the system (i.e. of all reactions). This is typically not required. E.g. consider the following system:
```@example jump_simulation_guide_2
@reaction_network begin
    k1, X1 --> X2
    k2, X2 --> X3
    k3, X3 --> 0
end
```
Here, the propensities of the `k1, X1 --> X2` and `k2, X2 --> X3` reactions (`k1*X1*X2` and `k2*X2` respectively) do not depend on the amount of $X3$ in the system. Hence, the value of their propensities are unchanged by the occurrence of the `k3, X3 --> 0` reaction. Performant jump simulation methods have several ways to determine which rates require recomputing after the occurrence of a reaction, which improves their computational performance. Many of these requires a dependency graphs, which track which other reactions' propensities must be recomputed after the occurrence of a given reaction. Catalyst automatically builds such dependency graphs, which means that all JumpProcesses SSAs can be used without any additional inputs.

A full list of jump simulation method implemented by JumpProcesses can be found [here](https://docs.sciml.ai/JumpProcesses/stable/jump_types/#Jump-Aggregators-for-Exact-Simulation). Generally, `RSSA()` (the rejection SSA method [^4][^5]) is recommended for small models, with `RSSACR()` (the rejection SSA with composition-rejection method [^6]) typically being more performant for larger models. For models that are simulated a large number of times, it can be worthwhile to try a few different jump simulation methods to determine which one is most performant in each given case.

## [Hybrid simulations](@id jump_simulation_guide_hybrid_simulations)
For some models, copy numbers may vary greatly between different species. E.g. consider a genetic promoter which can either be in an inactive form ($Pᵢ$) or an active form ($Pₐ$). The active promoter produces a molecule ($M$):
```@example jump_simulation_guide_3
using Catalyst # hide
promoter = @reaction_network begin
    (kA,kI), Pᵢ <--> Pₐ
    p, Pₐ --> Pₐ + M
    d, M --> ∅
end
```
Let us simulate this model and consider the copy numbers of each individual component:
```@example jump_simulation_guide_3
using JumpProcesses, Plots # hide
u0 = [:Pᵢ => 1, :Pₐ => 0, :M => 10000]
tspan = (0.0, 5000.0)
p = [:kA => 0.05, :kI => .01, :p => 500.0, :d => 0.0005]

dprob = DiscreteProblem(promoter, u0, tspan, p)
jprob = JumpProblem(promoter, dprob, RSSA())
sol = solve(jprob, SSAStepper())
sol = solve(jprob, SSAStepper(); seed = 1234) # hide

plt_1 = plot(sol; idxs=2)
plt_2 = plot(sol; idxs=3)
plot(plt_1, plt_2, size=(1200, 400))
```
We note that the copy numbers $Pₐ$ are highly stochastic. Meanwhile, $M$ exists in so large copy numbers that its trajectory can be considered approximately deterministic. For models like this one, jump simulations are required to capture the stochastic behaviour of the promoter. However, the large number of jumps generated by $M$ makes such simulations excessively expensive. Here, *hybrid simulations* can be used. Hybrid simulations employ different simulation strategies for their different components. In our example, the state of the promoter could be simulated using a jump simulation, while $M$ could be simulated using an ODE.

Hybrid simulations for Catalyst models are currently not supported, though they are supported by the lower-level solvers in JumpProcesses and OrdinaryDiffEq. It is work in progress to interface Catalyst to such lower-level solvers, and to provide better optimised hybrid methods. If you require such simulations, please [raise an issue](https://github.com/SciML/Catalyst.jl/issues) and we can notify you of the current state of our implementation, or give suggestions on how one can directly interface with the lower-level hybrid solver interface.

## [Simulation parallelisation](@id jump_simulation_guide_parallelisation)
When multiple simulations are carried out, simulations can be improved by running these in parallel. We have [previously described](@ref ode_simulation_performance_parallelisation) how to do so for ODE simulation. Jump simulations can be parallelised in the same manner, with the exception that GPU parallelisation is currently not supported.

---
## References
[^1]: [L. Marchetti, C. Priami, V. H. Thanh, *Simulation Algorithms for Computational Systems Biology*, Springer (2017).](https://link.springer.com/book/10.1007/978-3-319-63113-4)
[^2]: [D. T. Gillespie, *A general method for numerically simulating the stochastic time evolution of coupled chemical reactions*, Journal of Computational Physics (1976).](https://www.sciencedirect.com/science/article/abs/pii/0021999176900413)
[^3]: [D. T. Gillespie, *Exact Stochastic Simulation of Coupled Chemical Reactions*, The Journal of Physical Chemistry (1977).](https://pubs.acs.org/doi/10.1021/j100540a008)
[^4]: [V. H. Thanh, C. Priami and R. Zunino, *Efficient rejection-based simulation of biochemical reactions with stochastic noise and delays*, Journal of Chemical Physics (2014).](https://pubmed.ncbi.nlm.nih.gov/25296793/)
[^5]: [V. H. Thanh, R. Zunino and C. Priami, *On the rejection-based algorithm for simulation and analysis of large-scale reaction networks*, Journal of Chemical Physics (2015).](https://pubmed.ncbi.nlm.nih.gov/26133409/)
[^6]: [V. H. Thanh, R. Zunino, and C. Priami, *Efficient constant-time complexity algorithm for stochastic simulation of large reaction networks*, IEEE/ACM Transactions on Computational Biology and Bioinformatics (2017).](https://pubmed.ncbi.nlm.nih.gov/26890923/)
