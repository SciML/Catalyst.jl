# [Ensemble/Monte Carlo Simulations](@id ensemble_simulations)
In many contexts, a single model is re-simulated under similar conditions. Examples include:
- Performing Monte Carlo simulations of a stochastic model to gain insight in its behaviour.
- Scanning a model's behaviour for different parameter values and/or initial conditions.

While this can be handled using `for` loops, it is typically better to first create an `EnsembleProblem`, and then perform an ensemble simulation. Advantages include a more concise interface and the option for [automatic simulation parallelisation](@ref ode_simulation_performance_parallelisation). Here we provide a short tutorial on how to perform parallel ensemble simulations, with a more extensive documentation being available [here](https://docs.sciml.ai/DiffEqDocs/stable/features/ensemble/).

## [Monte Carlo simulations using unmodified conditions](@id ensemble_simulations_monte_carlo)
We will first consider Monte Carlo simulations where the simulation conditions are identical in-between simulations. First, we declare a [simple self-activation loop](@ref basic_CRN_library_self_activation) model
```@example ensemble
using Catalyst
sa_model = @reaction_network begin
    v0 + hill(X,v,K,n), ∅ --> X
    deg, X --> ∅
end
u0 = [:X => 10.0]
tspan = (0.0, 1000.0)
ps = [:v0 => 0.1, :v => 2.5, :K => 40.0, :n => 4.0, :deg => 0.01]
nothing # hide
```
We wish to simulate it as an SDE. Rather than performing a single simulation, however, we want to perform multiple ones. Here, we first create a normal `SDEProblem`, and use it as the single input to a `EnsembleProblem` (`EnsembleProblem` are created similarly for ODE and jump simulations, but the `ODEProblem` or `JumpProblem` is used instead).
```@example ensemble
using StochasticDiffEq
sprob = SDEProblem(sa_model, u0, tspan, ps)
eprob = EnsembleProblem(sprob)
nothing # hide
```
Next, the `EnsembleProblem` can be used as input to the `solve` command. Here, we use exactly the same inputs that we use for single simulations, however, we add a `trajectories` argument to denote how many simulations we wish to carry out. Here we perform 10 simulations:
```@example ensemble
sols = solve(eprob, STrapezoid(); trajectories = 10)
nothing # hide
```
Finally, we can use our ensemble simulation solution as input to `plot` (just like normal simulations):
```@example ensemble
using Plots
plot(sols)
```
Here, each simulation is displayed as an individual trajectory.
!!! note
    While not used here, the [`la` plotting option](@ref simulation_plotting_options) (which modifies line transparency) can help improve the plot visual when a large number of (overlapping) lines are plotted.

Various convenience functions are available for analysing and plotting ensemble simulations (a full list can be found [here]). Here, we use these to first create an `EnsembleSummary` (retrieving each simulation's value at time points `0.0, 1.0, 2.0, ... 1000.0`). Next, we use this as an input to the `plot` command, which automatically plots the mean $X$ activity across the ensemble, while also displaying the 5% and 95% quantiles as the shaded area:
```@example ensemble
e_sumary = EnsembleAnalysis.EnsembleSummary(sols, 0.0:1.0:1000.0)
plot(e_sumary)
```

## [Ensemble simulations using varying simulation conditions](@id ensemble_simulations_varying_conditions)
Previously, we assumed that each simulation used the same initial conditions and parameter values. If this is not the case (when e.g. performing a parameter scan), a `prob_func` optional argument must be supplied to the `EnsembleProblem`, this describes how the problem should be modified for each individual simulation of the ensemble.

Here, we first create an `ODEProblem` of our previous self-activation loop:
```@example ensemble
using OrdinaryDiffEqTsit5
oprob = ODEProblem(sa_model, u0, tspan, ps)
nothing # hide
```
Next, we wish to simulate the model for a range of initial conditions of $X$`. To do this we create a problem function, which takes the following arguments:
- `prob`: The problem given to our `EnsembleProblem` (which is the problem that `prob_func` modifies in each iteration).
- `i`: The number of this specific Monte Carlo iteration in the interval `1:trajectories`.
- `repeat`: The iteration of the repeat of the simulation. Typically `1`, but potentially higher if [the simulation re-running option](https://docs.sciml.ai/DiffEqDocs/stable/features/ensemble/#Building-a-Problem) is used.

Here we will use the following problem function (utilising [remake](@ref simulation_structure_interfacing_problems_remake)), which will provide a uniform range of initial concentrations of $X$:
```@example ensemble
function prob_func(prob, i, repeat)
    remake(prob; u0 = [:X => i * 5.0])
end
nothing # hide
```
Next, we create our `EnsembleProblem`, and simulate it 10 times:
```@example ensemble
eprob = EnsembleProblem(oprob; prob_func)
sols = solve(eprob, Tsit5(); trajectories = 10)
plot(sols)
```
Here we see that the deterministic model (unlike the stochastic one), only activates for some initial conditions (while other tends to an inactive state).

The `EnsembleProblem` constructor accepts a few additional optional options (`output_func`, `reduction`, `u_init`, and `safetycopy`), which are described in more detail [here](https://docs.sciml.ai/DiffEqDocs/stable/features/ensemble/#Building-a-Problem). These can be used to e.g. rerun an individual simulation which does not fulfil some condition, or extract and save only some selected information from a simulation (rather than the full simulation).
