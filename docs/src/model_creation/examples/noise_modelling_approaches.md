# [Approaches for modelling system noise](@id noise_modelling_approaches)
Catalyst's primary tools for modelling stochasticity include the creation of `SDEProblem`s or `JumpProblem`s from reaction network models. However, other approaches for incorporating model noise exist, some of which will be discussed here. We will first consider *intrinsic* and *extrinsic* noise. These are well-established terms, both of which we will describe below (however, to our knowledge, no generally agreed-upon definition of these terms exists)[^1]. Finally, we will demonstrate a third approach, the utilisation of a noisy input process to an otherwise deterministic system. This approach is infrequently used, however, as it is encountered in the literature, we will demonstrate it here as well. 

We note that these approaches can all be combined. E.g. an intrinsic noise model (using an SDE) can be combined with extrinsic noise (using randomised parameter values), while also feeding a noisy input process into the system.

!!! note
    Here we use intrinsic and extrinsic noise as descriptions of two of our modelling approaches. It should be noted that while these are established terminologies for noisy biological systems[^1], our use of these terms to describe different approaches for modelling noise is only inspired by this terminology, and nothing that is established in the field.  Please consider the [references](@ref noise_modelling_approaches_references) for more information on intrinsic and extrinsic noise.

## [The repressilator model](@id noise_modelling_approaches_model_intro)
For this tutorial we will use the oscillating [repressilator](@ref basic_CRN_library_repressilator) model.
```@example noise_modelling_approaches
using Catalyst
repressilator = @reaction_network begin
    hillr(Z,v,K,n), ∅ --> X
    hillr(X,v,K,n), ∅ --> Y
    hillr(Y,v,K,n), ∅ --> Z
    d, (X, Y, Z) --> ∅
end
```

## [Using intrinsic noise](@id noise_modelling_approaches_model_intrinsic)
Generally, intrinsic noise is randomness inherent to a system itself. This means that it cannot be controlled for, or filtered out by, experimental settings. Low-copy number cellular systems, were reaction occurs due to the encounters of molecules due to random diffusion, is an example of intrinsic noise. In practise, this can be modelled exactly through [SDE](@ref simulation_intro_SDEs) (chemical Langevin equations) or [jump](@ref simulation_intro_jumps) (stochastic chemical kinetics) simulations. 

In Catalyst, intrinsic noise is accounted for whenever an `SDEProblem` or `JumpProblem` is created and simulated. Here we will model intrinsic noise through SDEs, which means creating an `SDEProblem` using the standard approach.
```@example noise_modelling_approaches
u0 = [:X => 45.0, :Y => 20.0, :Z => 20.0]
tend = 200.0
ps = [:v => 10.0, :K => 20.0, :n => 3, :d => 0.1]
sprob = SDEProblem(repressilator, u0, tend, ps)
nothing # hide
```
Next, to illustrate the system's noisiness, we will perform multiple simulations. We do this by [creating an `EnsembleProblem`](@ref ensemble_simulations_monte_carlo). From it, we perform, and plot, 4 simulations.
```@example noise_modelling_approaches
using StochasticDiffEq, Plots
eprob_intrinsic = EnsembleProblem(sprob)
sol_intrinsic = solve(eprob_intrinsic, ImplicitEM(); trajectories = 4)
plot(sol_intrinsic; idxs = :X)
```
Here, each simulation is performed from the same system using the same settings. Despite this, due to the noise, the individual trajectories are different.

## [Using extrinsic noise](@id noise_modelling_approaches_model_extrinsic)
Next, we consider extrinsic noise. This is randomness caused by stochasticity external to, yet affecting, a system. Examples could be different bacteria experiencing different microenvironments or cells being in different parts of the cell cycle. This is noise which (in theory) can be controlled for experimentally (e.g. by ensuring a uniform environment). Whenever a specific source of noise is intrinsic and extrinsic to a system may depend on how one defines the system itself (this is a reason why giving an exact definition of these terms is difficult).

In Catalyst we can model extrinsic noise by letting the model parameters be probability distributions. Here, at the beginning of each simulation, random parameter values are drawn from their distributions. Let us imagine that our repressilator circuit was inserted into a bacterial population. Here, while each bacteria would have the same circuit, their individual number of e.g. ribosomes (which will be random) might affect the production rates (which while constant within each bacteria, might differ between the individuals).

Again we will perform ensemble simulation. Instead of creating an `SDEProblem`, we will create an `ODEProblem`, as well as a [problem function](@ref ensemble_simulations_varying_conditions) which draws random parameter values for each simulation. Here we have implemented the parameter's probability distributions as [normal distributions](https://en.wikipedia.org/wiki/Normal_distribution) using the [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) package.
```@example noise_modelling_approaches
using Distributions
p_dists = Dict([:v => Normal(10.0, 2.0), :K => Normal(20.0, 5.0), :n => Normal(3, 0.2), :d => Normal(0.1, 0.02)])
function prob_func(prob, i, repeat)
    p = [par => rand(p_dists[par]) for par in keys(p_dists)]
    return remake(prob; p)
end
nothing # hide
```
Next, we again perform 4 simulations. While the individual trajectories are performed using deterministic simulations, the randomised parameter values create heterogeneity across the ensemble.
```@example noise_modelling_approaches
using OrdinaryDiffEqDefault
oprob = ODEProblem(repressilator, u0, tend, ps)
eprob_extrinsic = EnsembleProblem(oprob; prob_func)
sol_extrinsic = solve(eprob_extrinsic; trajectories = 4)
plot(sol_extrinsic; idxs = :X)
```
We note that a similar approach can be used to also randomise the initial conditions. In a very detailed model, the parameter values could fluctuate during a single simulation, something which could be implemented using the approach from the next section.

## [Using a noisy input process](@id noise_modelling_approaches_model_input_noise)
Finally, we will consider the case where we have a deterministic system, but which is exposed to a noisy input process. One example could be a [light sensitive system, where the amount of experienced sunlight is stochastic due to e.g. variable cloud cover](@ref functional_parameters_circ_rhythm). Practically, this can be considered as extrinsic noise, however, we will generate the noise using a different approach from in the previous section. Here, we pre-simulate a random process in time, which we then feed into the system as a functional, time-dependent, parameter. A more detailed description of functional parameters can be found [here](@ref time_dependent_parameters).

We assume that our repressilator has an input, which corresponds to the $K$ value that controls $X$'s production. First we create a function, `make_K_series`, which creates a randomised time series representing $K$'s value over time. 
```@example noise_modelling_approaches
using DataInterpolations
function make_K_series(; K_mean = 20.0, n = 500, θ = 0.01)
    t_samples = range(0.0, stop = tend, length = n)
    K_series = fill(K_mean, n)
    for i = 2:n
        K_series[i] = K_series[i-1] + (rand() - 0.5) - θ*(K_series[i-1] - K_mean)
    end
    return LinearInterpolation(K_series, t_samples)
end
plot(make_K_series())
```
Next, we create an updated repressilator model, where the input $K$ value is modelled as a time-dependent parameter.
```@example noise_modelling_approaches
@parameters (K_in::typeof(make_K_series()))(..)
K_in = K_in(default_t())
repressilator_Kin = @reaction_network begin
    hillr(Z,v,$K_in,n), ∅ --> X
    hillr(X,v,K,n), ∅ --> Y
    hillr(Y,v,K,n), ∅ --> Z
    d, (X, Y, Z) --> ∅
end
```
Finally, we will again perform ensemble simulations of our model. This time, at the beginning of each simulation, we will use `make_K_series` to generate a new $K$, and set this as the `K_in` parameter's value.
```@example noise_modelling_approaches
function prob_func_Kin(prob, i, repeat)
    p = [ps; :K_in => make_K_series()]    
    return ODEProblem(repressilator_Kin, prob.u0, prob.tspan, p)
end
oprob = ODEProblem(repressilator_Kin, u0, tend, [ps; :K_in => make_K_series()])
eprob_inputnoise = EnsembleProblem(oprob; prob_func = prob_func_Kin)
sol_inputnoise = solve(eprob_inputnoise; trajectories = 4)
plot(sol_inputnoise; idxs = :X)
```
Like in the previous two cases, this generates heterogeneous trajectories across our ensemble.


## [Investigating the mean of noisy oscillations](@id noise_modelling_approaches_model_noisy_oscillation_mean)
Finally, we will observe an interesting phenomenon for ensembles of stochastic oscillators. First, we create ensemble simulations with a larger number of trajectories.
```@example noise_modelling_approaches
sol_intrinsic = solve(eprob_intrinsic, ImplicitEM(); trajectories = 200)
sol_extrinsic = solve(eprob_extrinsic; trajectories = 200)
```
Next, we can use the `EnsembleSummary` interface to plot each ensemble's mean activity (as well as 5% and 95% quantiles) over time:
```@example noise_modelling_approaches
e_sumary_intrinsic = EnsembleAnalysis.EnsembleSummary(sols, 0.0:1.0:tend)
e_sumary_extrinsic = EnsembleAnalysis.EnsembleSummary(sols, 0.0:1.0:tend)
plot(e_sumary_intrinsic; label = "Intrinsic noise", idxs = 1)
plot!(e_sumary_extrinsic; label = "Extrinsic noise", idxs = 1)
```
Here we can see that, over time, the systems' mean $X$ activity reaches a constant level around $30$.

This is a well-known phenomenon (especially in circadian biology[^2]). Here, as stochastic oscillators evolve from a common initial condition the mean behaves as a damped oscillator. This can be caused by two different phenomena:
- The individual trajectories are themselves damped.
- The individual trajectories's phases get de-synchronised.
However, if we only observe the mean behaviour (and not the single trajectories), we cannot know which of these cases we are encountering. Here, by checking the single-trajectory plots from the previous sections, we note that this is due to trajectory de-synchronisation. Stochastic oscillators have often been cited as a reason for the importance to study cellular systems at the *single-cell* level, and not just in bulk.

---
## [References](@id noise_modelling_approaches_references)
[^1]: [Michael B. Elowitz, Arnold J. Levine, Eric D. Siggia, Peter S. Swain, *Stochastic Gene Expression in a Single Cell*, Science (2002).](https://www.science.org/doi/10.1126/science.1070919)
[^2]: [Qiong Yang, Bernardo F. Pando, Guogang Dong, Susan S. Golden, Alexander van Oudenaarden, *Circadian Gating of the Cell Cycle Revealed in Single Cyanobacterial Cells*, Science (2010).](https://www.science.org/doi/10.1126/science.1181759)