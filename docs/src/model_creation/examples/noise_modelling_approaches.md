# [Approaches for modelling system noise](@id noise_modelling_approaches)
Catalyst's primary tools for modelling stochasticity includes the creation of `SDEProblem`s or `JumpProblem`s from reaction network models. However, other approaches for incorporating noise in a model exists, some of which will be discussed here. We will first consider *intrinsic* and *extrinsic* noise. These are well-established terms, both of which we will describe in detail below (however, to our knowledge, no generally agreed upon definition of these exists)[^1]. The final approach, the utilisation of a noisy input process to an otherwise deterministic system, is less frequently used. However, as it has been used in the literature, we will demonstrate it here as well. 

Finally, we note that these approaches can all be combined. E.g. modelling intrinsic noise (using an SDE) can be combined extrinsic noise (using randomised parameter values) while also feeding a noisy input process into the system.

!!! note
    Here we use intrinsic and extrinsic noise as description of two of our modelling approaches. It should be noted that while these are established terminologies for noisy biological systems[^1], our use of these terms to describe different approaches for modelling noise is only inspired by this terminology, and nothing that is established in the field.

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
Generally, intrinsic noise is randomness inherent to a system itself. This means that it cannot be controlled for, or filtered out by, experimental settings. Low-copy number cellular systems, were reaction occurs due to the encounters of molecules due to random diffusion, is an example of intrinsic noise. In practise, this can be precisely modelled through [SDE](@ref simulation_intro_SDEs) (chemical Langevin equations) or [jump](@ref simulation_intro_jumps) (stochastic chemical kinetics) simulations. 

In Catalyst, intrinsic noise is accounted for whenever a `SDEProblem` or `JumpProblem` is created and simulated. Here we will model intrinsic noise through SDEs, which means creating a `SDEProblem` using the standard approach.
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
Next, we consider extrinsic noise. This is randomness caused by stochasticity external to, yet affecting, a system. Examples could be different bacteria experiencing different microenvironments or cells being in different part of the cell cycle. Here, what noise is intrinsic and extrinsic to a system may depend on how one defines the system itself (which is a reason why exact definitions of these terms is difficult, please consider the [references](@ref noise_modelling_approaches_references) for more information).

In Catalyst we can modelled extrinsic noise by letting the model parameters be probability distributions. Here, at the beginning of each simulation, random parameter values are drawn from their distributions. Let us imagine that our repressilator circuit was inserted into an bacterial population. Here, while each bacteria would have the same circuit, their individual number of e.g. ribosomes (which will be random) might affect the production rates (which while constant within each bacteria, might differ between the individuals).

Again we will perform ensemble simulation. Instead of creating an `SDEProblem`, we will create an `ODEProblem`, as well as a [problem function](@ref ensemble_simulations_varying_conditions) which draws random parameter values for each simulations (here we have implemented the parameter's probability distributions using the [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) package).
```@example noise_modelling_approaches
using Distributions
p_dists = Dict([:v => Normal(10.0, 2.0), :K => Normal(20.0, 5.0), 
           :n => Normal(3, 0.2), :d => Normal(0.1, 0.02)])
function prob_func(prob, i, repeat)
    p = [par => rand(p_dists[par]) for par in keys(p_dists)]
    return remake(prob; p)
end
nothing # hide
```
Next, we again performs 4 simulations. While the individual trajectories are performed using deterministic simulations, the randomised parameter values creates heterogeneity across the ensemble.
```@example noise_modelling_approaches
using OrdinaryDiffEqDefault
oprob = ODEProblem(repressilator, u0, tend, ps)
eprob_extrinsic = EnsembleProblem(oprob; prob_func)
sol_extrinsic = solve(eprob_extrinsic; trajectories = 4)
plot(sol_extrinsic; idxs = :X)
```
We note that a similar approach can be used to also randomise the initial conditions. Finally, in a very detailed model the parameter values could fluctuate across the simulations, something which could be implemented using the approach from the next section.

## [Using a noisy input process](@id noise_modelling_approaches_model_input_noise)
Finally, we will consider the case where we have a deterministic system, but which is exposed to a some noisy input process. One example could be a [light sensitive system, where the amount of experienced sunlight is stochastic due to e.g. variable cloud cover](@ref functional_parameters_circ_rhythm). Practically, this can be considered extrinsic noise, but the model is created using a different approach from teh previous section. Here, we will pre-simulate a random process in time, which we then feed into the system as a function, time-dependent, parameter. The technical aspects of this approach is described in more details [here](@ref time_dependent_parameters).

He we assume that our repressilator have an input, which correspond to the $K$ value controlling $X$'s production. First we create a function, `make_K_series`, which creates a randomised time series representing $K$'s value over time. 
```@example noise_modelling_approaches
using DataInterpolations
function make_K_series(;K_mean = 20.0, n = 500, θ = 0.01)
    t_samples = range(0.0, stop = tend, length = n)
    K_series = fill(K_mean, n)
    for i = 2:n
        K_series[i] = K_series[i-1] + (rand()-0.5) - θ*(K_series[i-1] - K_mean)
    end
    return LinearInterpolation(K_series, t_samples)
end
plot(make_K_series())
```
Next, we create an updated repressilator model, where the input $K$ value is modelled as a time-dependent parameter.
```@example noise_modelling_approaches
@parameters (K_in::typeof(make_K_series()))(..)
K_in = K_in(default_t())
repressilator = @reaction_network begin
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
    # return remake(prob; p)
    return ODEProblem(repressilator, prob.u0, prob.tspan, p)
end
oprob = ODEProblem(repressilator, u0, tend, [ps; :K_in => make_K_series()])
eprob_inputnoise = EnsembleProblem(oprob; prob_func = prob_func_Kin)
sol_inputnoise = solve(eprob_inputnoise; trajectories = 4)
plot(sol_inputnoise; idxs = :X)
```
Like in the previous two cases, this generates heterogeneous trajectories across our ensemble.


## [Investigating the mean of noisy oscillations](@id noise_modelling_approaches_model_noisy_oscillation_mean)
Finally, we will consider what happens to the mean activity of $X$ as the heterogeneous trajectories diverge from their initial condition. First we create ensemble simulations with a larger number of trajectories.
```@example noise_modelling_approaches
sol_intrinsic = solve(eprob_intrinsic, ImplicitEM(); trajectories = 100)
sol_extrinsic = solve(eprob_extrinsic; trajectories = 100)
```
Next we can use the `EnsembleSummary` interface to plot each ensemble's mean activity (as well as 5% and 95% quantiles) over time:

```@example noise_modelling_approaches
e_sumary_intrinsic = EnsembleAnalysis.EnsembleSummary(sols, 0.0:1.0:tend)
e_sumary_extrinsic = EnsembleAnalysis.EnsembleSummary(sols, 0.0:1.0:tend)
plot(e_sumary_intrinsic; label = "Intrinsic noise")
plot(e_sumary_extrinsic; label = "Extrinsic noise")
```
Here we can see that while the individual simulations (as seen in the previous sections) maintain their amplitudes, the mean amplitude across the ensemble tends towards zero. This is a well-known phenomenon in circadian biology (and a reason for considering single-cell observations, i.e. by only observing the mean it is impossible to tell whether the individual trajectories experience de-synchronisation or dampening).

---
## [References](@id noise_modelling_approaches_references)
[^1]: [Michael B. Elowitz, Arnold J. Levine, Eric D. Siggia, Peter S. Swain, *Stochastic Gene Expression in a Single Cell*, Science (2002).](https://www.science.org/doi/10.1126/science.1070919)