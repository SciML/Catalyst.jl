# [Approaches for modelling system noise](@id noise_modelling_approaches)
Catalyst's primary tools for modelling stochasticity includes the creation of `SDEProblem`s or `JumpProblem`s from reaction network models. However, other approaches for incorporating noise in a model exists, some of which will be discussed here. We will first consider *intrinsic* and *extrinsic* noise. These are well-established terms, both of which we will describe in detail below (however, to our knowledge, no generally agreed upon definition of these exists). The final approach, the utilisation of a noisy input process to an otherwise deterministic system, is less frequently used. However, as it has been used in the literature, we will demonstrate it here as well.

Finally, we note that these approaches can all be combined. E.g. modelling intrinsic noise using an SDE can be combined with random parameter values, and a randomised input process can be added to this.

!!! note
    Here we use intrinsic and extrinsic noise as description of two of our modelling approaches. It should be noted that these are established terminologies for noisy biological systems [REF]. However, our use of these terms to describe different approaches for modelling noise is only inspired by this terminology, and nothing that is established in the field.

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
Generally, intrinsic noise is noise inherent to the system itself. This means that it cannot be controlled for, or filtered out by, experimental settings. Low-copy number cellular systems, were reaction occurs due to the encounters of molecules due to random diffusion, is an example of intrinsic noise. In practise, this can be precisely modelled through [SDE](@ref simulation_intro_SDEs) (chemical Langevin equations) or [jump](@ref simulation_intro_jumps) (stochastic chemical kinetics) simulations. 

Here we will model intrinsic noise through SDEs. First we create a `SDEProblem` using the standard approach.
```@example noise_modelling_approaches
u0 = [:X => 45.0, :Y => 20.0, :Z => 20.0]
tend = 200.0
ps = [:v => 10.0, :K => 20.0, :n => 3, :d => 0.1]
sprob = SDEProblem(repressilator, u0, tend, ps)
nothing # hide
```
Next, to illustrate the system's noisiness, we will perform multiple simulations. We do this by [creating an `EnsembleProblem`](@ref ensemble_simulations_monte_carlo). From it, we make, and plot, 4 simulations.
```@example noise_modelling_approaches
using StochasticDiffEq, Plots
eprob_intrinsic = EnsembleProblem(sprob)
sol_intrinsic = solve(eprob_intrinsic, ImplicitEM(); trajectories = 4)
plot(sol_intrinsic; idxs = :X)
```

## [Using extrinsic noise](@id noise_modelling_approaches_model_extrinsic)

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

```@example noise_modelling_approaches
using OrdinaryDiffEqDefault
oprob = ODEProblem(repressilator, u0, tend, ps)
eprob_extrinsic = EnsembleProblem(oprob; prob_func)
sol_extrinsic = solve(eprob_extrinsic; trajectories = 4)
plot(sol_extrinsic; idxs = :X)
```

## [Using a noisy input process](@id noise_modelling_approaches_model_input_noise)

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

## [Investigating the mean of noisy oscillations](@id noise_modelling_approaches_model_noisy_oscillation_mean)

```@example noise_modelling_approaches

```

```@example noise_modelling_approaches
```