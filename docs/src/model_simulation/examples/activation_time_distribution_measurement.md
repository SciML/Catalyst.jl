# [Measuring the Distribution of System Activation Times](@id activation_time_distribution_measurement)

In this example we will consider a model which, while initially inactive, activates in response to an input. The model is *stochastic*, causing the activation times to be *random*. By combining events, callbacks, and stochastic ensemble simulations, we will measure the probability distribution of the activation times (so-called [*first passage times*](https://en.wikipedia.org/wiki/First-hitting-time_model)).

Our model will be a version of the [simple self-activation loop](@ref basic_CRN_library_self_activation) (the ensemble simulations of which we have [considered previously](@ref ensemble_simulations_monte_carlo)). Here, we will consider the activation threshold parameter ($K$) to be activated by an input (at an input time $t = 0$). Before the input, $K$ is very large (essentially keeping the system inactive). After the input, it is reduced to a lower value (which permits the system to activate). We will model this using two additional parameters ($Kᵢ$ and $Kₐ$, describing the pre and post-activation values of $K$, respectively). Initially, $K$ will [default to](@ref dsl_advanced_options_default_vals) $Kᵢ$. Next, at the input time ($t = 0$), an event will change $K$'s value to $Kᵢ$.
```@example activation_time_distribution_measurement
using Catalyst
sa_model = @reaction_network begin
    @parameters Kᵢ Kₐ K=Kᵢ
    @discrete_events [0.0] => [K ~ Kₐ]
    v0 + hill(X,v,K,n), 0 --> X
    deg, X --> 0
end
```
Next, to perform stochastic simulations of the system we will create an `SDEProblem`. Here, we will need to assign parameter values to $Kᵢ$ and $Kₐ$, but not to $K$ (as its value is controlled by its default and the event). Also note that we start the simulation at a time $t = -200 < 0$. This ensures that by the input time ($t = 0$), the system has (more or less) reached its (inactive) steady state distribution. It also means that the activation time can be measured exactly as the simulation time at the time of activation (as this will be the time from the input at $t = 0$).
```@example activation_time_distribution_measurement
u0 = [:X => 10.0]
tspan = (-200.0, 2000.0)
ps = [:v0 => 0.1, :v => 2.5, :Kᵢ => 1000.0, :Kₐ => 40.0, :n => 3.0, :deg => 0.01]
sprob = SDEProblem(sa_model, u0, tspan, ps)
nothing # hide
```
We can now create a simple `EnsembleProblem` and perform an ensemble simulation (as described [here](@ref ensemble_simulations)). Please note that the system has an event which modifies its parameters, hence we must add the `safetycopy = true` argument to `EnsembleProblem` (else, subsequent simulations would start with $K = Kₐ$).
```@example activation_time_distribution_measurement
using Plots, StochasticDiffEq
eprob = EnsembleProblem(sprob; safetycopy = true)
esol = solve(eprob, ImplicitEM(); trajectories = 10)
plot(esol)
```
Here we see how, after the input time, the system (randomly) switches from the inactive state to the active one (several examples of this, bistability-based, activation have been studied in the literature, both in models and experiments[^1][^2]).

Next, we wish to measure the distribution of these activation times. First we will create a [callback](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/) which terminates the simulation once it has reached a threshold. This both ensures that we do not have to expend unnecessary computer time on the simulation after its activation, and also enables us to measure the activation time as the final time point of the simulation. Here we will use a [discrete callback](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/#SciMLBase.DiscreteCallback). By looking at the previous simulations, we determine $X = 100$ as a suitable activation threshold. We use the `terminate!` function to terminate the simulation once this value has been reached.
```@example activation_time_distribution_measurement
condition(u, t, integrator) = integrator[:X] > 100.0
affect!(integrator) = terminate!(integrator)
callback = DiscreteCallback(condition, affect!)
nothing # hide
```
Next, we will perform our ensemble simulation. By default, for each simulation, these save the full trajectory. Here, however, we are only interested in the activation time. This enables us to utilise an [*output function*](https://docs.sciml.ai/DiffEqDocs/dev/features/ensemble/#Building-a-Problem). This will be called at the end of every single simulation and determine what to save (it can also be used to potentially rerun individual simulations, however, we will not use this feature here). Here we create an output function which saves only the simulation's final time point. We also make it throw a warning if the simulation reaches the end of the simulation time frame ($t = 2000$) (this indicates a simulation that never activated, the occurrences of such simulation would cause us to underestimate the activation times). Finally, just like previously, we must set `safetycopy = true`.
```@example activation_time_distribution_measurement
function output_func(sol, i)
    (sol.t[end] == tspan[2]) && @warn "A simulation did not activate during the given time span."
    return (sol.t[end], false)
end
eprob = EnsembleProblem(sprob; output_func, safetycopy = true)
esol = solve(eprob, ImplicitEM(); trajectories = 250, callback)
nothing # hide
```
Finally, we can plot the distribution of activation times. For this, we will use the [`histogram`](https://docs.juliaplots.org/latest/series_types/histogram/) function (with the `normalize = true` argument to create a probability density function). An alternative we also recommend is [StatsPlots.jl](https://docs.juliaplots.org/latest/generated/statsplots/)'s `density` function (which creates a smoothed histogram that is also easier to combine with other plots). The input to `density` is the activation times (which our output function has saved to `esol.u`).
```@example activation_time_distribution_measurement
histogram(esol.u; normalize = true, label = "Activation time distribution", xlabel = "Activation time")
```
Here we that the activation times take some form of long tail distribution (for non-trivial models like this one, it is generally not possible to identify the activation times as any known statistical distribution).

---

## References

[^1]: [David Frigola, Laura Casanellas, José M. Sancho, Marta Ibañes, *Asymmetric Stochastic Switching Driven by Intrinsic Molecular Noise*, PLoS One (2012).](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0031407)
[^2]: [Christian P Schwall et al., *Tunable phenotypic variability through an autoregulatory alternative sigma factor circuit*, Molecular Systems Biology (2021).](https://www.embopress.org/doi/full/10.15252/msb.20209832)
