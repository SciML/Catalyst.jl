# [Bayesian Parameter Fitting for ODEs using Turing.jl](@id turing_parameter_fitting)
```@raw html
<details><summary><strong>Environment setup and package installation</strong></summary>
```
The following code sets up an environment for running the code on this page.
```julia
using Pkg
Pkg.activate(".")
Pkg.add("Catalyst")
Pkg.add("Distributions")
Pkg.add("OrdinaryDiffEqDefault")
Pkg.add("Plots")
Pkg.add("StatsPlots")
Pkg.add("SymbolicIndexingInterface")
Pkg.add("Turing")
```
```@raw html
</details>
```
```@raw html
<details><summary><strong>Quick-start example</strong></summary>
```
The following code provides a minimal example of how to infer parameter posteriors from data usings the [Turing.jl](https://github.com/TuringLang/Turing.jl) package.
```julia
# Create reaction network model (a bistable switch).
using Catalyst
sir = @reaction_network begin
    γ, S + I --> 2I
    ν, I --> R
end

# Generate some (synthetic) data for the fitting procedure.
using Distributions, OrdinaryDiffEqDefault, Plots
t_measurement = 1.0:100.0
u0 = [:S => 999.0, :I => 1.0, :R => 0.0]
p_true = [:γ => 0.0003, :ν => 0.1]
oprob_true = ODEProblem(sir, u0, t_measurement[end], p_true)
sol_true = solve(oprob_true)
I_observed = sol_true(t_measurement; idxs = :I)
σI = 10.0
I_observed = [rand(Normal(I, σI)) for I in I_observed]
plot(sol_true; label = ["S (true)" "I (true)" "R (true)"])
plot!(t_measurement, I_observed; label = "I (measured)", color = 1, seriestype = :scatter)

# Create a Turing modelThe `x ~ Distribution(...)` notation both defines priors on the parameters
# which posterior we wish to infer, and the likelihoods of the observables (which are function arguments).
# If `x` is an undefined parameter, `x` becomes a prior for an estimated parameter, else,
# it is interpreted as an observable likelihood.
using Turing
@model function sir_model(I_observed, oprob_base, setp_oop, saveat)
    # Defines the parameters we wish to estimate and their prior distributions.
    γ ~ LogUniform(0.00001, 0.001)
    ν ~ LogUniform(0.01, 1.0)
    σI ~ LogUniform(0.1, 100.0)

    # Simulate the model for parameter values γ, ν. Saves the solution at the measurement times.
    p = setp_oop(oprob_base, [γ, ν])
    oprob_base = remake(oprob_base; p)
    sol = solve(oprob_base; saveat, verbose = false, maxiters = 10000)

    # If simulation was unsuccessful, the likelihood is -Inf.
    if !SciMLBase.successful_retcode(sol)
        Turing.@addlogprob! -Inf
        return nothing
    end

    # Computes the likelihood of the observations.
    for idx in eachindex(I_observed)
        I_observed[idx] ~ Normal(sol[:I][idx], σI)
    end
end

# Estimates parameter priors through Markov Chain Monte Carlo.
import SymbolicIndexingInterface # Required for `setp_oop`.
n_steps = 1000
n_chains = 4
setp_oop = SymbolicIndexingInterface.setp_oop(oprob_true, [:γ, :ν])
model = sir_model(I_observed, oprob_true, setp_oop, t_measurement)
chain = sample(model, NUTS(), MCMCThreads(), n_steps, n_chains; progress = false)

# Plots the resulting chains and posteriors.
using StatsPlots
plot(chain)
```
```@raw html
</details>
```
  \

In [previous](@ref petab_parameter_fitting) [sections](@ref optimization_parameter_fitting) we have described how to fit parameter values to data. These sections demonstrated so-called [*frequentist approaches*](https://en.wikipedia.org/wiki/Frequentist_inference), i.e. where we attempt to find the single solution which fits the data best. An alternative is to instead use a [*Bayesian approach*](https://en.wikipedia.org/wiki/Bayesian_statistics)[^1]. Here, we recognise that we cannot plausibly find the true solution, and instead attempts to accurate quantify our knowledge of the true solution as a *probability distribution* in solution space. Another hallmark of Bayesian approaches is the incorporation of *prior knowledge*, i.e., we assume that we have some prior "guess" of what the solution might be, which is incorporated in the inference of the solution distribution (also known as the *posterior distribution*). Generally, the posterior distribution cannot be computed directly, however, so called *Markov Chain Monte Carlo* (MCMC) methods can be used to generate samples from them. Here, a correctly constructed MCMC that is run for enough steps should provide a good estimate of the posterior distribution. Bayesian inference is advantageous in that it provides more information of the solution (as opposed to frequentist approaches), which also can be used as some form of identifiability analysis (@ref structural_identifiability). On the flips side, these methods are often more complicated and computationally intensive.

In Julia, Bayesian inference is primarily carried out using the [Turing.jl](https://github.com/TuringLang/Turing.jl) package. In this tutorial, we will give a brief introduction on Turing, and how to combine it with Catalyst to perform Bayesian inference on model parameters to data. A more through introduction to Turing can be found in its [documentation](https://turinglang.org/docs/getting-started/index.html). Finally, we note that PEtab (the primary package for fitting the parameters of ODEs) offer [direct support for Bayesian inference](https://sebapersson.github.io/PEtab.jl/stable/inference/). However, if you have an inverse problem that cannot be encoded using PEtab, using a Turing-based workflow (as described below) is an alternative approach.

While in frequentist parameter fitting, we can use a cost function based on likelihoods (as created by e.g. [PEtab.jl](@ref petab_parameter_fitting)) or something else (e.g. sum of square distances), Bayesian inference is heavily associated with likelihoods. Here, we will assign each unknown parameter a prior distribution. Next, we will assign each observed quantity an observation likelihood. Using these, the probability of observing our given data can be computed for any potential parameter set. Next, the posterior distribution of the true parameter set can be computed.

## [Inferring parameter posterior distributions for an ODE model using Turing](@id turing_parameter_fitting_basic_example)
For this example, we will consider a simple [SIR model of an infectious disease](@ref basic_CRN_library_sir).
```@example turing_paramfit
using Catalyst
sir = @reaction_network begin
    γ, S + I --> 2I
    ν, I --> R
end
```
From an initial condition where only a small fraction of the population is in the infected state, the model exhibits a peak of infections, after which the epidemic subsides.
```@example turing_paramfit
using OrdinaryDiffEqDefault, Plots
u0 = [:S => 999.0, :I => 1.0, :R => 0.0]
p_true = [:γ => 0.0005, :ν => 0.1]
oprob_true = ODEProblem(sir, u0, 100.0, p_true)
sol_true = solve(oprob_true)
plot(sol_true; label = ["S (true)" "I (true)" "R (true)"], lw = 4)
```
From this simulation, we generate synthetic measurements of the epidemic (on which we will demonstrate a Turing-based inference workflow). We make $25$, evenly spaced, measurements of the number of infected individuals ($I$). To these measurements we add some normally distributed noise (with mean at the true value and a standard deviation of $20$). 
```@example turing_paramfit
using Distributions
σI = 20.0
t_measurement = 4.0:4:100.0
I_observed = sol_true(t_measurement; idxs = :I)
I_observed = [rand(Normal(I, σI)) for I in I_observed]
plot!(t_measurement, I_observed; label = "I (measured)", color = 2, seriestype = :scatter)
```

Next, we are ready to create a Turing model/likelihood function (from which posteriors can be estimated). The Turing model have similarities to [the loss function utilised for normal parameter fitting workflows](@ref optimization_parameter_fitting_basics), but with a few differences:
- The declaration is prepended with the `@model` macro (which enables some specialised notation).
- The function have no input parameter values. However, its input should contain all observed quantities (and also other potential structures used within the likelihood function).
- All quantities which posteriors we wish to infer (in our case the parameters) are declared within the function (together with their priors).
- The function does not return a loss value. Instead it uses a special notation to compute likelihood of any observables (and from these Turing can compute a total likelihood of each parameter set).

Here, we declare our parameters on the form `p ~ Distribution(...)` where the left-hand side is the parameter and the right-hand side is its prior distribution (any distribution defined within the [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) package can be used). The likelihood of observing each observable is defined in a similar manner, i.e `o ~ Distribution(...)`. The only difference is that the observables have their value already declared at the time of the `o ~ Distribution(...)` notation, while parameters are declared through* the `p ~ Distribution(...)` notation.

In our case, we first declare each parameter and their priors. Next we simulate the SIR model for a specific parameter set. Then, we compute the likelihood of observing out observables given the simulation.
```@example turing_paramfit
using Turing, SymbolicIndexingInterface
setp_oop = SymbolicIndexingInterface.setp_oop(oprob_true, [:γ, :ν])
@model function sir_likelihood(I_observed, oprob_base, setp_oop, saveat)
    # Defines the parameters we wish to estimate and their prior distributions.
    γ ~ LogUniform(0.00001, 0.001)
    ν ~ LogUniform(0.01, 1.0)
    σI ~ LogUniform(0.1, 100.0)

    # Simulate the model for parameter values γ, ν. Saves the solution at the measurement times.
    p = setp_oop(oprob_base, [γ, ν])
    oprob_base = remake(oprob_base; p)
    sol = solve(oprob_base; saveat, verbose = false, maxiters = 10000)

    # If simulation was unsuccessful, the likelihood is -Inf.
    if !SciMLBase.successful_retcode(sol)
        Turing.@addlogprob! -Inf
        return nothing
    end

    # Computes the likelihood of the observations.
    for idx in eachindex(I_observed)
        I_observed[idx] ~ Normal(sol[:I][idx], σI)
    end
end
nothing # hide
```

Some specific comments regarding how we have declared the model above:
- Like for [normal parameter fitting](@ref optimization_parameter_fitting_basics), we use the `maxiters = 10000` (to prevent spending long time simulating unfeasible parameter sets) and `verbose = false` (to prevent unnecessary printing of warning messages) arguments to `solve`.
- Again, we need to handle parameter sets where the model cannot be successfully simulated. Here, we use `Turing.@addlogprob! -Inf` to set an non-existent likelihood, and `return nothing` to finish further evaluation of the specific parameter set.
- Just like for normal parameter fitting we wish to [fit parameters on a log scale](@ref optimization_parameter_fitting_log_scale). Here we do so by declaring log-scaled prior distributions.
- Here we assume that we (correctly) know that the noise is normally distributed. However, we assume that we *do not know the standard deviation*. Instead, we make the standard deviation a third parameter (which value we infer as part of the inference process). More complicated noise formulas can be used (and is sometimes even advisable[^2]).
- Like for normal parameter fitting, we use [SymbolicIndexingInterface](https://github.com/SciML/SymbolicIndexingInterface.jl)'s`setp_oop` function to update the parameter values in each step (as this works well with [*automatic differentiation*](@ref optimization_parameter_fitting_AD)).

Finally, we can estimate the posterior distributions of all parameters. First we generate a Turing model from the likelihood function declared above. Here, we also provide the likelihood's input values (the observables and all other required inputs, such as the base `ODEProblem`). Next, we use Turing's `sample` function to compute the posteriors.
```@example turing_paramfit
n_steps = 1000
n_chains = 4
sir_model = sir_likelihood(I_observed, oprob_true, setp_oop, t_measurement)
chain = sample(sir_model, NUTS(), MCMCThreads(), n_steps, n_chains; progress = false)
nothing # hide
```
Here, `sample`'s input is:
- The Turing model.
- The sampling method (here we use the [*No U-Turn Sampler*](https://www.jmlr.org/papers/volume15/hoffman14a/hoffman14a.pdf), a type of *Hamiltonian Monte Carlo sampler*(https://en.wikipedia.org/wiki/Hamiltonian_Monte_Carlo))
- The parallelisation approach (here we use `MCMCThreads`, which parallelises through [*multithreading*](https://en.wikipedia.org/wiki/Multithreading_(computer_architecture))).
- The number of steps to take in each chain
- The number of parallel chain.
- An optional argument disabling the printing of updates through the MCMC computation. 

Additional options are discussed [here](@ref turing_parameter_fitting_other_options) and in [Turing's documentation](@ref https://turinglang.org/).

Turing contain a plotting interface for plotting the results:
```@example turing_paramfit
using StatsPlots
plot(chain)
```
Here, for each parameter, the left-hand side is shows the MCMC chains, and the right-hand side plot the posterior distributions (one for each of the four chains).

### [Encoding non-negativity in observables formulas](@id turing_parameter_fitting_nonnegative_observables)
In biology, most quantities are non-negative, which is information that we wish to incorporate in our inference problem. This holds both for priors (i.e. we know that the inferred parameters are non-negative) and observables (i.e. we know that the observed quantities are non-negative). This can be encoded either by:
- Using a prior/observation formula where the distribution is strictly non-negative.
- Truncating the prior/observation formula distribution at zero.

In the previous example we used the former approach for our parameter priors - i.e. `γ ~ LogUniform(0.00001, 0.001)` is a strictly non-negative distribution, implying the information that `γ` is non-negative. However, for our observations we used `I_observed[idx] ~ Normal(sol[:I][idx], σI)`. Here, `Normal(sol[:I][idx], σI)` suggests that negative number of infected cases can be observed. While the way with which we generated our synthetic data mean that this actually could happen, for real applications, we might want to encode non-negativity in the distribution. A simple alternative is to truncate the distribution at zero using the [`truncated`](https://juliastats.org/Distributions.jl/stable/truncate/#Distributions.truncated) function. Here we can use
```julia
I_observed[idx] ~ truncated(Normal(sol[:I][idx], σI), 0.0, Inf)
```
to create a version of our normal distribution that is truncated at *0* and infinity.

### [Accessing posterior information](@id turing_parameter_fitting_basic_example_output_interfacing)
Say that we want to sample a parameter set from the computed posterior distribution. Here, we can use the following syntax to sample a single vector with the values of `γ`, `ν`, and `σI` (in that order):
```@example turing_paramfit
collect(chain.value[rand(1:n_steps), 1:3, rand(1:n_chains)])
```

We can use this to e.g. draw $10$ random parameter sets from the posterior distribution, simulate the model for these parameter sets, and plot the resulting ensemble simulation. For this, we will create an `EnsembleProblem` from our `ODEProblem` using the approach described [here](@ref @id ensemble_simulations_varying_conditions).
```@example turing_paramfit
function prob_func(prob, _, _)
    γ, ν = collect(chain.value[rand(1:n_steps), 1:2, rand(1:n_chains)])
    remake(prob; p = [:γ => γ, :ν => ν])
end
eprob = EnsembleProblem(oprob_true; prob_func)
sols = solve(eprob; trajectories = 10)
plot(sols; color = [1 2 3])
```
Here, we can see that all parameter sets sampled from the posterior yields very similar simulations.


---
## [Citations](@id turing_citations)
If you use this functionality in your research, [in addition to Catalyst](@ref doc_index_citation), please cite the following paper to support the authors of the Turing.jl package.
```
@article{10.1145/3711897,
    author = {Fjelde, Tor Erlend and Xu, Kai and Widmann, David and Tarek, Mohamed and Pfiffer, Cameron and Trapp, Martin and Axen, Seth D. and Sun, Xianda and Hauru, Markus and Yong, Penelope and Tebbutt, Will and Ghahramani, Zoubin and Ge, Hong},
    title = {Turing.jl: a general-purpose probabilistic programming language},
    year = {2025},
    publisher = {Association for Computing Machinery},
    address = {New York, NY, USA},
    url = {https://doi.org/10.1145/3711897},
    doi = {10.1145/3711897},
    note = {Just Accepted},
    journal = {ACM Trans. Probab. Mach. Learn.},
    month = feb,
}
```

---
## References
[^1]: [Andrew Gelman et al. *Bayesian Data Analysis*, CRC Press (2013).](https://sites.stat.columbia.edu/gelman/book/BDA3.pdf)
[^2]: [Johannes H. Proost *Combined proportional and additive residual error models in population pharmacokinetic modelling*, European Journal of Pharmaceutical Sciences (2017).](https://www.sciencedirect.com/science/article/pii/S092809871730249X)