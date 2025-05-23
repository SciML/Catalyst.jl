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
Pkg.add("Turing")
```
```@raw html
</details>
```
```@raw html
<details><summary><strong>Quick-start example</strong></summary>
```
The following code provides a brief example of how to simulate the chemical master equation using the [FiniteStateProjection.jl](https://github.com/SciML/FiniteStateProjection.jl) package.
```julia
# Create reaction network model (a bistable switch).
using Catalyst
sir = @reaction_network begin
    γ, S + I --> 2I
    ν, I --> R
end

# Generate some (synthetic) data for the fitting procedure).
using Distributions, OrdinaryDiffEq, Plots
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
    for idx in eachindex(data)
        I_observed[idx] ~ Normal(sol[:I][idx], σI)
    end
end

# Estimates parameter priors through Markov Chain MOnte Carlo.
using StatsPlots
n_steps = 1000
n_chains = 4
setp_oop = ModelingToolkit.setp_oop(oprob_base, [:γ, :ν])
model = sir_model(I_observed, oprob_true, setp_oop, t_measurement)
chain = sample(model, NUTS(), MCMCSerial(), n_steps, n_chains; progress = false)
plot(chain)
```
```@raw html
</details>
```
  \

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
using OrdinaryDiffEq, Plots
u0 = [:S => 999.0, :I => 1.0, :R => 0.0]
p_true = [:γ => 0.0003, :ν => 0.1]
oprob_true = ODEProblem(sir, u0, t_measurement[end], p_true)
sol_true = solve(oprob_true)
plot(sol_true; label = ["S (true)" "I (true)" "R (true)"])
```
Next, we generate some synthetic measurements from this epidemic. We make $100$, evenly spaced, measurements of the number of infected individuals ($I$). To these measurements we add some normally distributed noise (with mean at the true value and a standard deviation of $10$). 
```@example turing_paramfit
using Distributions
σI = 10.0
t_measurement = 1.0:100.0
I_observed = sol_true(t_measurement; idxs = :I)
I_observed = [rand(Normal(I, σI)) for I in I_observed]
plot!(t_measurement, I_observed; label = "I (measured)", color = 1)
```

Next, we are ready to create a Turing model (from which posteriors can be estimated). The Turing model have similarities to [the loss function utilised for normal parameter fitting workflows](@ref optimization_parameter_fitting_basics), but with a few differences:
- The declaration is prepended with the `@model` macro (which enables some specialised notation).
- The function have no input parameter values. However, its input should contain all observables (and also other potential structures used within the loss function).
- All quantities which posteriors we wish to infer (in our case parameters) are declared within the function (together with their priors).
- The function does not return a loss value. Instead it uses a special notation to compute likelihood of all observables (and from these Turing can compute a total likelihood of each parameter set).

Here, we declare our parameters on the form `p ~ Distribution(...)` where the left-hand side is the parameter and the right-hand side is its prior distribution (any distribution defined within the [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) package can be used). The likelihood of observing each observable is defined in a similar manner, i.e `o ~ Distribution(...)`. The only difference is that the observables are provided as input values to the Turing model (and the parameters are fully declared within it).

In our case, we first declare each parameter and their prior. Next we simulate the SIR model for a specific parameter set. Then, we compute the likelihood of observing out observables given the simulation.
```@example turing_paramfit
using Turing
setp_oop = ModelingToolkit.setp_oop(oprob_base, [:γ, :ν])
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
    for idx in eachindex(data)
        I_observed[idx] ~ Normal(sol[:I][idx], σI)
    end
end
nothing # hide
```
Some specific comments regarding how we have declared the model above:
- Like for [normal parameter fitting](@ref optimization_parameter_fitting_basics), we use the `maxiters = 10000` (to prevent spending long time simulating unfeasible parameter sets) and `verbose = false` (to prevent unnecessary printing of warning messages) arguments to `solve`.
- Again, we need to handle parameter sets where the model cannot be successfully simulated. Here, we use `Turing.@addlogprob! -Inf` to set an non-existent likelihood, and `return nothing` to finish further evaluation of the specific parameter set.
- Here we assume that we (correctly) know that the noise is normally distributed. However, we assume that we *do not know the standard deviation*. Instead, we make the standard deviation a third parameter (which value we infer as part of the inference process). Here we assume that the noise do not scale with $I$, however, [*affine* noise is commonly used]().
- We create am `setp_oop` function to update the parameter values in each step (as this works well with [*automatic differentiation*](@ref optimization_parameter_fitting_AD)).

Finally, we can estimate the posterior distributions of all parameters. First we generate a Turing model from the likelihood function declared above. Here, we also provide the likelihood's input values (the observables and all other required inputs, such as the base `ODEProblem`). Next, we use Turing's `sample` function to compute the posteriors.
```@example turing_paramfit
n_steps = 1000
n_chains = 4
sir_model = sir_likelihood(I_observed, oprob_true, setp_oop, t_measurement)
chain = sample(sir_model, NUTS(), MCMCSerial(), n_steps, n_chains; progress = false)
nothing # hide
```
Here, `sample`'s input is:
- The Turing model.
- Something.
- Something else.
- The number of steps to take in each chain
- The number of parallel chain.
- An optional argument that we do not wish to print progression updates. 

Additional options are discussed [here](@ref turing_parameter_fitting_other_options) and in [Turing's documentation](@ref https://turinglang.org/).

Turing contain a plotting interface for plotting the results:
```@example turing_paramfit
using StatsPlots
plot(chain)
```
Here, for each parameter, the left-hand side is shows the MCMC chain, and the right-hand side plot the posterior distribution.

### [Accessing posterior information](@id turing_parameter_fitting_basic_example_output_interfacing)

## [Encoding non-negativity in observables formulas](@id turing_parameter_fitting_nonnegative_observables)

## [Other useful options](@id turing_parameter_fitting_other_options)