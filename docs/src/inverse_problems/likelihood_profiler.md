# [Practical Identifiability Analysis using LikelihoodProfiler.jl](@id likelihood_profiler)
```@raw html
<details><summary><strong>Environment setup and package installation</strong></summary>
```
The following code sets up an environment for running the code on this page.
```julia
using Pkg
Pkg.activate(; temp = true) # Creates a temporary environment, which is deleted when the Julia session ends.
Pkg.add("Catalyst")
Pkg.add("")
```
```@raw html
</details>
```
```@raw html
<details><summary><strong>Quick-start example</strong></summary>
```
The following code provides a brief example of how to compute profile likelihood diagrams using the [LikelihoodProfiler.jl](https://github.com/insysbio/LikelihoodProfiler.jl) package. Likelihood profiles is one of the most common approaches for performing [*practical identifiability analysis*](https://en.wikipedia.org/wiki/Identifiability_analysis).
```julia


# Create model (a logistic growth model).
using Catalyst
log_growth = @reaction_network begin
    r, X --> 2X
    2r/K, 2X --> X
end

# Generate some (synthetic) data for the fitting procedure.
using Distributions, OrdinaryDiffEqDefault, Plots
t_measurement = 1.0:5:150.0
u0 = [:X => 1.0]
p_true = [:r => 0.1, :K => 100.0]
oprob_true = ODEProblem(log_growth, u0, t_measurement[end], p_true)
sol_true = solve(oprob_true)
X_true = sol_true(t_measurement; idxs = :X)
σ = 5.0
X_observed = [rand(Normal(X, σ)) for X in X_true]
plot(sol_true; label = ["S (true)" "I (true)" "R (true)"])
plot!(t_measurement, X_observed; label = "X (measured)", color = 1, seriestype = :scatter)

# Recover the ground-truth parameter values from the data using PEtab.
# Using Optimization.jl-created `OptimizationProblem`s as a base for LikelihoodProfiler is also possible.
using DataFrames, Optim, PEtab
observables = [PEtabObservable(:obs_X, :X, σ)]
pest = [PEtabParameter(:r), PEtabParameter(:K)]
measurements = DataFrame(obs_id = "obs_X", time = t_measurement, measurement = X_observed)
model = PEtabModel(log_growth, observables, measurements, pest; speciemap = u0)
petab_prob = PEtabODEProblem(model)
petab_sol = calibrate_multistart(petab_prob, Optim.IPNewton(), 3)

# Compute and plot the likelihood profiler.
using LikelihoodProfiler, OptimizationLBFGSB
pl_prob = ProfileLikelihoodProblem(petab_sol, petab_prob)
pl_method = OptimizationProfiler(optimizer = LBFGSB(), stepper = FixedStep(; initial_step = 5e-3))
pl_sol = LikelihoodProfiler.solve(pl_prob, pl_method)
plot(pl_sol) # The plot shows the parameter log-scaled.
```
```@raw html
</details>
```
  \
  

During parameter fitting, parameter values are inferred from data. *Identifiability* is a concept which describes the extent to which parameter values can be uniquely determined from the data[1]. I.e., under some circumstances, multiple parameter sets can all yield good model-to-data fits. Under these circumstances, it is not possible to determine which parameter set corresponds to the true system dynamics. We term this as the model parameters being *non-identifiable*.

Identifiability can be divided into *structural* and *practical* identifiability. Structural identifiability considers only the mathematical model, and which parameters are and are not inherently identifiable due to model structure. Structural identifiability and methods for assessing it are described in the [corresponding tutorial](@ref structural_identifiability). Here we will instead describe practical identifiability, which describes our ability to identify unique parameter sets given the available data. Specifically, we will demonstrate how to compute *likelihood profiles* (a primary tool for assessing practical identifiability)[^2] using the [LikelihoodProfiler.jl](https://github.com/insysbio/LikelihoodProfiler.jl) package.

## [Introduction to practical identifiability](@id likelihood_profiler_practical_identifiability)
Before fitting parameters, we typically carry out structural identifiability analysis to show that the parameters are, in theory, identifiable give infintie quantities of noiseless data. However, available data is much more limited. Here, practical identifiability is typically carried out *after* we have fitted our parameters, aiming to determine which parameter values were identifiable (and can hence be trusted).

Let us consider the following [logistic growth model](@ref basic_CRN_library_logistic_growth):
```@example likelihood_profiler_pract_ident
using Catalyst
log_growth = @reaction_network begin
    r, X --> 2X
    2r/K, 2X --> X
end
```
From it we will generate two different data sets, one of higher quality than the other:
```@example likelihood_profiler_pract_ident
function generate_data(u0, σ, tend)
    t_measurement = range(1.0, tend; length = 50)
    p_true = Dict([:r => 0.1, :K => 100.0])
    oprob_true = ODEProblem(log_growth, u0, t_measurement[end], p_true)
    sol_true = solve(oprob_true)
    X_true = sol_true(t_measurement; idxs = :X)
    X_observed = [rand(Normal(X, σ)) for X in X_true]
    return t_measurement, X_observed, sol_true
end
u0_poor = [:X => 10.0]; σ_poor = 10.0; tend_poor = 75.0;
u0_good = [:X => 1.0]; σ_good = 1.0; tend_good = 150.0;
t_measurement_poor, X_observed_poor, sol_poor = generate_data(u0_poor, σ_poor, tend_poor)
t_measurement_good, X_observed_good, sol_good = generate_data(u0_good, σ_good, tend_good)
nothing # hide
```
Plotting the datasets we can se than one has less noise, and also capture more of the system dynamics. Using this dataset, we should be better able to recover true parameter values.
```@example likelihood_profiler_pract_ident
function plot_data(t_measurement, X_observed, sol_true)
    plot(sol_true; label = "X (true)", lw = 5)
    plot!(t_measurement, X_observed; label = "X (measured)", ms = 5, color = 1, seriestype = :scatter, alpha = 0.7)
end
plot(
    plot_data(t_measurement_poor, X_observed_poor, sol_poor),
    plot_data(t_measurement_good, X_observed_good, sol_good);
    size = (1200, 400)
)
```

Next, we use [PEtab.jl](@ref petab_parameter_fitting) to create two optimisation problems for recovering the two model parameters ($r$ and $K$) from the two synthetic data sets.
```@example likelihood_profiler_pract_ident
function make_lossfun(t_measurement, X_observed, u0, σ)
    observables = [PEtabObservable(:obs_X, :X, σ)]
    pest = [PEtabParameter(:r), PEtabParameter(:K)]
    measurements = DataFrame(obs_id = "obs_X", time = t_measurement, measurement = X_observed)
    model = PEtabModel(log_growth, observables, measurements, pest; speciemap = u0)
    return PEtabODEProblem(model)
end
petab_prob_poor = make_lossfun(t_measurement_poor, X_observed_poor, u0_poor, σ_poor)
petab_prob_good = make_lossfun(t_measurement_good, X_observed_good, u0_good, σ_good)
nothing # hide
```

From there, we could apply some optimisation algorithm to find the best-fitting parameter sets. However, when considering identifiability, we not only want to find the single best-fitting parameter set, but rather assess the space of all parameter sets which can generate good model-to-data fits. Utilising that there are only two parameters, we will [extract the loss function from the PEtab problem](@ref ), and plot it across parameter space:
```@example likelihood_profiler_pract_ident
r_range = logrange(0.1/5, 5*0.1; length = 50)
K_range = logrange(100.0/5, 5*100.0; length = 50)
losses_poor = [log10(petab_prob_poor.nllh([log10(r), log10(K)])) for r in r_range, K in K_range]
losses_good = [log10(petab_prob_good.nllh([log10(r), log10(K)])) for r in r_range, K in K_range]
plot(
    heatmap(r_range, K_range, losses_poor'; xscale = :log10, yscale = :log10, title = "Poor data"),
    heatmap(r_range, K_range, losses_good'; xscale = :log10, yscale = :log10, title = "Good data");
    size = (1200, 400)
)
```
Here, we see distinct types of loss landscapes for the two parameter sets. For the high-quality data set, the loss landscape exhibits a single peak corresponding to the ground-truth parameter set. Here, non-ground-truth parameters sets yields poor model-to-data fits and will be discarded. For the low-quality data, however, we see a "banana" of values all corresponding to good model-to-data fits. These do include the ground-truth parameter set, but also alternative ones. Since our loss function do not allow us to distinguish the true parameters from alternative ones, we say that these are non-identifiable. When fitting high-dimensional parameter sets, computing the full likelihood landscape is unfeasible, instead, so-called "likelihood profiles" are used. These are described in the next section and are essentially projections of the loss landscape onto individual parameters.

  
## [Computing profile likelihood diagrams using LikelihoodProfiler](@id likelihood_profiler_basics)

Here, we will use a combination of the [PEtab.jl](https://github.com/sebapersson/PEtab.jl) and LikelihoodProfiler.jl packages to generate a parameter fitting problem and then compute the corresponding likelihood profiles. The LikelihoodProfiler package [can also be used](https://insysbio.github.io/LikelihoodProfiler.jl/stable/rosenbrock/) in conjecture with parameter fitting problems created using e.g. [Optimization.jl](@ref optimization_parameter_fitting). For a more detailed documentation, please consider [LikelihoodProfiler's documentation](https://insysbio.github.io/LikelihoodProfiler.jl/stable/).

First, we will generate a synthetic dataset for the logistic growth model considered in the previous section:
```@example likelihood_profiler_basics
# Declare the model.
using Catalyst
log_growth = @reaction_network begin
    r, X --> 2X
    2r/K, 2X --> X
end

# generate synthetic data.
using Distributions, OrdinaryDiffEqDefault, Plots
t_measurement = 1.0:5:150.0
u0 = [:X => 1.0]
p_true = [:r => 0.1, :K => 100.0]
oprob_true = ODEProblem(log_growth, u0, t_measurement[end], p_true)
sol_true = solve(oprob_true)
X_true = sol_true(t_measurement; idxs = :X)
σ = 5.0
X_observed = [rand(Normal(X, σ)) for X in X_true]
plot(sol_true; label = ["S (true)" "I (true)" "R (true)"])
plot!(t_measurement, X_observed; label = "X (measured)", color = 1, seriestype = :scatter)
```

Next, we will use a [standard PEtab.jl workflow](@ref petab_parameter_fitting) to generate a parameter fitting optimisation problem for fitting the parameters to this data.
```@example likelihood_profiler_basics
using DataFrames, Optim, PEtab
observables = [PEtabObservable(:obs_X, :X, σ)]
pest = [PEtabParameter(:r), PEtabParameter(:K)]
measurements = DataFrame(obs_id = "obs_X", time = t_measurement, measurement = X_observed)
model = PEtabModel(log_growth, observables, measurements, pest; speciemap = u0)
petab_prob = PEtabODEProblem(model)
petab_sol = calibrate_multistart(petab_prob, Optim.IPNewton(), 3)
nothing # hide
```
Now, we use the LikelihoodProfiler package to compute the likelihood profile. This consists of three steps, (1) Generating a `LikelihoodProblem`, (2) Declaring the profiling method we want to use, and (3) computing the profiles using `solve`. It is also possible to create likelihood profiles using non-PEtab-based workflows (as described in e.g. the [LikelihoodProfiler documentation](https://insysbio.github.io/LikelihoodProfiler.jl/stable/)). In subsequent sections we provide more details on designating the profiling method, while here we will use a simple LBFGSB-based `OptimizationProfiler` method.
```@example likelihood_profiler_basics
using LikelihoodProfiler, OptimizationLBFGSB
pl_prob = ProfileLikelihoodProblem(petab_sol, petab_prob)
meth_opt = OptimizationProfiler(optimizer = LBFGSB(), stepper = FixedStep(; initial_step = 5e-3))
pl_sol = LikelihoodProfiler.solve(pl_prob, meth_opt)
nothing # hide
```
Finally, we plot the profiles using the `plot` function.
```@example likelihood_profiler_basics
plot(pl_sol)
```
Here we see, for each parameter, the corresponding likelihood profile. These describes, for each parameter, the optimal loss function value that can be achieved by optimising the full parameter set conditioned on that parameter having that specific value. I.e. the parameter $r$'s likelihood profiles value at $r = 1.0$ describes the best loss function value conditioned on $r = 1.0$. It is achieved by solving the optimisation problem while keeping $r$ fixed at that value. The full profile is achieved be re-solving the optimisation problem at all potential values of $r$. Here, if a large range of potential parameter values have a similar value to the optimum, that suggests that that parameter can attain a wide range of values while still achieved a good fit (suggesting non-identifiability). For *likelihood-based* optimisation problems (such as e.g. PEtab generates), the a likelihood threshold value can be computed such that the 95% cofindence interval for that parameter's value correspond to the intersection of the likelihood profile with that threshold. This give a direct nummeric meassure of identifiability. This threshold is plotted automatically by Likelihood in the profile. Finally, we note that, when creating PEtab-based likelihood profiles, by default the log-sclaed parameter is plotted.

  
### [Likelihood profile plotting options](@id likelihood_profiler_plotting)
  
### [Available likelihood profiling methods and options](@id likelihood_profiler_methods)


---
## [Citation](@id structural_identifiability_citation)
If you use this functionality in your research, please cite the following paper to support the authors of the LikelihoodProfiler package:
```bibtex
@article{Borisov2026, 
  doi = {10.21105/joss.09501}, 
  url = {https://doi.org/10.21105/joss.09501}, 
  year = {2026}, 
  publisher = {The Open Journal}, 
  volume = {11}, 
  number = {117}, 
  pages = {9501}, 
  author = {Borisov, Ivan and Demin, Aleksander and Metelkin, Evgeny}, 
  title = {LikelihoodProfiler.jl: Unified profile-likelihood workflows for identifiability and confidence intervals}, 
  journal = {Journal of Open Source Software} 
}
```

Similarly, if you are using the PEtab interface, please also cite the following paper:
```bibtex
@article{PEtabBioinformatics2025,
  title={PEtab.jl: advancing the efficiency and utility of dynamic modelling},
  author={Persson, Sebastian and Fr{\"o}hlich, Fabian and Grein, Stephan and Loman, Torkel and Ognissanti, Damiano and Hasselgren, Viktor and Hasenauer, Jan and Cvijovic, Marija},
  journal={Bioinformatics},
  volume={41},
  number={9},
  pages={btaf497},
  year={2025},
  publisher={Oxford University Press}
}
``` 


---
## References


[^1]: [Simpson M.J., Baker R.E., *Parameter Identiﬁability, ParameterEstimation, and Model Prediction forDifferential Equation Models*, SIAM Review (2026).](https://epubs.siam.org/doi/epdf/10.1137/24M1667968)
[^2]: [Raue A., et. al., *Structural and practical identifiability analysis of partially observed dynamical models by exploiting the profile likelihood*, Bioinformatics (2009).](https://pubmed.ncbi.nlm.nih.gov/19505944/)