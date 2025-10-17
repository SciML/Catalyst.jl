# [Global Sensitivity Analysis](@id global_sensitivity_analysis)

*Global sensitivity analysis* (GSA) is used to study the sensitivity of a function's outputs with respect to its input[^1]. Within the context of chemical reaction network modelling it is primarily used for two purposes:

- When fitting a model's parameters to data, it can be applied to the cost function of the optimisation problem. Here, GSA helps determine which parameters do, and do not, affect the model's fit to the data. This can be used to identify parameters that are less relevant to the observed data.
- [When measuring some system behaviour or property](@ref behaviour_optimisation), it can help determine which parameters influence that property. E.g. for a model of a biofuel-producing circuit in a synthetic organism, GSA could determine which system parameters have the largest impact on the total rate of biofuel production.

GSA can be carried out using the [GlobalSensitivity.jl](https://github.com/SciML/GlobalSensitivity.jl) package. This tutorial contains a brief introduction of how to use it for GSA on Catalyst models, with [GlobalSensitivity providing a more complete documentation](https://docs.sciml.ai/GlobalSensitivity/stable/).

### [Global vs local sensitivity](@id global_sensitivity_analysis_global_vs_local_sensitivity)

A related concept to global sensitivity is *local sensitivity*. This, rather than measuring a function's sensitivity (with regards to its inputs) across its entire (or large part of its) domain, measures it at a specific point. This is equivalent to computing the function's gradients at a specific point in phase space, which is an important routine for most gradient-based optimisation methods (typically carried out through [*automatic differentiation*](https://en.wikipedia.org/wiki/Automatic_differentiation)). For most Catalyst-related functionalities, local sensitivities are computed using the [SciMLSensitivity.jl](https://github.com/SciML/SciMLSensitivity.jl) package. While certain GSA methods can utilise local sensitivities, this is not necessarily the case.

While local sensitivities are primarily used as a subroutine of other methodologies (such as optimisation schemes), it also has direct uses. E.g., in the context of fitting parameters to data, local sensitivity analysis can be used to, at the parameter set of the optimal fit, determine the cost function's sensitivity to the system parameters.

## [Basic example](@id global_sensitivity_analysis_basic_example)

We will consider a simple [SEIR model of an infectious disease](https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology). This is an expansion of the classic [SIR model](@ref basic_CRN_library_sir) with an additional *exposed* state, $E$, denoting individuals who are latently infected but currently unable to transmit their infection to others.

```@example gsa_1
using Catalyst
seir_model = @reaction_network begin
    10^β, S + I --> E + I
    10^a, E --> I
    10^γ, I --> R
end
```

We will study the peak number of infected cases's ($max(I(t))$) sensitivity to the system's three parameters. We create a function which simulates the system from a given initial condition and measures this property:

```@example gsa_1
using OrdinaryDiffEqDefault

u0 = [:S => 999.0, :I => 1.0, :E => 0.0, :R => 0.0]
p_dummy = [:β => 0.0, :a => 0.0, :γ => 0.0]
oprob_base = ODEProblem(seir_model, u0, (0.0, 10000.0), p_dummy)

function peak_cases(p)
    ps = [:β => p[1], :a => p[2], :γ => p[3]]
    oprob = remake(oprob_base; p = ps)
    sol = solve(oprob; maxiters = 100000, verbose = false)
    SciMLBase.successful_retcode(sol) || return Inf
    return maximum(sol[:I])
end
nothing # hide
```

Now, GSA can be applied to our `peak_cases` function using GlobalSensitivity's `gsa` function. It takes 3 mandatory inputs:

- The function for which we wish to carry out GSA.
- A method with which we wish to carry out GSA.
- A domain on which we carry out GSA. This is defined by a vector, which contains one two-valued Tuple for each parameter. These Tuples contain a lower and an upper bound for their respective parameter's value.

E.g., here we carry out GSA using [Morris's method](https://en.wikipedia.org/wiki/Morris_method):

```@example gsa_1
using GlobalSensitivity
global_sens = gsa(peak_cases, Morris(), [(-3.0,-1.0), (-2.0,0.0), (-2.0,0.0)])
nothing # hide
```

on the domain $10^β ∈ (-3.0,-1.0)$, $10^a ∈ (-2.0,0.0)$, $10^γ ∈ (-2.0,0.0)$ (which corresponds to $β ∈ (0.001,0.1)$, $a ∈ (0.01,1.0)$, $γ ∈ (0.01,1.0)$). The output of `gsa` varies depending on which GSA approach is used. GlobalSensitivity implements a range of methods for GSA. Below, we will describe the most common ones, as well as how to apply them and interpret their outputs.

!!! note
    We should make a couple of notes about the example above:
    - Here, we write our parameters on the forms $10^β$, $10^a$, and $10^γ$, which transforms them into log-space. As [previously described](@ref optimization_parameter_fitting_log_scale), this is advantageous in the context of inverse problems such as this one.
    - For GSA, where a function is evaluated a large number of times, it is ideal to write it as performant as possible. Hence, we initially create a base `ODEProblem`, and then apply the [`remake`](@ref simulation_structure_interfacing_problems_remake) function to it in each evaluation of `peak_cases` to generate a problem which is solved for that specific parameter set.
    - Again, as [previously described in other inverse problem tutorials](@ref optimization_parameter_fitting_basics), when exploring a function over large parameter spaces, we will likely simulate our model for unsuitable parameter sets. To reduce time spent on these, and to avoid excessive warning messages, we provide the `maxiters = 100000` and `verbose = false` arguments to `solve`.
    - As we have encountered in [a few other cases](@ref optimization_parameter_fitting_basics), the `gsa` function is not able to take parameter inputs of the map form usually used for Catalyst. Hence, as a first step in `peak_cases` we convert the parameter vector to this form. Next, we remember that the order of the parameters when we e.g. evaluate the GSA output, or set the parameter bounds, corresponds to the order used in `ps = [:β => p[1], :a => p[2], :γ => p[3]]`.

## [Sobol's method-based global sensitivity analysis](@id global_sensitivity_analysis_sobol)

The most common method for GSA is [Sobol's method](https://en.wikipedia.org/wiki/Variance-based_sensitivity_analysis). This can be carried out using:

```@example gsa_1
global_sens = gsa(peak_cases, Sobol(), [(-3.0,-1.0), (-2.0,0.0), (-2.0,0.0)]; samples = 500)
nothing # hide
```

Note: when `Sobol()` is used as the method, the `samples` argument must also be used.

Sobol's method computes so-called *Sobol indices*, each measuring some combination of input's effect on the output. Here, when `Sobol()` is used, the *first order*, *second order*, and *total order* Sobol indices are computed. These can be accessed through the following fields:

- `global_sens.S1`: A vector where the i'th element is the output's sensitivity to variations in the i'th input.
- `global_sens.S2`: A matrix where element i-j contains the output's sensitivity to simultaneous variations in the i'th and j'th inputs.
- `global_sens.ST`: A vector where the i'th element is the output's sensitivity to any simultaneous variation of any combination of inputs that contain the i'th input. While only the first and second-order (and the total) Sobol indices are computed, the total order index compounds the information contained in Sobol indices across all orders.

We can plot the first-order Sobol indices to analyse their content:

```@example gsa_1
using Plots
bar(["β", "a", "γ"], global_sens.S1; group = ["β", "a", "γ"], fillrange = 1e-3)
```

Here, we see that $β$ has a relatively low effect on the peak in infected cases, as compared to $a$ and $γ$. Plotting the total order indices suggests the same:

```@example gsa_1
bar(["β", "a", "γ"], global_sens.ST; group = ["β", "a", "γ"], fillrange = 1e-3)
```

GlobalSensitivity implements several versions of Sobol's method, and also provides several options. These are described [here](https://docs.sciml.ai/GlobalSensitivity/stable/methods/sobol/). Specifically, it is often recommended to, due to its quick computation time, use the related extended Fourier amplitude sensitivity test (EFAST) version. We can run this using:

```@example gsa_1
global_sens = gsa(peak_cases, eFAST(), [(-3.0,-1.0), (-2.0,0.0), (-2.0,0.0)]; samples = 500)
nothing # hide
```

It should be noted that when EFAST is used, only the first and total-order Sobol indices are computed (and not the second-order ones).

## [Morris's method-based global sensitivity analysis](@id global_sensitivity_analysis_morris)

An alternative to using Sobol's method is to use [Morris's method](https://en.wikipedia.org/wiki/Morris_method). The syntax is similar to previously (however, the `samples` argument is no longer required):

```@example gsa_1
global_sens = gsa(peak_cases, Morris(), [(-3.0,-1.0), (-2.0,0.0), (-2.0,0.0)])
nothing # hide
```

Morris's method computes, for parameter samples across parameter space, their *elementary effect* on the output. Next, the output's sensitivity with respect to each parameter is assessed through various statistics on these elementary effects. In practice, the following two fields are considered:

- `global_sens.means_star` (called $μ*$): Measures each parameter's influence on the output. A large $μ*$ indicates a parameter to which the output is sensitive.
- `global_sens.variances`: Measures the variance of each parameter's influence on the output. A large variance suggests that a parameter's influence on the output is highly dependent on other parameter values.

We can check these values for our example:

```@example gsa_1
mean_star_plot = bar(["β" "a" "γ"], global_sens.means_star; labels=["β" "a" "γ"], title="μ*")
variances_plot = bar(["β" "a" "γ"], global_sens.variances; labels=["β" "a" "γ"], title="σ²")
plot(mean_star_plot, variances_plot)
```

As previously, we note that the peak number of infected cases is more sensitive to $a$ and $γ$ than to $β$.

!!! note
    The syntax for plotting the output using Sobol's and Morris's methods is slightly different. The reason is that `global_sens.means_star` and `global_sens.variances` (for Morris's method) are 1x3 Matrices, while for Sobol's method, `global_sens.S1` and `global_sens.ST` are length-3 vectors.

Generally, Morris's method is computationally less intensive, and has easier to interpret output, as compared to Sobol's method. However, if computational resources are available, Sobol's method is more comprehensive.

## [Other global sensitivity analysis methods](@id global_sensitivity_analysis_other_methods)

GlobalSensitivity also implements additional methods for GSA, more details on these can be found in the [package's documentation](https://docs.sciml.ai/GlobalSensitivity/stable/).

## [Global sensitivity analysis for non-scalar outputs](@id global_sensitivity_analysis_nonscalars)

Previously, we have demonstrated GSA on functions with scalar outputs. However, it is also possible to apply it to functions with vector outputs. Let us consider our previous function, but where it provides both the peak number of exposed *and* infected individuals:

```@example gsa_1
function peak_cases_2(p)
    ps = [:β => p[1], :a => p[2], :γ => p[3]]
    oprob = remake(oprob_base; p = ps)
    sol = solve(oprob; maxiters = 100000, verbose = false)
    SciMLBase.successful_retcode(sol) || return Inf
    return [maximum(sol[:E]), maximum(sol[:I])]
end
nothing # hide
```

We can apply `gsa` to this function as previously:

```@example gsa_1
global_sens = gsa(peak_cases_2, Morris(), [(-3.0,-1.0), (-2.0,0.0), (-2.0,0.0)])
nothing # hide
```

however, each output field is now a multi-row matrix, containing one row for each of the outputs. E.g., we have

```@example gsa_1
global_sens.means_star
```

Here, the function's sensitivity is evaluated with respect to each output independently. Hence, GSA on `peak_cases_2` is equivalent to first carrying out GSA on a function returning the peak number of exposed individuals, and then on one returning the peak number of infected individuals.

---

## [Citations](@id global_sensitivity_analysis_citations)

If you use this functionality in your research, [in addition to Catalyst](@ref doc_index_citation), please cite the following paper to support the authors of the GlobalSensitivity.jl package:

```bibtex
@article{dixit2022globalsensitivity,
  title={GlobalSensitivity. jl: Performant and Parallel Global Sensitivity Analysis with Julia},
  author={Dixit, Vaibhav Kumar and Rackauckas, Christopher},
  journal={Journal of Open Source Software},
  volume={7},
  number={76},
  pages={4561},
  year={2022}
}
```

---

## References

[^1]: [Saltelli, A et al. *Global Sensitivity Analysis. The Primer*, Wiley (2008).](http://www.andreasaltelli.eu/file/repository/A_Saltelli_Marco_Ratto_Terry_Andres_Francesca_Campolongo_Jessica_Cariboni_Debora_Gatelli_Michaela_Saisana_Stefano_Tarantola_Global_Sensitivity_Analysis_The_Primer_Wiley_Interscience_2008_.pdf)
