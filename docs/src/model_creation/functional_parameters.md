# [Inputs and time-dependent (or functional) parameters](@id time_dependent_parameters)

Catalyst supports the usage of "functional parameters". In practice, these are parameters that are given by (typically) time-dependent functions (they can also depend on e.g. species values, as discussed [here](@ref functional_parameters_sir)). They are a way to inject custom functions into models. Functional parameters can be used when rates depend on real data, or to represent complicated functions (which use e.g. `for` loops or random number generation). Here, the function's values are declared as a data interpolation (which interpolates discrete samples to a continuous function). This is then used as the functional parameter's value in the simulation. This tutorial first shows how to create time-dependent functional parameters, and then gives an example where the functional parameter depends on a species value.

An alternative approach for representing complicated functions is by [using `@register_symbolic`](@ref dsl_description_nonconstant_rates_function_registration).

## [Basic example](@id functional_parameters_basic_example)

Let us first consider an easy, quick-start example (the next section will discuss what is going on in more detail). We will consider a simple [birth-death model](@ref basic_CRN_library_bd), but where the birth rate is determined by an input parameter (for which the value depends on time). First, we [define the input parameter programmatically](@ref programmatic_CRN_construction), and its values across all time points using the [DataInterpolations.jl](https://github.com/SciML/DataInterpolations.jl) package. In this example we will use the input function $pIn(t) = (2 + t)/(1 + t)$. Finally, we plot the input function, demonstrating how while it is defined at discrete points, DataInterpolations.jl generalises this to a continuous function.
```@example functional_parameters_basic_example
using Catalyst, DataInterpolations, Plots
t = default_t()
tend = 10.0
ts = collect(0.0:0.05:tend)
spline = LinearInterpolation((2 .+ ts) ./ (1 .+ ts), ts)
@parameters (pIn::typeof(spline))(..)
plot(spline)
```
Next, we create our model, [interpolating](@ref dsl_advanced_options_symbolics_and_DSL_interpolation) the input parameter into the `@reaction_network` declaration.
```@example functional_parameters_basic_example
bd_model = @reaction_network begin
    $pIn(t), 0 --> X
    d, X --> 0
end
```
Finally, we can simulate our model as normal (but where we set the value of the `pIn` parameter to our interpolated data).
```@example functional_parameters_basic_example
using OrdinaryDiffEqDefault
u0 = [:X => 0.5]
ps = [:d => 2.0, :pIn => spline]
oprob = ODEProblem(bd_model, u0, tend, ps)
sol = solve(oprob)
plot(sol)
```
!!! note
    For this simple example, $(2 + t)/(1 + t)$ could have been used directly as a reaction rate (or written as a normal function), technically making the functional parameter approach unnecessary. However, here we used this function as a simple example of how discrete data can be made continuous using DataInterpolations, and then have its values inserted using a (functional) parameter.


## [Inserting a customised, time-dependent, input](@id functional_parameters_circ_rhythm)

Let us now go through everything again, but providing some more details. Let us first consider the input parameter. We have previously described how a [time-dependent rate can model a circadian rhythm](@ref dsl_description_nonconstant_rates_time). For real applications, due to e.g. clouds, sunlight is not a perfect sine wave. Here, a common solution is to take real sunlight data from some location and use in the model. Here, we will create synthetic (noisy) data as our light input:
```@example functional_parameters_circ_rhythm
using Plots
tend = 120.0
ts = collect(0.0:1.0:tend)
light = sin.(ts/6) .+ 1
light = [max(0.0, l - rand()) for l in light]
plot(ts, light; seriestype = :scatter, label = "Experienced light")
```
Now this input is only actually defined at the sample points, making it incompatible with a continuous ODE simulation. To enable this, we will use the DataInterpolations package to create an interpolated version of this data, which forms the actual input:
```@example functional_parameters_circ_rhythm
using DataInterpolations
interpolated_light = LinearInterpolation(light, ts)
plot(interpolated_light)
```
We are now ready to declare our model. We will consider a protein with an active and an inactive form ($Pₐ$ and $Pᵢ$) where the activation is driven by the presence of sunlight. In this example we we create our model using the [programmatic approach](@ref programmatic_CRN_construction). Do note the special syntax we use to declare our input parameter, where we both designate it as a generic function and its type as the type of our interpolated input. Also note that, within the model, we mark the input parameter (`light_in`) as a function of `t`.
```@example functional_parameters_circ_rhythm
using Catalyst
t = default_t()
in_type = typeof(interpolated_light)
@parameters kA kD (light_in::in_type)(..)
@species Pₐ(t) Pᵢ(t)
rxs = [
    Reaction(kA*light_in(t), [Pᵢ], [Pₐ]),
    Reaction(kD, [Pₐ], [Pᵢ])
]
@named rs = ReactionSystem(rxs, t)
rs = complete(rs)
```
Now we can simulate our model. Here, we use the interpolated data as the input parameter's value.
```@example functional_parameters_circ_rhythm
using OrdinaryDiffEqDefault
u0 = [Pᵢ => 1.0, Pₐ => 0.0]
ps = [kA => 1.5, kD => 1.0, light_in => interpolated_light]
oprob = ODEProblem(rs, u0, tend, ps)
sol = solve(oprob)
plot(sol)
```

### [Interpolating the input into the DSL](@id functional_parameters_circ_rhythm_dsl)

It is possible to use time-dependent inputs when creating models [through the DSL](@ref dsl_description) as well. However, it can still be convenient to declare the input parameter programmatically as above. Next, we can [interpolate](@ref dsl_advanced_options_symbolics_and_DSL_interpolation) it into our DSL-declaration (ensuring to also make it a function of `t`):
```@example functional_parameters_circ_rhythm
rs_dsl = @reaction_network rs begin
    (kA*$light_in(t), kD), Pᵢ <--> Pₐ
end
```
We can confirm that this model is identical to our programmatic one (and should we wish to, we can simulate it using identical syntax).
```@example functional_parameters_circ_rhythm
rs == rs_dsl
```

## [Non-time functional parameters](@id functional_parameters_sir)

Previously we have demonstrated functional parameters that are functions of time. However, functional parameters can be functions of any variable (however, currently, more than one argument is not supported). Here we will demonstrate this using a [SIR model](@ref basic_CRN_library_sir), but instead of having the infection rate scale linearly with the number of infected individuals, we instead assume we have measured data of the infection rate (as dependent on the number of infected individuals) and wish to use this instead. Normally we use the following infection reaction in the SIR model:
```julia
@reaction k1, S + I --> 2I
```
For ODE models, this would give the same equations as
```julia
@reaction k1*I, S --> I
```
Due to performance reasons (especially for jump simulations) the former approach is *strongly* encouraged. Here, however, we will assume that we have measured real data of how the number of infected individuals affects the infection rate, and wish to use this in our model, i.e. something like this:
```julia
@reaction k1*i_rate(I), S --> I
```
This is a case where we can use a functional parameter, whose value depends on the value of $I$.

We start by declaring the functional parameter that describes how the infection rate depends on the number of infected individuals. We also plot the measured infection rate, and compare it to the theoretical rate usually used in the SIR model.
```@example functional_parameters_sir
using DataInterpolations, Plots
I_grid = collect(0.0:5.0:100.0)
I_measured = 300.0 *(0.8*rand(length(I_grid)) .+ 0.6) .* I_grid ./ (300 .+ I_grid)
I_rate = LinearInterpolation(I_measured, I_grid)
plot(I_rate; label = "Measured infection rate")
plot!(I_grid, I_grid; label = "Normal SIR infection rate")
```
Next, we create our model (using the DSL approach).
```@example functional_parameters_sir
using Catalyst
@parameters (inf_rate::typeof(I_rate))(..)
sir = @reaction_network rs begin
    k1*$inf_rate(I), S --> I
    k2, I --> R
end
nothing # hide
```
Finally, we can simulate our model.
```@example functional_parameters_sir
using OrdinaryDiffEqDefault
u0 = [:S => 99.0, :I => 1.0, :R => 0.0]
ps = [:k1 => 0.002, :k2 => 0.01, :inf_rate => I_rate]
oprob = ODEProblem(sir, u0, 250.0, ps)
sol = solve(oprob)
plot(sol)
```
!!! note
    When declaring a functional parameter of time, it is easy to set its grid values (i.e. ensure they range from the first to the final time point). For Functional parameters that depend on species concentrations it is trickier, and one must make sure that any potential input-species values that can occur during the simulation are represented in the interpolation.

### [Using different data interpolation approaches](@id functional_parameters_interpolation_algs)

Up until now we have used [linear interpolation](https://en.wikipedia.org/wiki/Linear_interpolation) of our data. However, DataInterpolations.jl [supports other interpolation methods](https://docs.sciml.ai/DataInterpolations/stable/methods/). To demonstrate these we here generate a data set, and then show the linear, cubic, and constant interpolations:
```@example functional_parameters_interpolation_algs
using DataInterpolations, Plots
xs = collect(0.0:1.0:10.0)
ys = xs ./ (5*rand(length(xs)) .+ xs)
spline_linear = LinearInterpolation(ys, xs)
spline_cubuc = CubicSpline(ys, xs)
spline_const = ConstantInterpolation(ys, xs)
plot(spline_linear)
plot!(spline_cubuc)
plot!(spline_const)
```
Finally, DataInterpolations.jl also allows various [extrapolation methods](https://docs.sciml.ai/DataInterpolations/stable/extrapolation_methods/).
