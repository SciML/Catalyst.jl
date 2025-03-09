# [Inputs and time-dependent (or functional) parameters](@id time_dependent_parameters)
Catalyst supports the usage of "functional parameters". In practice these are not parameters, but a way to inject custom functions into models. This can be used when rates depend on real data, or to represent complicated functions ( which uses e.g. `for` loops or random number generation). Here, the function's values are declared as a data interpolation, which is then used as the functional parameter's value in the simulation. On this page, we first show how to create time-dependent functional parameters, and then give an example where the functional parameter depends on a species value.

## [Basic example](@id functional_parameters_basic_example)
Let us first consider an easy, quick-start example. We will consider a simple [birth-death model](@ref basic_CRN_library_bd), but where the birth rate is determined by an input parameter (which value depends on time). First, we [define the input parameter programmatically](@ref programmatic_CRN_construction), and its values across all time values using the [DataInterpolations.jl](https://github.com/SciML/DataInterpolations.jl) package. In this example we will use the input function $pIn(t) = (2 + t)/(1 + t)$.
```@example functional_parameters_basic_example
using Catalyst, DataInterpolations
t = default_t()
tend = 10.0
ts = collect(0.0:0.01:tend)
spline = LinearInterpolation((2 .+ ts) ./ (1 .+ ts), ts)
@parameters (pIn::typeof(spline))(..)
nothing # hide
```
Next, we create our model, [interpolating](@ref dsl_advanced_options_symbolics_and_DSL_interpolation) the input parameter into it (making it a function of `t`).
```@example functional_parameters_basic_example
input = pIn(default_t())
bd_model = @reaction_network begin
    $input, 0 --> X
    d, X --> 0
end
```
Finally, we can simulate our model as normal, but setting the value of the `pIn` parameter to our interpolated data.
```@example functional_parameters_basic_example
using OrdinaryDiffEqDefault, Plots
u0 = [:X => 0.5]
ps = [:d => 2.0, pIn => spline]
oprob = ODEProblem(bd_model, u0, tend, ps)
sol = solve(oprob)
plot(sol)
```
!!! note
    For this simple example, $(2 + t)/(1 + t)$ could have been used directly as a reaction rate, technically making the functional parameter approach unnecessary.

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
We are now ready to declare our model. We will consider a protein with an active and an inactive form ($Pₐ$ and $Pᵢ$) where the activation is driven by the presence of sunlight. In this example we we create our model using the [programmatic approach](@ref programmatic_CRN_construction). Do note the special syntax we use to declare our input parameter, where we both designate it as a generic function and its type as the type of our interpolated input. Also note that, within the model, we mark the input parameter (`light_in`) a function of `t`.
```@example functional_parameters_circ_rhythm
using Catalyst
t = default_t()
in_type = typeof(interpolated_light)
@parameters kA kD (light_in::in_type)(..)
@species Pₐ(t) Pᵢ(t)
rxs = [
    Reaction(light_in(t)*kA, [Pᵢ], [Pₐ]),
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
It is possible to use time-dependent inputs when creating models [through the DSL](@ref dsl_description) as well. However, it can still be convenient to declare the input parameter programmatically as above. Next, we form an expression of it as a function of time, and then [interpolate](@ref dsl_advanced_options_symbolics_and_DSL_interpolation) it into our DSL-declaration:
```@example functional_parameters_circ_rhythm
input = light_in(t)
rs_dsl = @reaction_network rs begin
    ($input*kA, kD), Pᵢ <--> Pₐ
end
```
We can confirm that this model is identical to our programmatic one (and should we wish to, we can simulate it using identical syntax syntax).
```@example functional_parameters_circ_rhythm
rs == rs_dsl
```

## [Non-time functional parameters](@id functional_parameters_sir)
Previously we have demonstrated functional parameters that are a function of time. However, functional parameters can be functions of any variable (however, currently, more than one argument is not supported). Here we will demonstrate this using a [SIR model](@ref basic_CRN_library_sir), but instead of having the infection rate scale linearly with the number of infected individuals, we instead assume we have measured data of the infection rate (as dependent on the number of infected individuals) and wish to use this instead. Normally we use the following infection reaction in the SIR model:
```julia
@reaction k1, S + I --> 2I
```
In practise, this is identical to
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
I_meassured = 300.0 *(0.8*rand(length(I_grid)) .+ 0.6) .* I_grid ./ (300 .+ I_grid)
I_rate = LinearInterpolation(I_meassured, I_grid)
plot(I_rate; label = "Meassured infection rate")
plot!(I_grid, I_grid; label = "Normal SIR infection rate")
```
Next, we create our model (using the DSL-approach). As `I_rate` will be a function of $I$ we will need to declare this species first as well.
```@example functional_parameters_sir
using Catalyst
@species I(default_t())
inf_rate = I_rate(I)
sir = @reaction_network rs begin
    k1*$inf_rate, S --> I
    k2, I --> R
end
nothing # hide
```
Finally, we can simualte our model.
```@example functional_parameters_sir
using OrdinaryDiffEqDefault
u0 = [:S => 99.0, :I => 1.0, :R => 0.0]
ps = [:k1 => 0.002, :k2 => 0.01, inf_rate => I_rate]
oprob = ODEProblem(sir, u0, 250.0, ps)
sol = solve(oprob)
plot(sol)
```
!!! note
    When declaring a functional parameter of time, it is easy to set its grid values (i.e. ensure they range from the first to the final time point). For Functional parameters that depend on species concentrations it is trickier, and one must make sure that any potential input-species values that can occur during the simulation are represented in the interpolation. In this example it is fairly straightfor

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
Finally, DataInterpolation.jl also allows various [extrapolation methods](https://docs.sciml.ai/DataInterpolations/stable/extrapolation_methods/).