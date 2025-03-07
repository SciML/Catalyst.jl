# [Inputs and time-dependent (or functional) parameters](@id time_dependent_parameters)

## [Basic example](@id time_dependent_parameters_basic_example)
Let us first consider an easy, quick-start example. We will consider a simple [birth-death model](@ref basic_CRN_library_bd), but where the birth rate is determined by an input parameter (which value depends on time). First, we [define the input parameter programmatically](@ref programmatic_CRN_construction), and its values across all time values using the [DataInterpolations.jl](https://github.com/SciML/DataInterpolations.jl) package. In this example we will use the input function $pIn(t) = (2 + t)/(1 + t)$.
```@example time_dependent_parameters_basic_example
using Catalyst, DataInterpolations
t = default_t()
tend = 10.0
ts = collect(0.0:0.01:tend)
spline = LinearInterpolation((2 .+ ts) ./ (1 .+ ts), ts)
@parameters (pIn::typeof(spline))(..)
nothing # hide
```
Next, we create our model, [interpolating](@ref dsl_advanced_options_symbolics_and_DSL_interpolation) the input parameter into it (making it a function of `t`).
```@example time_dependent_parameters_basic_example
input = pIn(default_t())
bd_model = @reaction_network begin
    $input, 0 --> X
    d, X --> 0
end
```
Finally, we can simulate our model as normal, but setting the value of the `pIn` parameter to our interpolated data.
```@example time_dependent_parameters_basic_example
using OrdinaryDiffEq, Plots
u0 = [:X => 0.5]
ps = [:d => 2.0, pIn => spline]
oprob = ODEProblem(bd_model, u0, tend, ps)
sol = solve(oprob)
plot(sol)
```

## [Inserting a customised, time-dependent, input] (@id time_dependent_parameters_circ_rhythm)
Let us now go through everything again, but providing some more details. Let us first consider the input parameter. We have previously described how a [time-dependent rate can model a circadian rhythm](@ef dsl_description_nonconstant_rates_time). For real applications, due to e.g. clouds, sunlight is not a perfect sine wave. Here, a common solution is to take real sunlight data from some location and us in ones model. We will instead create some synthetic noisy data as our light input:
```@example time_dependent_parameters_circ_rhythm
using Plots
tend = 120.0
ts = collect(0.0:1.0:tend)
light = sin.(ts/6) .+ 1
light = [max(0.0, l - rand()) for l in light]
plot(ts, light; seriestype = :scatter, label = "Experienced light")
```
Now this input is only actually be defined at the sample points, making it incompatible with a continuous ODE simulation. To enable this, we will use the DataInterpolations package to create an interpolated version of this data, which forms the actual input:
```@example time_dependent_parameters_circ_rhythm
using DataInterpolations
interpolated_light = LinearInterpolation(light, ts)
plot(interpolated_light)
```
Finally, we are ready to declare our model. We will consider a protein with an active and an inactive from ($Pa$ and $Pi$) where the activation is driven by the presence of sunlight. In this example we we create our model using the [programmatic approach](@ref programmatic_CRN_construction). Do note the special syntax we use to declare our input parameter, where we designate its type as the type of our interpolated input. Also note that, within the model, we mark the input parameter (`light_in`) a function of `t`.
```@example time_dependent_parameters_circ_rhythm
using Catalyst
t = defaul_t()
in_type = typeof(interpolated_light)
@parameters kA kD (light_in::in_type)(..)
@species Pa(t) Pi(t)
rxs = [
    Reaction(kA*light_in(t), [Pi], [Pa]),
    Reaction(kD, [Pa], [Pi])
]
@named rs = ReactionSystem(rxs, t)
rs = complete(rs)
```
Now we can simulate our model our. Here, we use the interpolated data as the input parameter;s value.
```@example time_dependent_parameters_circ_rhythm
using OrdinaryDiffEqDefault
u0 = [Pi => 1.0, Pa => 0.0]
ps = [kA => 1.5, kB => 1.0, pIn => interpolated_light]
oprob = ODEProblem(rs, u0, tend, ps)
sol = solve()
```

### [Interpolating the input into the DSL] (@id time_dependent_parameters_circ_rhythm_dsl)
It is possible to use time-dependent inputs when creating models [through the DSL](@ref dsl_description) as well. However, it can still be convenient to declare the input parameter programmatically as above. Next, we form an expression of it as a function of time, and then [interpolate](@ref https://github.com/SciML/DataInterpolations.jl) it into our DSL-declaration:
```@example time_dependent_parameters_circ_rhythm
input = pIn(t)
rs_dsl = @reaction_network rs begin
    (kA*$input, kI), Pi <--> Pa
end
```
When can confirm that this model is identical to our programmatic one (and should we wish to, we can simulate it using identical syntax syntax).
```@example time_dependent_parameters_circ_rhythm
rs == rs_dsl
```
