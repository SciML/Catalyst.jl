# [The Catalyst DSL - Advanced Features and Options](@id dsl_advanced_options)

Within the Catalyst DSL, each line can represent either *a reaction* or *an option*. The [previous tutorial](@ref dsl_description) described how to create reactions. This one will instead describe how to use options. These are typically used to supplant the model with additional information. Examples include the declaration of initial condition/parameter default values, or the creation observables or events. 

All options are designated begins with a symbol starting with `@`, followed by its input. E.g. the `@observables` options allows for the generation of observables. Each option can only be used once within each declaration of `@reaction_network`. A full list of options can be found [here], with most (but not all) being described in more detail below.

 This tutorial will also describe some additional advanced DSL features that does not include using an option.

## [Explicit specification of network species and parameters](@id dsl_advanced_options_declaring_species_and_parameters)
Previously, we mentioned that 

### [Setting default values for species and parameters](@id dsl_advanced_options_default_vals)

### [Setting parametric initial conditions](@id dsl_advanced_options_parametric_initial_conditions)

### [Specifying non-species variables](@id dsl_advanced_options_variables)

### [Designating metadata for species and parameters](@id dsl_advanced_options_species_and_parameters_metadata)

## [Setting reaction metadata](@id dsl_advanced_options_)

## [Naming reaction networks](@id dsl_advanced_options_)

## [Creating observables](@id dsl_advanced_options_)

## [Creating events](@id dsl_advanced_options_)

## [Specifying non-time independent variables](@id dsl_advanced_options_)


## Explicit specification of network species and parameters
Recall that the `@reaction_network` macro automatically designates symbols used
in the macro as either parameters or species, with symbols that appear as a
substrate or product being species, and all other symbols becoming parameters
(i.e. those that only appear within a rate expression and/or as [stoichiometric coefficients](@ref parametric_stoichiometry)). Sometimes, one might want to manually override this default
behavior for a given symbol. E.g one might want something to be considered as a
species, even if it only appears within a rate expression. In the following
network
```@example dsl_1
rn = @reaction_network begin
  k*X, Y --> 0
end
```
`X` (as well as `k`) will be considered a parameter.

By using the `@species` and `@parameters` options within the `@reaction_network`
macro, one can manually declare that specified symbols should be considered a
species or parameter. E.g in:
```@example dsl_1
rn = @reaction_network begin
  @species X(t) Y(t)
  k*X, Y --> 0
end
```
`X` and `Y` are set as species. Please note that when declaring species using
the `@species` option, their dependant variable (almost always `t`) also needs
to be designated. Similarly in
```@example dsl_1
rn = @reaction_network begin
  @parameters k
  k*X, Y --> 0
end
```
both `X` and `k` will be considered as parameters. It is also possible to use
both options simultaneously, allowing users to fully specify which symbols are
species and/or parameters:
```@example dsl_1
rn = @reaction_network begin
  @species X(t) Y(t)
  @parameters k
  k*X, Y --> 0
end
```
Here, `X` and `Y` are designated as species and `k` as a parameter.

The lists provided to the `@species` and `@parameters` options do not need to be extensive. Any symbol that appears in neither list will use the default option as determined by the macro. E.g. in the previous example, where we only want to change the default designation of `X` (making it a species rather than a parameter), we can simply write:
```@example dsl_1
rn = @reaction_network begin
  @species X(t)
  k*X, Y --> 0
end
```

Finally, note that the `@species` and `@parameters` options can also be used in
`begin ... end` block form, allowing more formatted lists of species/parameters:
```@example dsl_1
rn = @reaction_network begin
  @parameters begin
      d1
      d2
  end
  @species begin
      X1(t)
      X2(t)
  end
  d2, X2 --> 0
  d1, X1 --> 0
end
```
This can be especially useful when declaring default values for clarity of model
specification (see the next section).

## [Setting default values for initial conditions and parameters](@id dsl_description_defaults)
When using the `@species` and ` @parameters` macros to declare species and/or
parameters, one can also provide default initial conditions for each species and
values for each parameter:
```@example dsl_1
rn = @reaction_network begin
  @species X(t)=1.0
  @parameters p=1.0 d=0.1
  p, 0 --> X
  d, X --> ∅
end
```
This system can now be simulated without providing initial condition or
parameter vectors to the DifferentialEquations.jl solvers:
```@example dsl_1
using DifferentialEquations, Plots
u0 = []
tspan = (0.0, 10.0)
p = []
oprob = ODEProblem(rn, u0, tspan, p)
sol = solve(oprob)
plot(sol)
```

When providing default values, it is possible to do so for only a subset of the
species or parameters, in which case the rest can be specified when constructing
the problem type to solve:
```@example dsl_1
rn = @reaction_network begin
  @species X(t)
  @parameters p=1.0 d
  p, 0 --> X
  d, X --> 0
end

u0 = [:X => 1.0]
tspan = (0.0, 10.0)
p = [:d => .1]
oprob = ODEProblem(rn, u0, tspan, p)
sol = solve(oprob)
plot(sol)
```

Finally, default values can be overridden by passing mapping vectors to the
DifferentialEquations.jl problem being constructed. Only those initial conditions
or parameters for which we want to change their value from the default will need to be passed
```@example dsl_1
u0 = [:X => 1.0]
tspan = (0.0, 10.0)
p = [:p => 2.0, :d => .1]   # we change p to 2.0
oprob = ODEProblem(rn, u0, tspan, p)
sol = solve(oprob)
plot(sol)
```

## [Setting initial conditions that depend on parameters](@id dsl_description_parametric_initial_conditions)
It is possible to set the initial condition of one (or several) species so that they depend on some system parameter. This is done in a similar way as default initial conditions, but giving the parameter instead of a value. When doing this, we also need to ensure that the initial condition parameter is a variable of the system:
```@example dsl_1
rn = @reaction_network begin
  @parameters X0
  @species X(t)=X0
  p, 0 --> X
  d, X --> ∅
end
```
We can now simulate the network without providing any initial conditions:
```@example dsl_1
u0 = []
tspan = (0.0, 10.0)
p = [:p => 2.0, :d => .1, :X0 => 1.0]
oprob = ODEProblem(rn, u0, tspan, p)
sol = solve(oprob)
plot(sol)
```

## Naming the generated `ReactionSystem`
ModelingToolkit uses system names to allow for compositional and hierarchical
models. To specify a name for the generated `ReactionSystem` via the
[`@reaction_network`](@ref) macro, just place the name before `begin`:
```@example dsl_1
rn = @reaction_network production_degradation begin
  p, ∅ --> X
  d, X --> ∅
end
ModelingToolkit.nameof(rn) == :production_degradation
```

## Including non-species variables
Non-species state variables can be specified in the DSL using the `@variables`
macro. These are declared similarly to species. For example,
```@example dsl_1
rn_with_volume = @reaction_network begin
  @variables V(t)
  k*V, 0 --> A
end
```
creates a network with one species
```@example dsl_1
species(rn_with_volume)
```
and one non-species
```@example dsl_1
nonspecies(rn_with_volume)
```
giving two state variables, always internally ordered by species and then
nonspecies:
```@example dsl_1
states(rn_with_volume)
```

`rn_with_volume` could then be extended with constraint equations for how `V(t)`
evolves in time, see the [associated tutorial](@ref constraint_equations).

## Specifying alternative time variables and/or extra independent variables
While the DSL defaults to allowing `t` as the time variable, one can use the
`@ivs` macro to specify an alternative independent variable. For example, to
make `s` the default time variable one can say
```@example dsl_1
rn_with_s = @reaction_network begin
    @ivs s
    @variables V(s)
    @species B(s)
    k, A + V*B --> C
end
show(stdout, MIME"text/plain"(), rn_with_s)  # hide
```
where we see all states are now functions of `s`.

Similarly, if one wants states to be functions of more than one independent
variable, for example to encode a spatial problem, one can list more than one
variable, i.e. `@ivs t x y`. Here the first listed independent variable is
always chosen to represent time. For example,
```@example dsl_1
rn_with_many_ivs = @reaction_network begin
    @ivs s x
    @variables V1(s) V2(s,x)
    @species A(s) B(s,x)
    k, V1*A --> V2*B + C
end
show(stdout, MIME"text/plain"(), rn_with_many_ivs)  # hide
```
Here again `s` will be the time variable, and any inferred species, `C` in this
case, are made functions of both variables, i.e. `C(s, x)`.

## [Interpolation of Julia variables](@id dsl_description_interpolation_of_variables)
The DSL allows Julia variables to be interpolated for the network name, within
rate constant expressions, or for species/stoichiometry within reactions. Using
the lower-level symbolic interface we can then define symbolic variables and
parameters outside of the macro, which can then be used within expressions in
the DSL (see the [Programmatic Construction of Symbolic Reaction Systems](@ref programmatic_CRN_construction)
tutorial for details on the lower-level symbolic interface). For example,
```@example dsl_1
@parameters k α
@variables t
@species A(t)
spec = A
par = α
rate = k*A
name = :network
rn = @reaction_network $name begin
    $rate*B, 2*$spec + $par*B --> $spec + C
  end
```
As the parameters `k` and `α` were pre-defined and appeared via interpolation,
we did not need to declare them within the `@reaction_network` macro,
i.e. they are automatically detected as parameters:
```@example dsl_1
parameters(rn)
```
as are the species coming from interpolated variables
```@example dsl_1
species(rn)
```

!!! note
    When using interpolation, expressions like `2$spec` won't work; the
    multiplication symbol must be explicitly included like `2*$spec`.
