# [Coupled non-reaction network equations to models](@id coupled_models)
```@raw html
<details><summary><strong>Environment setup and package installation</strong></summary>
```
The following code sets up an environment for running the code on this page.
```julia
using Pkg
Pkg.activate(; temp = true) # Creates a temporary environment, which is deleted when the Julia session ends.
Pkg.add("Catalyst")
Pkg.add("OrdinaryDiffEqTsit5")
Pkg.add("Plots")
```
```@raw html
</details>
```
  \

Non-reaction model components can be inserted directly in a Catalyst model. Here we will briefly describe the simplest case: adding an ODE to a model declared through the `@reaction_network` DSL. The equation is added using the `@equations` option, after which the equation is written (with `D(V)` denoting differential with respect to time).
```julia
using Catalyst, OrdinaryDiffEqDefault, Plots

# Create model with a ODE for the variable `V` (e.g. denoting Volume).
rn = @reaction_network begin
    @parameters Vmax # Parameters appearing in equations must be explicitly declared as such.
    @equations D(V) ~ X * (Vmax - V) # `~` (not `=`!) separates equation left and right-hand sides.
    (p,d), 0 <--> X
end

# The model can be simulated using normal syntax.
u0 = [:X => 2.0, :V => 0.5] # Must give an initial condition for `V` as well as `X`.
ps = [:p => 2.0, :d => 0.25, :Vmax => 2.0]
oprob = ODEProblem(rn, u0, 10.0, ps)
sol = solve(oprob)
plot(sol)
```
Other types of model components (algebraic equations and brownian and poissonian processes) are also possible, as are non-DSL approaches for adding these to a model.
```@raw html
</details>
```
  \

In many applications one has additional algebraic or differential equations for
non-chemical species that can be coupled to a chemical reaction network model.
Catalyst supports coupled differential and algebraic equations, and currently
allows conversion of such coupled systems to ModelingToolkitBase ODE, SDE, and
nonlinear `System`s.

In this tutorial we'll illustrate how to make use of coupled (i.e.
ODE/algebraic) equations. Let's consider a model of a cell with
volume $V(t)$ that grows at a rate $\lambda$. For now we'll assume the cell can
grow indefinitely. We'll also keep track of one protein $P(t)$, which is
produced at a rate proportional $V$ and can be degraded.

## [Coupling ODE via the DSL](@id coupled_models_dsl)
The easiest way to include ODEs and algebraic equations is to just include them
when using the DSL to specify a model. Here we include an ODE for $V(t)$ along
with degradation and production reactions for $P(t)$:
```@example ceq1
using Catalyst, OrdinaryDiffEqTsit5, Plots

rn = @reaction_network growing_cell begin
    # the growth rate
    @parameters λ = 1.0

    # assume there is no protein initially
    @species P(t) = 0.0

    # set the initial volume to 1.0
    @variables V(t) = 1.0

    # the reactions
    V, 0 --> P
    1.0, P --> 0

    # the coupled ODE for V(t)
    @equations begin
        D(V) ~ λ * V
    end
end
```
We can now create an `ODEProblem` from our model and solve it to see how $V(t)$
and $P(t)$ evolve in time:
```@example ceq1
oprob = ODEProblem(rn, [], (0.0, 1.0))
sol = solve(oprob, Tsit5())
plot(sol)
```

## Coupling ODEs via directly building a `ReactionSystem`
As an alternative to the previous approach, we could have also constructed our
`ReactionSystem` all at once using the symbolic interface:
```@example ceq2
using Catalyst, OrdinaryDiffEqTsit5, Plots

t = default_t()
D = default_time_deriv()

@parameters λ = 1.0
@variables V(t) = 1.0
@species P(t) = 0.0
eq = D(V) ~ λ * V
rx1 = @reaction $V, 0 --> $P
rx2 = @reaction 1.0, $P --> 0
@named growing_cell = ReactionSystem([rx1, rx2, eq], t)
growing_cell = complete(growing_cell)

oprob = ODEProblem(growing_cell, [], (0.0, 1.0))
sol = solve(oprob, Tsit5())
plot(sol)
```


## [Coupling ODEs via extending a system](@id coupled_models_extending_systems)

Finally, we could also construct our model by using compositional modeling. Here
we create separate `ReactionSystem`s, one with the reaction part of the model, and one with
the differential equation part. Next, we extend the first `ReactionSystem` with the second one. Let's
begin by creating these two systems.

Here, to create differentials with respect to time (for our differential
equations), we must import the time differential operator from Catalyst. We do
this through `D = default_time_deriv()`. Here, `D(V)` denotes the differential
of the variable `V` with respect to time.

```@example ceq2b
using Catalyst, OrdinaryDiffEqTsit5, Plots

t = default_t()
D = default_time_deriv()

# set the growth rate to 1.0
@parameters λ = 1.0

# set the initial volume to 1.0
@variables V(t) = 1.0

# build the ReactionSystem for dV/dt
eq = [D(V) ~ λ * V]
@named osys = ReactionSystem(eq, t)

# build the ReactionSystem with no protein initially
rn = @network_component begin
    @species P(t) = 0.0
    $V,   0 --> P
    1.0, P --> 0
end
```
Notice, here we interpolated the variable `V` with `$V` to ensure we use the
same symbolic unknown variable in the `rn` as we used in building `osys`. See
the doc section on [interpolation of variables](@ref
dsl_advanced_options_symbolics_and_DSL_interpolation) for more information. We
also use `@network_component` instead of `@reaction_network` as when merging
systems together Catalyst requires that the systems have not been marked as
`complete` (which indicates to Catalyst that a system is finalized).

We can now merge the two systems into one complete `ReactionSystem` model using
[`ModelingToolkitBase.extend`](@ref):
```@example ceq2b
@named growing_cell = extend(osys, rn)
growing_cell = complete(growing_cell)
```

We see that the combined model now has both the reactions and ODEs as its
equations. To solve and plot the model we proceed like normal
```@example ceq2b
oprob = ODEProblem(growing_cell, [], (0.0, 1.0))
sol = solve(oprob, Tsit5())
plot(sol)
```

## [Coupled reaction networks with differential equations](@id coupled_models_diffeqs)
When declaring models using the DSL, differential equations can be inserted directly into the model through the `@equations` option. I.e. in the following model we we have a protein in an inactive and active form ($Xᵢ$ and $Xₐ$, respectively). The conversion to the active form is promoted by the presence of a nutrient ($N$). If the nutrient is then degraded linearly with the concentration of $Xₐ$, then this can be describe by the differential equation
```math
\mathbf{\frac{dN}{dt} = - d Xₐ(t) N(t)}.
```
This is declared with the following model
```@example coupled_diff_eqs
using Catalyst
rs = @reaction_network begin
    @parameters d
    @equations D(N) ~ -d*Xₐ*N
    (kₐ*N, kᵢ), Xᵢ <--> Xₐ
end
```
Here we note that:
- Symbols not occurring as species within the `@equations` option are by default assumed to be [variables](@ref coupled_models_variables). Hence, `d` must [explicitly be declared as a parameter](@ref dsl_advanced_options_declaring_species_and_parameters).
- Whenever the `D(X)` notation occurs within the `@equations` option, this is by default assumed to be the differential with respect to the default time independent variable (i.e. $dx/dt$). This mean you have to *take extra care* if using `D` as e.g. a species name (typically by [declaring a custom symbol as a differential operator](@ref faq_custom_differentials)).

The model can be simulated and plotted using normal syntax, providing an initial condition for `N` and a value for `d` through the normal initial condition and parameter value vectors.
```@example coupled_diff_eqs
using OrdinaryDiffEqDefault, Plots
u0 = [:Xᵢ => 1.0, :Xₐ => 0.0, :N => 1.0]
ps = [:kₐ => 2.0, :kᵢ => 0.5, :d => 1.0]
oprob = ODEProblem(rs, u0, 10.0, ps)
sol = solve(oprob)
plot(sol)
```

!!! note
    Species that occur within reactions cannot be subject to a differential operator within a `@equations` block.
!!! note
    When providing multiple equations to a model, write each equation on a separate line within a `@equations begin ... end` block.
!!! note
    Generally, the rules for declaring equations are the same as those used within the [ModelingToolkitBase.jl](https://github.com/SciML/ModelingToolkit.jl) package.

## [Coupled reaction networks with algebraic equations](@id coupled_models_algeqs)
Catalyst also permits the coupling*algebraic equations* (equations not containing differentials). In practice, these are handled similarly to how differential equations are handled, but with the following differences:
- The `structural_simplify = true` option must be provided to any `ODEProblem` (or other problem) generated from a model containing algebraic equations.
- For variables involved in algebraic equations, a *guess* is provided, rather than an initial condition.

A special case of algebraic equations (where a new variable is trivially produced by an algebraic expression) is *observables*. These are described in a separate section [here](@ref dsl_advanced_options_observables).

We will demonstrate using a system a species $X$ is produced and degraded, and can also bind/dissociate to form a dime $X2$. Here, if the binding/unbinding dynamics is much faster than the production/degradation dynamics, $X2$ can be eliminated through an algebraic equation. This creates the following model:
```@example coupled_eqs_alg_eq
using Catalyst
algebraic_crn = @reaction_network begin
    @parameters kB kD
    @equations kB*X^2 ~ kD*X2
    (mmr(X2,v,K),d), 0 <--> X
end
```
where we also made the production rate of $X$ depend on $X2$ through a Michaelis-Menten function.

We can now simulate the model. Here, at each time step, $X2$'s value can be computed from $X$'s. Simialrily for the initial condition, if we know $X$, we can compute $X2$. Hence, we will not provide an initial condition value for $X2$. However, we will have to provide an initial *guess* for $X2$'s value, from which an internal [nonlinear solve-call](@ref steady_state_solving_nonlinear) will be used to compute the full ODE simulation initial condition.
```@example coupled_eqs_alg_eq
u0 = [:X => 1.0]
ps = [:v => 2.0, :K => 1.0, :d => 0.2, :kB => 0.1, :kD => 0.4]
guesses = [:X2 => 1.0]
nothing # hide
```
Next, we provide `guesses` to our `ODEProblem` as an additional argument. Furthermore, we will use the `mtkcompile = true` argument, which is always required when simulating models containing algebraic equations. With these modifcations, the model can be simulated using standard workflows.
```@example coupled_eqs_alg_eq
using OrdinaryDiffEqDefault, Plots
oprob = ODEProblem(algebraic_crn, u0, 10.0, ps; structural_simplify = true, guesses)
sol = solve(oprob)
plot(sol)
```

There is no requirement on which values are provided as guesses and which as initial conditions. I.e. if we know the value of $X2$ we can provide this as the initial condition while instead providing a guess for $X$'s value.
```@example coupled_eqs_alg_eq
u0 = [:X2 => 1.0]
guesses = [:X => 1.0]
oprob = ODEProblem(algebraic_crn, u0, 10.0, ps; structural_simplify = true, guesses)
nothing # hide
```

There is one exception to this, which is if the algebraic equation is formatted such that a variable is isolated on the equation left-hand side. I.e. if we were to declare the model using
```@example coupled_eqs_alg_eq
algebraic_crn_alt = @reaction_network begin
    @parameters kB kD
    @equations X2 ~ kB*X^2/kD
    (mmr(X2,v,K),d), 0 <--> X
end
```
Here, $X2$ can be fully eliminated from the system when the equations are generated, and we neither eed to provide an initial condition nor a guess for it:
```@example coupled_eqs_alg_eq
u0 = [:X => 1.0]
ps = [:v => 2.0, :K => 1.0, :d => 0.2, :kB => 0.1, :kD => 0.4]
oprob = ODEProblem(algebraic_crn_alt, u0, 10.0, ps; structural_simplify = true, guesses)
sol = solve(oprob)
nothing # hide
```
A side effect of this is that $X2$ by default is not plotted wit the solution, something which must be explicitly requested using the [`idxs` argument](@ref simulation_plotting_options):
```@example coupled_eqs_alg_eq
plot(sol; idxs = [:X, :X2])
```
In practise, we can check the equations that are generated by manually converting them to ODEs through `ode_model`:
```@example coupled_eqs_alg_eq
ode_model(algebraic_crn)
```
```@example coupled_eqs_alg_eq
ode_model(algebraic_crn_alt)
```

## [Notes on *species* and *variables*](@id coupled_models_variables)
Throughout this tutorial we have distinguishing between *species* and *variables*. Functionality, both these are *unknown quantities* that we are simulating or solving for in `solve` calls.  The differences are that
- Only *species* can participate are reactants in reactions.
- Only *variables* can be the subject of differentials in equations.

Furthermore, species and variables as declared using different syntax. Either when declaring them [programmatically](@ref programmatic_CRN_construction), or through DSL options, the use `@species` and `@variables`, respectively. I.e. if we want to designate a [default initial condition value](@ref dsl_advanced_options_default_vals) for $N$ in the [nutrient model](@ref coupled_models_diffeqs) we use
```@example coupled_eqs_variables
using Catalyst # hide
rs = @reaction_network begin
    @parameters d
    @variables N = 1.0
    @equations D(N) ~ -d*Xₐ*N
    (kₐ*N, kᵢ), Xᵢ <--> Xₐ
end
nothing # hide
```

To retrieve a either a models species, variables, or all unknowns, use the `species`, `nonspecies`, and `unknowns` functions:
```@example coupled_eqs_variables
species(rs)
```
```@example coupled_eqs_variables
nonspecies(rs)
```
```@example coupled_eqs_variables
unknowns(rs)
```
