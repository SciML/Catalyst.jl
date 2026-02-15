# [Coupled non-reaction network equations to models](@id coupled_models)
```julia
using Pkg
Pkg.activate(".")
Pkg.add("Catalyst")
Pkg.add("OrdinaryDiffEqTsit5")
Pkg.add("Plots")
```
```@raw html
</details>
```
Non-reaction model components can be inserted directly in a Catalyst model. Here we will briefly describe the simplest case: adding an ODE to a model declared through the `@reaction_network` DSL. The equation is added using the `@equation` option, after which the equation is written (with `D(V)` denoting differential with respect to time).
```julia
using Catalyst, OrdinaryDiffEqDefault, Plots

# Create model with a ODE for the variable `V` (e.g. denoting Volume).
rn = @reaction_network begin
    @parameters Vmax # Parameters appearing in equations must be explicitly declared as such.
    @equations D(V) ~ X * (Vmax - V) # `~` denotes equality (=) in an equation.
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