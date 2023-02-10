# FAQs

## How to index solution objects using symbolic variables and observables?
One can directly use symbolic variables to index into SciML solution objects.
Moreover, observables can also be evaluated in this way. For example,
consider the system
```@example faq1
using Catalyst, DifferentialEquations, Plots
rn = @reaction_network ABtoC begin
  (k₊,k₋), A + B <--> C
end

# initial condition and parameter values
setdefaults!(rn, [:A => 1.0, :B => 2.0, :C => 0.0, :k₊ => 1.0, :k₋ => 1.0])
```
Let's convert it to a system of ODEs, using the conservation laws of the system
to eliminate two of the species:
```@example faq1
osys = convert(ODESystem, rn; remove_conserved=true)
show(osys) # hide
```
Notice the resulting ODE system has just one ODE
```@example faq1
equations(osys)
show(equations(osys)) # hide
```
while algebraic observables have been added for the two removed species (in
terms of the conservation law constants, `_ConLaw[1]` and `_ConLaw[2]`)
```@example faq1
observed(osys)
show(observed(osys))  # hide
```
Let's solve the system and see how to index the solution using our symbolic
variables
```@example faq1
oprob = ODEProblem(osys, [], (0.0, 10.0), [])
sol = solve(oprob, Tsit5())
```
Suppose we want to plot just species `C`, without having to know its integer
index in the state vector. We can do this using the symbolic variable `C`, which
we can get at in several ways
```@example faq1
sol[osys.C]
```
or
```@example faq1
@unpack C = osys
sol[C]
```
To evaluate `C` at specific times and plot it we can just do
```@example faq1
t = range(0.0, 10.0, length=101)
plot(t, sol(t, idxs = C), label = "C(t)", xlabel = "t")
```
If we want to get multiple variables we can just do
```@example faq1
@unpack A, B = osys
sol(t, idxs = [A, B])
```
Plotting multiple variables using the SciML plot recipe can be achieved
like
```@example faq1
plot(sol; idxs = [A, B])
```

## How to disable rescaling of reaction rates in rate laws?
As explained in the [Reaction rate laws used in simulations](@ref) section, for
a reaction such as `k, 2X --> 0`, the generated rate law will rescale the rate
constant, giving `k*X^2/2` instead of `k*X^2` for ODEs and `k*X*(X-1)/2` instead
of `k*X*(X-1)` for jumps. This can be disabled when directly `convert`ing a
[`ReactionSystem`](@ref). If `rn` is a generated [`ReactionSystem`](@ref), we can
do
```julia
osys = convert(ODESystem, rn; combinatoric_ratelaws=false)
```
Disabling these rescalings should work for all conversions of `ReactionSystem`s
to other `ModelingToolkit.AbstractSystem`s.

## How to use non-integer stoichiometric coefficients?
```julia
rn = @reaction_network begin
  k, 2.5*A --> 3*B
end
```
or directly via
```julia
@parameters k b
@variables t A(t) B(t) C(t) D(t)
rx1 = Reaction(k,[B,C],[B,D], [2.5,1],[3.5, 2.5])
rx2 = Reaction(2*k, [B], [D], [1], [2.5])
rx3 = Reaction(2*k, [B], [D], [2.5], [2])
@named mixedsys = ReactionSystem([rx1,rx2,rx3],t,[A,B,C,D],[k,b])
osys = convert(ODESystem, mixedsys; combinatoric_ratelaws=false)
```
Note, when using `convert(ODESystem, mixedsys; combinatoric_ratelaws=false)` the
`combinatoric_ratelaws=false` parameter must be passed. This is also true when
calling `ODEProblem(mixedsys,...; combinatoric_ratelaws=false)`. As described
above, this disables Catalyst's standard rescaling of reaction rates when
generating reaction rate laws, see also the [Reaction rate laws used in
simulations](@ref) section. Leaving this keyword out for systems with floating
point stoichiometry will give an error message.

For a more extensive documentation of non-integer stoichiometric coefficients, pelase read the [Symbolic Stochiometries](@ref parametric_stoichiometry) section.

## How to set default values for initial conditions and parameters?
How to set defaults when using the `@reaction_network` macro is desribed in detail [here](@ref dsl_description_defaults). Briefly, to do this, use the `@species` and
`@parameters` options:
```julia
using Catalyst
sir = @reaction_network sir begin
    @species S(t)=999.0 I(t)=1.0 R(t)=0.0
    @parameters β=1e-4 ν=0.01
    β, S + I --> 2I
    ν, I --> R
end
```

When directly constructing a `ReactionSystem`, we can set the symbolic values to
have the desired default values, and this will automatically be propagated
through to the equation solvers:
```julia
Plots, OrdinaryDiffEq
@parameters β=1e-4 ν=.01
@variables t S(t)=999.0 I(t)=1.0 R(t)=0.0
rx1 = Reaction(β, [S,I], [I], [1,1], [2])
rx2 = Reaction(ν, [I], [R])
@named sir = ReactionSystem([rx1, rx2], t)
oprob = ODEProblem(sir, [], (0.0, 250.0))
sol = solve(oprob, Tsit5())
plot(sol)
```

One can also build a mapping from parameter/species to value/initial
condition and pass these to the `ReactionSystem` via the `defaults` keyword
argument:
```julia
@parameters β ν
@variables t S(t) I(t) R(t)
rx1 = Reaction(β, [S,I], [I], [1,1], [2])
rx2 = Reaction(ν, [I], [R])
defs = [β => 1e-4, ν => .01, S => 999.0, I => 1.0, R => 0.0]
@named sir = ReactionSystem([rx1, rx2], t; defaults = defs)
```

Finally, default values can also be added after creating the system via the
`setdefaults!` command, like
```julia
sir = @reaction_network sir begin
    β, S + I --> 2I
    ν, I --> R
end
setdefaults!(sir, [:β => 1e-4, :ν => .01, :S => 999.0, :I => 1.0, :R => 0.0])
```

## How to specify initial conditions and parameters values for `ODEProblem` and other problem types?
To explicitly pass initial conditions and parameters we can use mappings from
Julia `Symbol`s corresponding to each variable/parameter to values, or from
ModelingToolkit symbolic variables to each variable/parameter. Using `Symbol`s
we have
```julia
rn = @reaction_network begin
    α, S + I --> 2I
    β, I --> R
end
u0 = [:S => 999.0, :I => 1.0, :R => 0.0]
p  = (:α => 1e-4, :β => .01)
op  = ODEProblem(rn, u0, (0.0,250.0), p)
sol = solve(op, Tsit5())
```
while using ModelingToolkit symbolic variables we have
```julia
@parameters α β
@variables t S(t) I(t) R(t)
u0 = [S => 999.0, I => 1.0, R => 0.0]
p  = (α => 1e-4, β => .01)
op  = ODEProblem(rn, u0, (0.0,250.0), p)
sol = solve(op, Tsit5())
```

## How to modify generated ODEs?
Conversion to other `ModelingToolkit.AbstractSystem`s allows the possibility to
modify the system with further terms that are difficult to encode as a chemical
reaction. For example, suppose we wish to add a forcing term, $10\sin(10t)$, to
the ODE for `dX/dt` above. We can do so as:
```julia
dXdteq = equations(osys)[1]
t      = get_iv(osys)
dXdteq = Equation(dXdteq.lhs, dXdteq.rhs + 10*sin(10*t))
@named osys2  = ODESystem([dXdteq], t, states(osys), parameters(osys))
oprob  = ODEProblem(osys2, u0map, tspan, pmap)
osol   = solve(oprob, Tsit5())
```
We can add $e^{-X}$ to $dX/dt$ as a forcing term by
```julia
dXdteq = equations(osys)[1]
@variables X(t)
dXdteq = Equation(dXdteq.lhs, dXdteq.rhs + exp(-X))
@named osys2  = ODESystem([dXdteq], t, states(osys), parameters(osys))
oprob  = ODEProblem(osys2, u0map, tspan, pmap)
osol   = solve(oprob, Tsit5())
```

## How to override mass action kinetics rate laws?
While generally one wants the reaction rate law to use the law of mass action,
so the reaction
```julia
rn = @reaction_network begin
  k, X --> ∅
end
```
occurs at the (ODE) rate ``d[X]/dt = -k[X]``, it is possible to override this by
using any of the following non-filled arrows when declaring the reaction: `<=`, `⇐`,
`⟽`, `=>`, `⇒`, `⟾`, `⇔`, `⟺` (`<=>` currently not possible due to Julia langauge technical reasons). This means that the reaction
```julia
rn = @reaction_network begin
  k, X => ∅
end
```
will occur at rate ``d[X]/dt = -k`` (which might become a problem since ``[X]``
will be degraded at a constant rate even when very small or equal to 0).

Note, stoichiometric coefficients are still included, i.e. the reaction
```julia
rn = @reaction_network begin
  k, 2*X ⇒ ∅
end
```
has rate ``d[X]/dt = -2 k``.

## [How to specify user-defined functions as reaction rates?](@id user_functions)
The reaction network DSL can "see" user-defined functions that work with
ModelingToolkit. e.g., this is should work
```julia
myHill(x) = 2.0*x^3/(x^3+1.5^3)
rn = @reaction_network begin
  myHill(X), ∅ → X
end
```
In some cases, it may be necessary or desirable to register functions with
Symbolics.jl before their use in Catalyst, see the discussion
[here](https://symbolics.juliasymbolics.org/dev/manual/functions/).
