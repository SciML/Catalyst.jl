# FAQs

## How to index solution objects using symbolic variables and observables?
One can directly use symbolic variables to index into SciML solution objects.
Moreover, observables can also be evaluated in this way. For example,
consider the system
```@example faq1
using Catalyst, OrdinaryDiffEq, Plots
rn = @reaction_network ABtoC begin
  (k₊,k₋), A + B <--> C
end
nothing    # hide
```
Let's convert it to a system of ODEs, using the conservation laws of the system
to eliminate two of the species:
```@example faq1
osys = convert(ODESystem, rn; remove_conserved = true)
osys = complete(osys)
```
Notice the resulting ODE system has just one ODE, while algebraic observables
have been added for the two removed species (in terms of the conservation law
constants, `Γ[1]` and `Γ[2]`)
```@example faq1
observed(osys)
```
Let's solve the system and see how to index the solution using our symbolic
variables
```@example faq1
u0 = [rn.A => 1.0, rn.B => 2.0, rn.C => 0.0]
ps = [rn.k₊ => 1.0, rn.k₋ => 1.0]
oprob = ODEProblem(osys, [], (0.0, 10.0), [])
sol = solve(oprob, Tsit5())
```
Suppose we want to plot just species `C`, without having to know its integer
index in the unknown vector. We can do this using the symbolic variable `C`, which
we can get at in several ways
```julia
sol[osys.C]
```
or
```julia
@unpack C = osys
sol[C]
```
To evaluate `C` at specific times and plot it we can just do
```julia
t = range(0.0, 10.0, length=101)
plot(t, sol(t, idxs = C), label = "C(t)", xlabel = "t")
```
If we want to get multiple variables we can just do
```julia
@unpack A, B = osys
sol(t, idxs = [A, B])
```
Plotting multiple variables using the SciML plot recipe can be achieved
like
```julia
plot(sol; idxs = [A, B])
```

## How to disable rescaling of reaction rates in rate laws?
As explained in the [Reaction rate laws used in simulations](@ref) section, for
a reaction such as `k, 2X --> 0`, the generated rate law will rescale the rate
constant, giving `k*X^2/2` instead of `k*X^2` for ODEs and `k*X*(X-1)/2` instead
of `k*X*(X-1)` for jumps. This can be disabled when directly `convert`ing a
[`ReactionSystem`](@ref). If `rn` is a generated [`ReactionSystem`](@ref), we can
do
```@example faq1
osys = convert(ODESystem, rn; combinatoric_ratelaws=false)
```
Disabling these rescalings should work for all conversions of `ReactionSystem`s
to other `ModelingToolkit.AbstractSystem`s.

## How to use non-integer stoichiometric coefficients?
```@example faq2
using Catalyst
rn = @reaction_network begin
  k, 2.5*A --> 3*B
end
```
or directly via
```@example faq2
t = default_t()
@parameters k b
@species A(t) B(t) C(t) D(t)
rx1 = Reaction(k,[B,C],[B,D], [2.5,1],[3.5, 2.5])
rx2 = Reaction(2*k, [B], [D], [1], [2.5])
rx3 = Reaction(2*k, [B], [D], [2.5], [2])
@named mixedsys = ReactionSystem([rx1, rx2, rx3], t, [A, B, C, D], [k, b])
osys = convert(ODESystem, mixedsys; combinatoric_ratelaws = false)
osys = complete(osys)
```
Note, when using `convert(ODESystem, mixedsys; combinatoric_ratelaws=false)` the
`combinatoric_ratelaws=false` parameter must be passed. This is also true when
calling `ODEProblem(mixedsys,...; combinatoric_ratelaws=false)`. As described
above, this disables Catalyst's standard rescaling of reaction rates when
generating reaction rate laws, see also the [Reaction rate laws used in
simulations](@ref) section. Leaving this keyword out for systems with floating
point stoichiometry will give an error message.

For a more extensive documentation of using non-integer stoichiometric
coefficients, please see the [Symbolic Stochiometries](@ref
parametric_stoichiometry) section.

## How to set default values for initial conditions and parameters?
How to set defaults when using the `@reaction_network` macro is described in
more detail [here](@ref dsl_description_defaults). There are several ways to do
this. Using the DSL, one can use the `@species` and `@parameters` options:
```@example faq3
using Catalyst
sir = @reaction_network sir begin
    @species S(t)=999.0 I(t)=1.0 R(t)=0.0
    @parameters β=1e-4 ν=0.01
    β, S + I --> 2I
    ν, I --> R
end
show(stdout, MIME"text/plain"(), sir) # hide
```

When directly constructing a `ReactionSystem`, we can set the symbolic values to
have the desired default values, and this will automatically be propagated
through to the equation solvers:
```@example faq3
using Catalyst, Plots, OrdinaryDiffEq
t = default_t()
@parameters β=1e-4 ν=.01
@species S(t)=999.0 I(t)=1.0 R(t)=0.0
rx1 = Reaction(β, [S, I], [I], [1,1], [2])
rx2 = Reaction(ν, [I], [R])
@named sir = ReactionSystem([rx1, rx2], t)
sir = complete(sir)
oprob = ODEProblem(sir, [], (0.0, 250.0))
sol = solve(oprob, Tsit5())
plot(sol)
```

One can also build a mapping from symbolic parameter/species to value/initial
condition and pass these to the `ReactionSystem` via the `defaults` keyword
argument:
```@example faq3
@parameters β ν
@species S(t) I(t) R(t)
rx1 = Reaction(β, [S,I], [I], [1,1], [2])
rx2 = Reaction(ν, [I], [R])
defs = [β => 1e-4, ν => .01, S => 999.0, I => 1.0, R => 0.0]
@named sir = ReactionSystem([rx1, rx2], t; defaults = defs)
nothing # hide
```

Finally, default values can also be added after creating the system via the
`setdefaults!` command and passing a `Symbol` based mapping, like
```@example faq3
sir = @reaction_network sir begin
    β, S + I --> 2I
    ν, I --> R
end
setdefaults!(sir, [:β => 1e-4, :ν => .01, :S => 999.0, :I => 1.0, :R => 0.0])
nothing # hide
```

## How to specify initial conditions and parameters values for `ODEProblem` and other problem types?
To explicitly pass initial conditions and parameters we can use mappings from
Julia `Symbol`s corresponding to each variable/parameter to their values, or
from ModelingToolkit symbolic variables/parameters to their values. Using
`Symbol`s we have
```@example faq4
using Catalyst, OrdinaryDiffEq
rn = @reaction_network begin
    α, S + I --> 2I
    β, I --> R
end
u0 = [:S => 999.0, :I => 1.0, :R => 0.0]
p  = (:α => 1e-4, :β => .01)
op1  = ODEProblem(rn, u0, (0.0, 250.0), p)
nothing  # hide
```
while using ModelingToolkit symbolic variables we have
```@example faq4
t = default_t()
@parameters α β
@species S(t) I(t) R(t)
u0 = [S => 999.0, I => 1.0, R => 0.0]
p  = (α => 1e-4, β => .01)
op2  = ODEProblem(rn, u0, (0.0, 250.0), p)
nothing  # hide
```
*Note,* while symbolic mappings as in the last example will work with *any*
`ModelingToolkit.AbstractSystem`, for example if one `convert`s `rn` to an
`ODESystem`, `Symbol`-based mappings only work when passing a `ReactionSystem`
directly into a problem type. That is, the following does not work
```@julia
osys = convert(ODESystem, rn)

# this fails
u0 = [:S => 999.0, :I => 1.0, :R => 0.0]
p  = (:α => 1e-4, :β => .01)
op  = ODEProblem(osys, u0, (0.0, 250.0), p)
```
In this case one must either use a symbolic mapping as was used to make `op2` in
the second example, or one can use the `symmap_to_varmap` function to convert the
`Symbol` mapping to a symbolic mapping. I.e. this works
```@example faq4
osys = convert(ODESystem, rn)
osys = complete(osys)

# this works
u0 = symmap_to_varmap(rn, [:S => 999.0, :I => 1.0, :R => 0.0])
p  = symmap_to_varmap(rn, (:α => 1e-4, :β => .01))
op  = ODEProblem(osys, u0, (0.0, 250.0), p)
nothing # hide
```

## How to include non-reaction terms in equations for a chemical species?
One method to add non-reaction terms into an ODE or algebraic equation for a
chemical species is to add a new (non-species) unknown variable that represents
those terms, let it be the rate of zero order reaction, and add a constraint
equation. I.e., to add a force of `(1 + sin(t))` to ``dA/dt`` in a system with
the reaction `k, A --> 0`, we can do
```@example faq5
using Catalyst
t = default_t()
@variables f(t)
rx1 = @reaction k, A --> 0
rx2 = @reaction $f, 0 --> A
eq = f ~ (1 + sin(t))
@named rs = ReactionSystem([rx1, rx2, eq], t)
rs = complete(rs)
osys = convert(ODESystem, rs)
```
In the final ODE model, `f` can be eliminated by using
`ModelingToolkit.structural_simplify`
```@example faq5
osyss = structural_simplify(osys)
full_equations(osyss)
```

## How to modify generated ODEs?
Conversion to other `ModelingToolkit.AbstractSystem`s allows the possibility to
modify the system with further terms that are difficult to encode as a chemical
reaction or a constraint equation. For example, an alternative method to the
previous question for adding a forcing function, $1 + \sin(t)$, to the ODE for
`dA/dt` is
```@example faq6
using Catalyst
rn = @reaction_network begin
    k, A --> 0
end
osys = convert(ODESystem, rn)
dAdteq = equations(osys)[1]
t      = ModelingToolkit.get_iv(osys)
dAdteq = Equation(dAdteq.lhs, dAdteq.rhs + 1 + sin(t))

# create a new ODESystem with the modified equation
@named osys2  = ODESystem([dAdteq], t)
```

## How to override mass action kinetics rate laws?
While generally one wants the reaction rate law to use the law of mass action,
so the reaction
```@example faq7
using Catalyst
rn = @reaction_network begin
  k, X --> ∅
end
convert(ODESystem, rn)
```
occurs at the (ODE) rate ``d[X]/dt = -k[X]``, it is possible to override this by
using any of the following non-filled arrows when declaring the reaction: `<=`,
`⇐`, `⟽`, `=>`, `⇒`, `⟾`, `⇔`, `⟺`. This means that the reaction
```@example faq7
rn = @reaction_network begin
  k, X => ∅
end
convert(ODESystem, rn)
```
will occur at rate ``d[X]/dt = -k`` (which might become a problem since ``[X]``
will be degraded at a constant rate even when very small or equal to 0).

Note, stoichiometric coefficients are still included, i.e. the reaction
```@example faq7
rn = @reaction_network begin
  k, 2*X ⇒ ∅
end
convert(ODESystem, rn)
```
has rate ``d[X]/dt = -2 k``.

## [How to specify user-defined functions as reaction rates?](@id user_functions)
The reaction network DSL can "see" user-defined functions that work with
ModelingToolkit. e.g., this is should work
```@example faq8
using Catalyst
myHill(x) = 2*x^3/(x^3+1.5^3)
rn = @reaction_network begin
  myHill(X), ∅ --> X
end
```
In some cases, it may be necessary or desirable to register functions with
Symbolics.jl before their use in Catalyst, see the discussion
[here](https://symbolics.juliasymbolics.org/stable/manual/functions/).
