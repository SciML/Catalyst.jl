# FAQs

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
end k
```
occurs at the (ODE) rate ``d[X]/dt = -k[X]``, it is possible to override this by
using any of the following non-filled arrows when declaring the reaction: `⇐`,
`⟽`, `⇒`, `⟾`, `⇔`, `⟺`. This means that the reaction
```julia
rn = @reaction_network begin
  k, X ⇒ ∅
end k
```
will occur at rate ``d[X]/dt = -k`` (which might become a problem since ``[X]``
will be degraded at a constant rate even when very small or equal to 0). 

Note, stoichiometric coefficients are still included, i.e. the reaction
```julia
rn = @reaction_network begin
  k, 2*X ⇒ ∅
end k
```
has rate ``d[X]/dt = -2 k``.

## [How to specify user defined functions as reaction rates?](@id user_functions)
The reaction network DSL can "see" user defined functions that work with
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


