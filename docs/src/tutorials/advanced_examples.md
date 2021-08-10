# Advanced Chemical Reaction Network Examples
For additional flexibility, we can convert the generated `ReactionSystem` first
to another `ModelingToolkit.AbstractSystem`, e.g., an `ODESystem`, `SDESystem`,
`JumpSystem`, etc. These systems can then be used in problem generation. Please
also see the [ModelingToolkit](http://mtk.sciml.ai/stable/) docs, which give
many options for optimized problem generation (i.e., generating dense or sparse
Jacobians with or without threading and/or parallelization), creating LaTeX
representations for systems, etc.

Note, when generating problems from other system types, `u0` and `p` must
provide vectors of `Pair`s that map each species or parameter to their numerical
value. E.g., for the Michaelis-Menten example above we'd use
```julia
rs = @reaction_network begin
  c1, X --> 2X
  c2, X --> 0
  c3, 0 --> X
end c1 c2 c3
p     = (1.0,2.0,50.)
tspan = (0.,4.)
u0    = [5.]   
osys  = convert(ODESystem, rs)
u0map = map((x,y) -> Pair(x,y), species(rs), u0)
pmap  = map((x,y) -> Pair(x,y), params(rs), p)
oprob = ODEProblem(osys, u0map, tspan, pmap)
sol   = solve(oprob, Tsit5())
```

#### Example: Disabling rescaling of reaction rates
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

#### Example: Modifying generated ODEs by adding forcing
Conversion to other `ModelingToolkit.AbstractSystem`s allows the possibility to
modify the system with further terms that are difficult to encode as a chemical
reaction. For example, suppose we wish to add a forcing term, $10\sin(10t)$, to
the ODE for `dX/dt` above. We can do so as:
```julia
dXdteq = equations(osys)[1]           
t      = get_iv(osys)    
dXdteq = Equation(dXdteq.lhs, dXdteq.rhs + 10*sin(10*t))   
osys2  = ODESystem([dXdteq], t, states(osys), parameters(osys))
oprob  = ODEProblem(osys2, u0map, tspan, pmap)
osol   = solve(oprob, Tsit5())
```
We can add $e^{-X}$ to $dX/dt$ as a forcing term by
```julia
dXdteq = equations(osys)[1]           
@variables X(t)
dXdteq = Equation(dXdteq.lhs, dXdteq.rhs + exp(-X))   
osys2  = ODESystem([dXdteq], t, states(osys), parameters(osys))
oprob  = ODEProblem(osys2, u0map, tspan, pmap)
osol   = solve(oprob, Tsit5())
```
