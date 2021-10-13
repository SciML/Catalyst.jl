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
pmap  = map((x,y) -> Pair(x,y), parameters(rs), p)
oprob = ODEProblem(osys, u0map, tspan, pmap)
sol   = solve(oprob, Tsit5())
```

