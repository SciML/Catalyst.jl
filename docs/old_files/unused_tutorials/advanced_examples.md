# Advanced Chemical Reaction Network Examples

For additional flexibility, we can convert the generated `ReactionSystem` first
to another `ModelingToolkit.AbstractSystem`, e.g., an `ODESystem`, `SDESystem`,
`JumpSystem`, etc. These systems can then be used in problem generation. Please
also see the [ModelingToolkit](http://docs.sciml.ai/ModelingToolkit/stable/) docs, which give
many options for optimized problem generation (i.e., generating dense or sparse
Jacobians with or without threading and/or parallelization), creating LaTeX
representations for systems, etc.

Note, when generating problems from other system types, `u0` and `p` must
provide vectors, tuples or dictionaries of `Pair`s that map each the symbolic
variables for each species or parameter to their numerical value. E.g., for the
Michaelis-Menten example above we'd use
```julia
rs = @reaction_network begin
  c1, X --> 2X
  c2, X --> 0
  c3, 0 --> X
end
p     = (:c1 => 1.0, :c2 => 2.0, :c3 => 50.)
pmap  = symmap_to_varmap(rs,p)   # convert Symbol map to symbolic variable map
tspan = (0.,4.)
u0    = [:X => 5.]
u0map = symmap_to_varmap(rs,u0)  # convert Symbol map to symbolic variable map
osys  = convert(ODESystem, rs)
oprob = ODEProblem(osys, u0map, tspan, pmap)
sol   = solve(oprob, Tsit5())
```
