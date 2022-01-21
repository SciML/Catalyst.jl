# Basic Chemical Reaction Network Examples

#### Example: Birth-Death Process

```julia
rs = @reaction_network begin
  c1, X --> 2X
  c2, X --> 0
  c3, 0 --> X
end c1 c2 c3
p = (:c1 => 1.0, :c2 => 2.0, :c3 => 50.) 
tspan = (0.,4.)
u0 = [:X => 5.]

# solve ODEs
oprob = ODEProblem(rs, u0, tspan, p)
osol  = solve(oprob, Tsit5())

# solve for Steady-States
ssprob = SteadyStateProblem(rs, u0, p)
sssol  = solve(ssprob, SSRootfind())

# solve SDEs
sprob = SDEProblem(rs, u0, tspan, p)
ssol  = solve(sprob, EM(), dt=.01)

# solve JumpProblem
u0 = [:X => 5]
dprob = DiscreteProblem(rs, u0, tspan, p)
jprob = JumpProblem(rs, dprob, Direct())
jsol = solve(jprob, SSAStepper())
```

#### Example: Michaelis-Menten Enzyme Kinetics

```julia
rs = @reaction_network begin
  c1, S + E --> SE
  c2, SE --> S + E
  c3, SE --> P + E
end c1 c2 c3
p = (:c1 => 0.00166, :c2 => 0.0001, :c3 => 0.1) 
tspan = (0., 100.)
u0 = [:S => 301., :E => 100., :SE => 0., :P => 0.]  

# solve ODEs
oprob = ODEProblem(rs, u0, tspan, p)
osol  = solve(oprob, Tsit5())

# solve JumpProblem
u0 = [:S => 301, :E => 100, :SE => 0, :P => 0] 
dprob = DiscreteProblem(rs, u0, tspan, p)
jprob = JumpProblem(rs, dprob, Direct())
jsol = solve(jprob, SSAStepper())
```