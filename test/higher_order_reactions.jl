using DiffEqBiological, DiffEqJump, DiffEqBase, OrdinaryDiffEq

k = rand()
# First order
r1 = Reaction(k,[1],((1,-1),(2,1)))
jump = build_jumps_from_reaction(r1)
@test  jump.rate([1,0],nothing,0.0) == k


# Second order
r2 = Reaction(k,[1,1],((1,-2),(2,1)))
jump2 = build_jumps_from_reaction(r2)
@test  jump2.rate([1,0],nothing,0.0) == 0
@test  jump2.rate([2,0],nothing,0.0) == k*2*1

# Thrid order
r3 = Reaction(k,[1,1,1],((1,-3),(2,1)))
jump3 = build_jumps_from_reaction(r3)
@test  jump3.rate([1,0],nothing,0.0) == 0
@test  jump3.rate([2,0],nothing,0.0) == 0
@test  jump3.rate([3,0],nothing,0.0) == k*3*2*1

# Mixed Third order
r21 = Reaction(k ,[1,1,2],((1,-2),(2,-1)))
jump21 = build_jumps_from_reaction(r21)
@test  jump21.rate([1,1],nothing,0.0) == 0
@test  jump21.rate([2,1],nothing,0.0) == k*2*1*1

# Mixed Third order different order
r12 = Reaction(k ,[2,1,1],((1,-2),(2,-1)))
jump12 = build_jumps_from_reaction(r12)
  jump12.rate([1,1],nothing,0.0) == 0
@test  jump12.rate([2,1],nothing,0.0) == k*2*1*1


# Solve the dimerization problem
prob = DiscreteProblem([1,0],(0.0,250.0))
jump_prob = GillespieProblem(prob,Direct(),r2)
sol = solve(jump_prob,Discrete())
@test find( x-> x!=0 ,[u[2] for u in sol.u]) == []
