using Catalyst, LinearAlgebra, DiffEqJump, OrdinaryDiffEq, Test

rs = @reaction_network begin
  (c1,c2), A + B <--> C + D
end c1 c2 

u0 = [5, 5, 0, 0]
u0 = [10, 15, 3, 1]
specmap = Num.(states(rs)) .=> u0
state_space = Catalyst.traverse_reactionsystem(rs, u0)
cmesys = Catalyst.ODESystem_cme(rs, u0)
prob = ODEProblem(cmesys, [], (0, 10.), [1., 2.])

sol = solve(prob, Tsit5())
arr = Array(sol)
@test all(sum(arr; dims=1) .≈ 1.) # sanity check that states sum to probability 1

# birth death (infinte)
rs = @reaction_network begin
  (c1, c2), X <--> 2X
end c1 c2
u0 = [2]

@test_throws ErrorException traverse_reactionsystem(rs, u0)
@test length(traverse_reactionsystem(rs, u0; truncation=[10])) == 11
