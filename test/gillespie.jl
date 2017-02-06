using BiologicalModels, JumpDiffEq, DiffEqBase, OrdinaryDiffEq

r1 = Reaction(1e-4,(1,2),((1,-1),(2,1)))
r2 = Reaction(0.01,[2],[(2,-1),(3,1)])

prob = DiscreteProblem([999,1,0],(0.0,250.0))
jump_prob = GillespieProblem(prob,Direct(),r1,r2)

srand(100)
sol = solve(jump_prob,Discrete())
#using Plots; plotly(); plot(sol)

nums = Int[]
@time for i in 1:1000
  jump_prob = GillespieProblem(prob,Direct(),r1,r2)
  sol = solve(jump_prob,Discrete())
  push!(nums,sol[end][3])
end
mean(nums)

f = function (t,u,du)
  du[4] = u[2]*u[3]/100000 - u[1]*u[2]/100000
end

prob = ODEProblem(f,[999.0,1.0,0.0,100.0],(0.0,250.0))
jump_prob = GillespieProblem(prob,Direct(),r1,r2)

sol = solve(jump_prob,Tsit5())
#plot(sol)
