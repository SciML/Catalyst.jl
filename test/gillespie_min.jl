sir_model = @min_reaction_network rnm begin
    0.1/1000, s + i --> 2i
    0.01, i --> r
end
addodes!(sir_model)
addjumps!(sir_model)

sir_prob = DiscreteProblem([999,1,0],(0.0,250.0))
sir_jump_prob = JumpProblem(sir_prob,Direct(),sir_model)
sir_sol = solve(sir_jump_prob,FunctionMap())

prob = DiscreteProblem([999,1,0],(0.0,250.0))
jump_prob = JumpProblem(prob,Direct(),sir_model)

Random.seed!(100)
sol = solve(jump_prob,FunctionMap())

nums = Int[]
@time for i in 1:1000
  jump_prob = JumpProblem(prob,Direct(),sir_model)
  sol = solve(jump_prob,FunctionMap())
  push!(nums,sol[end][3])
end
@test mean(nums) - 740 < 20

f = function (du,u,p,t)
  du[4] = u[2]*u[3]/100000 - u[1]*u[2]/100000
end

prob = ODEProblem(f,[999.0,1.0,0.0,100.0],(0.0,250.0))
jump_prob = JumpProblem(prob,Direct(),sir_model)

sol = solve(jump_prob,Tsit5())
