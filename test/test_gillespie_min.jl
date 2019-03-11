sir_model = @min_reaction_network rn begin
    0.1/1000, s + i --> 2i
    0.01, i --> r
end
addjumps!(sir_model)
sir_prob = DiscreteProblem([999,1,0],(0.0,250.0))
sir_jump_prob = JumpProblem(sir_prob,Direct(),sir_model)

sir_sol = solve(sir_jump_prob,FunctionMap())

using Plots; plot(sir_sol)

Random.seed!(1234)
nums = Int[]
@time for i in 1:100000
  sir_sol = solve(sir_jump_prob,FunctionMap())
  push!(nums,sir_sol[end][3])
end
println("Reaction DSL: $(mean(nums))")

rate = (u,p,t) -> (0.1/1000.0)*u[1]*u[2]
affect! = function (integrator)
  integrator.u[1] -= 1
  integrator.u[2] += 1
end
jump = ConstantRateJump(rate,affect!)

rate = (u,p,t) -> 0.01u[2]
affect! = function (integrator)
  integrator.u[2] -= 1
  integrator.u[3] += 1
end
jump2 = ConstantRateJump(rate,affect!)


prob = DiscreteProblem([999.0,1.0,0.0],(0.0,250.0))
jump_prob = JumpProblem(prob,Direct(),jump,jump2)
sol = solve(jump_prob,FunctionMap())

using Plots; plot(sol)

nums = Int[]
@time for i in 1:100000
  sol = solve(jump_prob,SSAStepper())
  push!(nums,sol[end][3])
end
println("Direct Jumps: $(mean(nums))")


using Gillespie

function F(x,parms)
  (S,I,R) = x
  (beta,gamma) = parms
  infection = beta*S*I
  recovery = gamma*I
  [infection,recovery]
end

x0 = [999,1,0]
nu = [[-1 1 0];[0 -1 1]]
parms = [0.1/1000.0,0.01]
tf = 250.0
Random.seed!(1234)

nums = Int[]
@time for i in 1:100000
  result = ssa(copy(x0),F,nu,parms,tf)
  data = ssa_data(result)
  push!(nums,data[:x3][end])
end
println("Gillespie: $(mean(nums))")
