using DifferentialEquations

infection = Reaction(0.1/1000,[1,2],[(1,-1),(2,1)])
recovery = Reaction(0.01,[2],[(2,-1),(3,1)])
sir_prob = DiscreteProblem([999,1,0],(0.0,250.0))
sir_jump_prob = GillespieProblem(sir_prob,Direct(),infection,recovery)

sir_sol = solve(sir_jump_prob,Discrete())

using Plots; plot(sir_sol)

srand(1234)
nums = Int[]
@time for i in 1:100000
  sir_sol = solve(sir_jump_prob,Discrete())
  push!(nums,sir_sol[end][3])
end
println("Reaction DSL: $(mean(nums))")


using DiffEqJump, DiffEqBase, OrdinaryDiffEq
using Base.Test

rate = (t,u) -> (0.1/1000.0)*u[1]*u[2]
affect! = function (integrator)
  integrator.u[1] -= 1
  integrator.u[2] += 1
end
jump = ConstantRateJump(rate,affect!;save_positions=(false,true))

rate = (t,u) -> 0.01u[2]
affect! = function (integrator)
  integrator.u[2] -= 1
  integrator.u[3] += 1
end
jump2 = ConstantRateJump(rate,affect!;save_positions=(false,true))


prob = DiscreteProblem([999.0,1.0,0.0],(0.0,250.0))
jump_prob = JumpProblem(prob,Direct(),jump,jump2)
sol = solve(jump_prob,Discrete(apply_map=false))

using Plots; plot(sol)

nums = Int[]
@time for i in 1:100000
  sol = solve(jump_prob,Discrete(apply_map=false))
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
srand(1234)

nums = Int[]
@time for i in 1:100000
  result = ssa(x0,F,nu,parms,tf)
  data = ssa_data(result)
  push!(nums,data[:x3][end])
end
println("Gillespie: $(mean(nums))")
