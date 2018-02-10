using DiffEqBiological, DiffEqJump, DiffEqBase, OrdinaryDiffEq, StochasticDiffEq
using Base.Test
using tmpMod

sir_model = @reaction_network rn begin
    0.1/1000, s + i --> 2i
    0.01, i --> r
end

prob = DiscreteProblem([999.0,3.0,0],(0.0,250.0))
jump_prob = JumpProblem(prob,Direct(),sir_model)

f = function (du,u,p,t)
  du[4] = u[2]*u[3]/100000 - u[1]*u[2]/100000
end

prob = ODEProblem(f,[999.0,3.0,0.0,100.0],(0.0,250.0))
jump_prob = JumpProblem(prob,Direct(),sir_model)

srand(230)
println("Add ODE to Gillespie Problem")
sol = solve(jump_prob,Tsit5())
@test 650 <= sol[end][3] <= 900

sir_model = @reaction_network rn begin
    0.1/1000, s + i --> 2i
    0.01, i --> r
    0.01d, d + i ⟾ d + r
end

prob = ODEProblem(f,[999.0,3.0,0.0,1.0],(0.0,250.0))
jump_prob = JumpProblem(prob,Direct(),sir_model)
sol = solve(jump_prob,Tsit5())

g = function (du,u,p,t)
  du[4] = 0.05u[4]
end

srand(331)
println("Turn Gillespie Problem into SDE")
prob = SDEProblem(f,g,[999.0,3.0,0.0,1.0],(0.0,250.0))
jump_prob = JumpProblem(prob,Direct(),sir_model)
sol = solve(jump_prob,SRIW1())
@test 650 <= sol[end][3] <= 1500

println("Now mix constant rate reactions")
jump_prob = JumpProblem(prob,Direct(),sir_model)
integrator = init(jump_prob,SRIW1())
solve!(integrator)
@test 650 <= sol[end][3] <= 1500
