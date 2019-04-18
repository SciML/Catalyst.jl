function tmp_sol2vec(sol,j)
    vect = Vector{Float64}(undef,length(sol.u))
    for i = 1:length(sol.u)
        vect[i] = sol.u[i][j]
    end
    return vect
end

ho_model1 = @reaction_network rn begin
    (1.0,1.0), X + 2X ↔ Y + X
end
ho_model2 = @reaction_network rn begin
    (X^3/6,Y * X), X + 2X ⟺ Y + X
end
ho_model3 = @reaction_network rn begin
    (X*(X-1)*(X-2)/6,Y * X), X + 2X ⟺ Y+ X
end
u0 = [1000.,1000.]
tspan = (0.,100.)

prob_det1 = ODEProblem(ho_model1, u0, tspan)
prob_det2 = ODEProblem(ho_model2, u0, tspan)
sol_det1 = solve(prob_det1,Tsit5(), atol=1e-8)
sol_det2 = solve(prob_det2,Tsit5(), atol=1e-8)
#@test tmp_sol2vec(sol_det1,1) == tmp_sol2vec(sol_det2,1)
#@test tmp_sol2vec(sol_det1,2) == tmp_sol2vec(sol_det2,2)
@test isapprox(tmp_sol2vec(sol_det1,1), tmp_sol2vec(sol_det2,1), atol=1e-7)
@test isapprox(tmp_sol2vec(sol_det1,2), tmp_sol2vec(sol_det2,2), atol=1e-7)

disc_prob =  DiscreteProblem([1000,1000],(0.,1.))
jump_prob1 = JumpProblem(disc_prob,Direct(),ho_model1)
jump_prob2 = JumpProblem(disc_prob,Direct(),ho_model2)
sol_jump1 = solve(jump_prob1,FunctionMap(),maxiters = 1e6)
sol_jump2 = solve(jump_prob2,FunctionMap(),maxiters = 1e6)

@test 0.95 < std(tmp_sol2vec(sol_jump1,1)) / std(tmp_sol2vec(sol_jump2,1)) < 1.05
@test 0.95 < std(tmp_sol2vec(sol_jump1,2)) / std(tmp_sol2vec(sol_jump2,2)) < 1.05
