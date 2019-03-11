sir_model = @min_reaction_network rn begin
    0.1/1000, s + i --> 2i
    p1, i --> r
end p1
addodes!(sir_model)
function sir_f(du,u,p,t)
    du[1] = -(0.1/1000)*u[1]*u[2]
    du[2] = (0.1/1000)*u[1]*u[2] - 0.01 * u[2]
    du[3] = p[1] * u[2]
end

for i = 1:100
    p1 = 0.01; p = [p1];
    u0 = [990. + 200 * rand(1)[1],10. + 2* rand(1)[1],0.]
    tspan = (0.,10. + 5 * rand(1)[1])

    prob1 = ODEProblem(sir_model,u0,tspan,p)
    sol1 = solve(prob1,Tsit5())
    prob2 = ODEProblem(sir_f,u0,tspan,p)
    sol2 = solve(prob2,Tsit5())
    @test sol1[end][1] == sol2[end][1]
end

model = @min_reaction_network rn begin
    (d,1000), X ↔ 0
end d
addsdes!(model)

## For Julialang/julia#28356
i = 1.0
prob = SDEProblem(model,[1000.0+i],(0.,200.),[i])
sol = solve(prob, EM(), dt = 0.01)
@test sol[end][1] < 2000
##

for i in [1., 2., 3., 4., 5.]
    prob = SDEProblem(model,[1000.0+i],(0.,200.),[i])
    sol = solve(prob, EM(), dt = 0.01)
    @test sol[end][1] < 2000
end

#Tests time dependence
equi_model = @min_reaction_network rn begin
    100*(tanh(t-100)+2), 0 --> X
    1, X --> 0
end
addodes!(equi_model)
prob1 = ODEProblem(equi_model,[100.],(0.,100.))
sol1 = solve(prob1,Tsit5())
prob2 = ODEProblem(equi_model,[100.],(0.,200.))
sol2 = solve(prob2,Tsit5())
@test 1.5*sol1[end][1] < sol2[end][1]
