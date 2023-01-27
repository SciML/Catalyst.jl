using Catalyst, Test, OrdinaryDiffEq, LinearAlgebra

rn = @reaction_network begin
    @parameters a1 a2 k1 k2 b1
    (a1, a2), C <--> 0
    (k1, k2), A + B <--> C
    b1, 0 <-- B
end

rn2 = @reaction_network begin
    @parameters a1 a2 k1 k2 b1
    a1, C --> 0
    a2, 0 --> C
    k1, A + B --> C
    k2, C --> A + B
    b1, B --> 0
end

ps = ones(5)
u0 = [10.0, 20.0, 0.0]
tspan = (0, 100.0)
oprob1 = ODEProblem(rn, u0, tspan, ps)
oprob2 = ODEProblem(rn2, u0, tspan, ps)
sol1 = solve(oprob1, Tsit5())
sol2 = solve(oprob2, Tsit5())

tv = range(tspan[1], tspan[2], length = 100)
@test norm(sol1 - sol2, Inf) ≈ 0
