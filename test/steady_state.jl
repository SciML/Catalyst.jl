model = @reaction_network SSTest begin
    (d_x, d_y), (x, y) --> 0
    mm(y, 2, 1), 0 --> x
    1., 0 --> y
end d_x d_y

u0 = [0.1, 0.1]
p = [1., 1.]

prob = SteadyStateProblem{true}(model, u0, p)
sol = solve(prob, SSRootfind())
@test sol.u ≈ [1, 1]

u0 = [0.1, 0.1]
prob = SteadyStateProblem(model, u0, p)
sol = solve(prob, SSRootfind())
@test sol.u ≈ [1, 1]
