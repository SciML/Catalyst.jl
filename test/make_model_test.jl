network1 = @reaction_network rn begin
    2.0, X + Y --> XY
    1.0, XY ← X + Y
    (1,X), Z ↔ Z1 + Z2
end

@test length(network1.f_func) == 6
@test typeof(network1.f) <: Function
@test length(network1.g_func) == 24
@test typeof(network1.g) <: Function
@test size(network1.p_matrix) == (6,4)
@test length(network1.params) == 0
@test length(network1.syms) == 6
@test typeof(network1.jumps[1]) <:ConstantRateJump
@test length(network1.jump_rate_expr) == 4
@test length(network1.jump_affect_expr) == 4

network2 = @reaction_network rn begin
    1.0, ∅ > X1
    X, 0 → X1
    (1.0,X1,X2), (X1,X2,X3) --> (X1,X2,X3)
    (1.0,X2,X1), (X2,X1,X3) --> X1
    ((1.0,2,3.), (P1,X2,X1)), X2 ↔ (X2+X4,X1,X2)
    hill(P2,2.,1,1), asdasdrqr < X1
    mm(X1,2.,1), X1 ⟾ X5
end

@test length(network2.f_func) == 6
@test typeof(network2.f) <: Function
@test length(network2.g_func) == 96
@test typeof(network2.g) <: Function
@test size(network2.p_matrix) == (6,16)
@test length(network2.params) == 0
@test length(network2.syms) == 6
@test typeof(network2.jumps[1]) <:ConstantRateJump
@test length(network2.jump_rate_expr) == 16
@test length(network2.jump_affect_expr) == 16

network3 = @reaction_network rn2 begin
    p3, X + Y --> XY
    p2, XY ← X + Y
    (p1,X), Z ↔ Z1 + Z2
end p1 p2 p3

@test length(network3.params) == 3
@test typeof(network3) == rn2
@test typeof(network1) == typeof(network2)

network4 = @reaction_network rn begin
    1, X + Y > XY
    p1, X + Y < XY
    2, Z --> 0
end p1
network5 = @reaction_network rn begin
    (1,p1), X + Y ↔ XY
    2Z, ∅ ⟽ Z
end p1

@test network4.f_func == network5.f_func
@test network4.g_func == network5.g_func
@test network4.p_matrix == network5.p_matrix
@test network4.params == network5.params
@test network4.jump_rate_expr == network5.jump_rate_expr
@test network4.jump_affect_expr == network5.jump_affect_expr

for i = 1:100
    u = 5*rand(4)
    du4 = 3*rand(4); du5 = du4;
    du4g = 2.5*rand(4,3); du5g = du4g;
    t = 9*rand(1)[1]
    p12 = 2*rand(1)[1]; p =[p12]

    @test network4.f(du4,u,p,t) == network5.f(du5,u,p,t)
    @test network4.g(du4g,u,p,t) == network5.g(du5g,u,p,t)
    @test network4.jumps[1].rate(u,p,t) == network5.jumps[1].rate(u,p,t)
    @test network4.jumps[2].rate(u,p,t) == network5.jumps[2].rate(u,p,t)
    @test network4.jumps[3].rate(u,p,t) == network5.jumps[3].rate(u,p,t)
end

network6 = @reaction_network begin
    p2, X + Y --> XY
    p1, XY ← X + Y
    (1,X), Z ↔ Z1 + Z2
end p1 p2

@test length(network6.f_func) == 6
@test length(network6.g_func) == 24
@test size(network6.p_matrix) == (6,4)
@test length(network6.params) == 2
@test length(network6.syms) == 6
@test typeof(network6) <: DiffEqBase.AbstractReactionNetwork
@test typeof(network6) == reaction_network

# test that can we construct a network that involves symengine constants as species names
symengnet = @reaction_network begin
    k1, I + E --> 2*A
    k2, 2I --> 2E
    k4*E*A, E + A ⇒ 0
end k1 k2 k4
realjac = DiffEqBiological.ExprValues[  
    :(-(k1 * E) - 2 * (k2 * I))  :(-(k1 * I))  0
    :(-(k1 * E) + 2 * (k2 * I))  :(-(k1 * I) - k4 * A)      :(-k4 * E)  
    :(2 * (k1 * E))              :(2 * (k1 * I) - k4 * A)   :(-k4 * E)
]
@assert realjac == jacobianexprs(symengnet)


# test that we can construct a network that involves a singular Jacobian
network7 = @reaction_network begin
    k1, I + E --> 2*B
    k2, 2I --> 2E
    k4*E*A, E + A ⇒ 0
end k1 k2 k4
network8 = @min_reaction_network begin
    k1, I + E --> 2*B
    k2, 2I --> 2E
    k4*E*A, E + A ⇒ 0
end k1 k2 k4
addodes!(network8, sparse_jac=true)