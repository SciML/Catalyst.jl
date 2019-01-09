using DiffEqBiological, Test, DiffEqBase

network1 = @reaction_network rn begin
    2.0, X + Y --> XY
    1.0, XY ← X + Y
    (1,X), Z ↔ Z1 + Z2
end

get_odefun!(network1)
get_sdefun!(network1)
gen_jumpfun!(network1)
@test length(network1.properties[:f_func]) == 6
@test typeof(network1.properties[:f]) <: Function
@test length(network1.properties[:g_func]) == 24
@test typeof(network1.properties[:g]) <: Function
@test size(network1.properties[:p_matrix]) == (6,4)
@test length(network1.params) == 0
@test length(network1.syms) == 6
@test typeof(network1.properties[:jumps][1]) <:ConstantRateJump
@test length(network1.properties[:jump_rate_expr]) == 4
@test length(network1.properties[:jump_affect_expr]) == 4

network2 = @reaction_network rn begin
    1.0, ∅ > X1
    X, 0 → X1
    (1.0,X1,X2), (X1,X2,X3) --> (X1,X2,X3)
    (1.0,X2,X1), (X2,X1,X3) --> X1
    ((1.0,2,3.), (P1,X2,X1)), X2 ↔ (X2+X4,X1,X2)
    hill(P2,2.,1,1), asdasdrqr < X1
    mm(X1,2.,1), X1 ⟾ X5
end
get_odefun!(network2)
get_sdefun!(network2)
gen_jumpfun!(network2)
@test length(network2.properties[:f_func]) == 6
@test typeof(network2.properties[:f]) <: Function
@test length(network2.properties[:g_func]) == 96
@test typeof(network2.properties[:g]) <: Function
@test size(network2.properties[:p_matrix]) == (6,16)
@test length(network2.:params) == 0
@test length(network2.syms) == 6
@test typeof(network2.properties[:jumps][1]) <:ConstantRateJump
@test length(network2.properties[:jump_rate_expr]) == 16
@test length(network2.properties[:jump_affect_expr]) == 16

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

get_odefun!(network4)
get_sdefun!(network4)
gen_jumpfun!(network4)
get_odefun!(network5)
get_sdefun!(network5)
gen_jumpfun!(network5)
@test network4.properties[:f_func] == network5.properties[:f_func]
@test network4.properties[:g_func] == network5.properties[:g_func]
@test network4.properties[:p_matrix] == network5.properties[:p_matrix]
@test network4.params == network5.params
@test network4.properties[:jump_rate_expr] == network5.properties[:jump_rate_expr]
@test network4.properties[:jump_affect_expr] == network5.properties[:jump_affect_expr]

for i = 1:100
    u = 5*rand(4)
    du4 = 3*rand(4); du5 = du4;
    du4g = 2.5*rand(4,3); du5g = du4g;
    t = 9*rand(1)[1]
    p12 = 2*rand(1)[1]; p =[p12]

    @test network4.properties[:f](du4,u,p,t) == network5.properties[:f](du5,u,p,t)
    @test network4.properties[:g](du4g,u,p,t) == network5.properties[:g](du5g,u,p,t)
    @test network4.properties[:jumps][1].rate(u,p,t) == network5.properties[:jumps][1].rate(u,p,t)
    @test network4.properties[:jumps][2].rate(u,p,t) == network5.properties[:jumps][2].rate(u,p,t)
    @test network4.properties[:jumps][3].rate(u,p,t) == network5.properties[:jumps][3].rate(u,p,t)
end

network6 = @reaction_network begin
    p2, X + Y --> XY
    p1, XY ← X + Y
    (1,X), Z ↔ Z1 + Z2
end p1 p2

get_odefun!(network6)
get_sdefun!(network6)
@test length(network6.properties[:f_func]) == 6
@test length(network6.properties[:g_func]) == 24
@test size(network6.properties[:p_matrix]) == (6,4)
@test length(network6.params) == 2
@test length(network6.syms) == 6
@test typeof(network6) <: DiffEqBase.AbstractReactionNetwork
@test typeof(network6) == reaction_network
