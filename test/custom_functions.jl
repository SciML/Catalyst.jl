using DiffEqBase, Catalyst, Random, Test
using ModelingToolkit: get_states, get_ps

using StableRNGs
rng = StableRNG(12345)

### Tests various cutom made functions ###
new_hill(x, v, k, n) = v*x^n/(k^n+x^n)
new_poly(x,p1,p2) = .5*p1*x^2
new_exp(x,p) = exp(-p*x)

custom_function_network_1 = @reaction_network begin
    hill(X1,v1,K1,2), X1 + Y1--> Z1
    mm(X2,v2,K2), X2 + Y2 --> Z2
    p1*X3^2+p2, X3 + Y3 --> Z3
    exp(-p3*Y4), X4 + Y4 --> Z4
    hillR(X5,v3,K3,2), X5 + Y5 --> Z5
    mmR(X6,v4,K4), X6 + Y6 --> Z6
    hillAR(X7,Y7,v5,K5,2), X7 + Y7 --> Z7
end v1 K1 v2 K2 p1 p2 p3 v3 K3 v4 K4 v5 K5

custom_function_network_2 = @reaction_network begin
    new_hill(X1,v1,K1,2), X1 + Y1 --> Z1
    v2*X2/(X2+K2), X2 + Y2 --> Z2
    2*new_poly(X3,p1,p2)+p2, X3 + Y3 --> Z3
    new_exp(Y4,p3), X4 + Y4--> Z4
    v3*(K3^2)/(K3^2+X5^2), X5 + Y5 --> Z5
    v4*K4/(X6+K4), X6 + Y6 --> Z6
    v5*(X7^2)/(K5^2+X7^2+Y7^2), X7 + Y7 --> Z7
end v1 K1 v2 K2 p1 p2 p3 v3 K3 v4 K4 v5 K5

f1 = ODEFunction(convert(ODESystem,custom_function_network_1),jac=true)
f2 = ODEFunction(convert(ODESystem,custom_function_network_2),jac=true)
g1 = SDEFunction(convert(SDESystem,custom_function_network_1))
g2 = SDEFunction(convert(SDESystem,custom_function_network_2))
for factor in [1e-2, 1e-1, 1e0, 1e1, 1e2]
    u0 = factor*rand(rng,length(get_states(custom_function_network_1)))
    p = factor*rand(rng,length(get_ps(custom_function_network_2)))
    t = rand(rng)
    @test all(abs.(f1(u0,p,t) .- f2(u0,p,t)) .< 10e-10)
    @test all(abs.(f1.jac(u0,p,t) .- f2.jac(u0,p,t)) .< 10e-10)
    @test all(abs.(g1(u0,p,t) .- g2(u0,p,t)) .< 10e-10)
end