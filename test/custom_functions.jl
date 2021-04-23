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

### Tests that the various notations gives identical results ###

# Michaelis-Menten function.
mm_network = @reaction_network begin
    (1.,1.), 0 ↔ X
    mm(X,v,K), 0 --> X1
    Mm(X,v,K), 0 --> X2
    MM(X,v,K), 0 --> X3
end v K
f_mm = ODEFunction(convert(ODESystem,mm_network),jac=true)

u0 = 10*rand(rng,length(get_states(mm_network)))
p = 10*rand(rng,length(get_ps(mm_network)))
t = 10*rand(rng)

f_mm_output = f_mm(u0,p,t)[2:end]
f_mm_jac_output = f_mm.jac(u0,p,t)[2:end,1]
@test (maximum(f_mm_output) - minimum(f_mm_output)) .< 100*eps()
@test (maximum(f_mm_jac_output) - minimum(f_mm_jac_output)) .< 100*eps()

# Repressing Michaelis-Menten function.
mmR_network = @reaction_network begin
    (1.,1.), 0 ↔ X
    mmR(X,v,K), 0 --> X1
    MmR(X,v,K), 0 --> X2
    MMR(X,v,K), 0 --> X3
end v K
f_mmR = ODEFunction(convert(ODESystem,mmR_network),jac=true)

u0 = 10*rand(rng,length(get_states(mmR_network)))
p = 10*rand(rng,length(get_ps(mmR_network)))
t = 10*rand(rng)

f_mmR_output = f_mmR(u0,p,t)[2:end]
f_mmR_jac_output = f_mmR.jac(u0,p,t)[2:end,1]
@test (maximum(f_mmR_output) - minimum(f_mmR_output)) .< 100*eps()
@test (maximum(f_mmR_jac_output) - minimum(f_mmR_jac_output)) .< 100*eps()

# Hill function.
hill_network = @reaction_network begin
    (1.,1.), 0 ↔ X
    hill(X,v,K,2), 0 --> X1
    Hill(X,v,K,2), 0 --> X2
end v K
f_hill = ODEFunction(convert(ODESystem,hill_network),jac=true)

u0 = 10*rand(rng,length(get_states(hill_network)))
p = 10*rand(rng,length(get_ps(hill_network)))
t = 10*rand(rng)

f_hill_output = f_hill(u0,p,t)[2:end]
f_hill_jac_output = f_hill.jac(u0,p,t)[2:end,1]
@test (maximum(f_hill_output) - minimum(f_hill_output)) .< 100*eps()
@test (maximum(f_hill_jac_output) - minimum(f_hill_jac_output)) .< 100*eps()

# Repressing Hill function.
hillR_network = @reaction_network begin
    (1.,1.), 0 ↔ X
    hillR(X,v,K,2), 0 --> X1
    HillR(X,v,K,2), 0 --> X2
end v K
f_hillR = ODEFunction(convert(ODESystem,hillR_network),jac=true)

u0 = 10*rand(rng,length(get_states(hillR_network)))
p = 10*rand(rng,length(get_ps(hillR_network)))
t = 10*rand(rng)

f_hillR_output = f_hillR(u0,p,t)[2:end]
f_hillR_jac_output = f_hillR.jac(u0,p,t)[2:end,1]
@test (maximum(f_hillR_output) - minimum(f_hillR_output)) .< 100*eps()
@test (maximum(f_hillR_jac_output) - minimum(f_hillR_jac_output)) .< 100*eps()

# Activation/repressing Hill function.
hillAR_network = @reaction_network begin
    (1.,1.), 0 ↔ (X,Y)
    hillC(X,Y,v,K,2), 0 --> X1
    HillC(X,Y,v,K,2), 0 --> X2
    hillAR(X,Y,v,K,2), 0 --> X3
    HillAR(X,Y,v,K,2), 0 --> X4
end v K
f_hillAR = ODEFunction(convert(ODESystem,hillAR_network),jac=true)

u0 = 10*rand(rng,length(get_states(hillAR_network)))
p = 10*rand(rng,length(get_ps(hillAR_network)))
t = 10*rand(rng)

f_hillAR_output = f_hillAR(u0,p,t)[3:end]
f_hillAR_jac_output = f_hillAR.jac(u0,p,t)[3:end,1]
@test (maximum(f_hillAR_output) - minimum(f_hillAR_output)) .< 100*eps()
@test (maximum(f_hillAR_jac_output) - minimum(f_hillAR_jac_output)) .< 100*eps()