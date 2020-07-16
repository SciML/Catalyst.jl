### Fetch required packages ###
using DiffEqBase, Catalyst, DiffEqJump, Random, Statistics, Test


### Declares a test network. ###
higher_order_network_1 = @reaction_network begin
    p,                                                          ∅ ⟼ X1
    r1,                                                         2X1 ⟼ 3X2
    mm(X1,r2,K),                                                3X2 ⟼ X3 + 2X4
    r3,                                                         X3 + 2X4 ⟼ 3X5 + 3X6
    r4*X2,                                                      3X5 + 3X6 ⟼ 3X5 + 2X7 + 4X8
    r5,                                                         3X5 + 2X7 + 4X8 ⟼ 10X9
    r6,                                                         10X9 ⟼ X10
    d,                                                          2X10 ⟼ ∅
end p r1 r2 K r3 r4 r5 r6 d


### Tests that deterministic and stochastic differential functions are identical. ###

higher_order_network_2 = @reaction_network begin
    p,                                                          ∅ ⟾ X1
    r1*X1^2/factorial(2),                                       2X1 ⟾ 3X2
    mm(X1,r2,K)*X2^3/factorial(3),                              3X2 ⟾ X3 + 2X4
    r3*X3*X4^2/factorial(2),                                    X3 + 2X4 ⟾ 3X5 + 3X6
    r4*X2*X5^3*X6^3/(factorial(3)*factorial(3)),                3X5 + 3X6 ⟾ 3X5 + 2X7 + 4X8
    r5*X5^3*X7^2*X8^4/(factorial(3)*factorial(2)*factorial(4)), 3X5 + 2X7 + 4X8 ⟾ 10X9
    r6*X9^10/factorial(10),                                     10X9 ⟾ X10
    d*X10^2/factorial(2),                                       2X10 ⟾ ∅
end p r1 r2 K r3 r4 r5 r6 d

f1 = ODEFunction(convert(ODESystem,higher_order_network_1),jac=true)
f2 = ODEFunction(convert(ODESystem,higher_order_network_2),jac=true)
g1 = SDEFunction(convert(SDESystem,higher_order_network_1))
g2 = SDEFunction(convert(SDESystem,higher_order_network_2))
for factor in [1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3]
    u0 = factor*rand(length(higher_order_network_1.states))
    p = factor*rand(length(higher_order_network_2.ps))
    t = rand()
    @test all(abs.(f1(u0,p,t) .- f2(u0,p,t)) .< 100*eps())
    @test all(abs.(f1.jac(u0,p,t) .- f2.jac(u0,p,t)) .< 100*eps())
    @test all(abs.(g1(u0,p,t) .- g2(u0,p,t)) .< 100*eps())
end

### Tests that the discrete jump systems are equal. ###
higher_order_network_3 = @reaction_network begin
    p,                                                  ∅ ⟼ X1
    r1*binomial(X1,2),                                  2X1 ⟾ 3X2
    mm(X1,r2,K)*binomial(X2,3),                         3X2 ⟾ X3 + 2X4
    r3*binomial(X3,1)*binomial(X4,2),                   X3 + 2X4 ⟾ 3X5 + 3X6
    r4*X2*binomial(X5,3)*binomial(X6,3),                3X5 + 3X6 ⟾ 3X5 + 2X7 + 4X8
    r5*binomial(X5,3)*binomial(X7,2)*binomial(X8,4),    3X5 + 2X7 + 4X8 ⟾ 10X9
    r6*binomial(X9,10),                                 10X9 ⟾ X10
    d*binomial(X10,2),                                  2X10 ⟾ ∅
end p r1 r2 K r3 r4 r5 r6 d

for factor in [1e-1, 1e0], repeat = 1:5
    u0 = rand(1:Int64(factor*100),length(higher_order_network_1.states))
    p = factor*rand(length(higher_order_network_3.ps))
    prob1 = JumpProblem(higher_order_network_1,DiscreteProblem(higher_order_network_1,u0,(0.,1000.),p),Direct())
    sol1 = solve(prob1,SSAStepper())
    prob2 = JumpProblem(higher_order_network_3,DiscreteProblem(higher_order_network_3,u0,(0.,1000.),p),Direct())
    sol2 = solve(prob2,SSAStepper())
    for i = 1:length(u0)
        vals1 = getindex.(sol1.u,i);
        vals2 = getindex.(sol1.u,i);
        (mean(vals2)>0.001) && @test 0.8 < mean(vals1)/mean(vals2) < 1.25
        (std(vals2)>0.001) && @test 0.8 < std(vals1)/std(vals2) < 1.25
    end
end
