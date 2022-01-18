using Catalyst, ModelingToolkit, OrdinaryDiffEq, LinearAlgebra

# naming tests
@parameters k
@variables t, A(t)
rx = Reaction(k, [A], nothing)
function rntest(rn, name)
    @test nameof(rn) == name
    @test isequal(species(rn)[1], ModelingToolkit.unwrap(A))
    @test isequal(parameters(rn)[1], ModelingToolkit.unwrap(k))
    @test reactions(rn)[1] == rx
end

function emptyrntest(rn, name)
    @test nameof(rn) == name
    @test numreactions(rn) == 0
    @test numspecies(rn) == 0
    @test numreactionparams(rn) == 0
end

rn = @reaction_network name begin
    k, A --> 0
end k
rntest(rn, :name)

name = :blah
rn = @reaction_network $name begin
    k, A --> 0
end k
rntest(rn, :blah)

rn = @reaction_network begin
    k, A --> 0
end k
rntest(rn, nameof(rn))

function makern(; name)
    @reaction_network $name begin
        k, A --> 0
    end k
end
@named testnet = makern()
rntest(testnet, :testnet)

rn = @reaction_network name
emptyrntest(rn, :name)

rn = @reaction_network $name
emptyrntest(rn, :blah)

# test variables that appear only in rates and aren't ps
# are categorized as species
rn = @reaction_network begin
    π*k*D*hill(B,k2,B*D*H,n), 3*A  --> 2*C
end k k2 n
@parameters k,k2,n
@variables t,A(t),B(t),C(t),D(t),H(t)
@test issetequal([A,B,C,D,H], species(rn))
@test issetequal([k,k2,n], parameters(rn))

# test interpolation within the DSL
@parameters k, k1, k2
@variables t, A(t), B(t), C(t), D(t)
AA = A
AAA = A^2 + B
rn = @reaction_network rn begin
    k*$AAA, C --> D
end k
rn2 = ReactionSystem([Reaction(k*AAA, [C], [D])], t; name=:rn)
@test rn == rn2

rn = @reaction_network rn begin
    k, $AA + C --> D
end k
rn2 = ReactionSystem([Reaction(k, [AA,C], [D])], t; name=:rn)
@test rn == rn2

BB = B; A2 = A
rn = @reaction_network rn begin
    (k1,k2), C + $A2 + $BB + $A2 <--> $BB + $BB
end k1 k2
rn2 = ReactionSystem([Reaction(k1, [C, A, B], [B], [1,2,1],[2]),
                      Reaction(k2, [B], [C, A, B], [2], [1,2,1])],
                     t; name=:rn)
@test rn == rn2

kk1 = k^2*A
kk2 = k1+k2
rn = @reaction_network rn begin
    α+$kk1*$kk2*$AA, 2*$AA + B --> $AA
end α
@parameters α
rn2 = ReactionSystem([Reaction(α+kk1*kk2*AA, [A, B], [A], [2, 1], [1])], t; name=:rn)
@test rn == rn2

@testset "make_reaction_system can be called from another module" begin
    ex = quote
        (Ka, Depot --> Central)
        (CL / Vc, Central --> 0)
    end
    # Line number nodes aren't ignored so have to be manually removed
    Base.remove_linenums!(ex)
    @test eval(Catalyst.make_reaction_system(ex, (:Ka, :Cl, :Vc))) isa ReactionSystem
end

# test defaults
rn = @reaction_network begin
    α, S + I --> 2I
    β, I --> R
end α β
p     = [.1/1000, .01]
tspan = (0.0,250.0)
u0    = [999.0,1.0,0.0]
op    = ODEProblem(rn, species(rn) .=> u0, tspan, parameters(rn) .=> p)
sol   = solve(op, Tsit5())  # old style 
setdefaults!(rn, [:S => 999.0, :I => 1.0, :R => 0.0, :α => 1e-4, :β => .01])
op = ODEProblem(rn, [], tspan, [])
sol2 = solve(op, Tsit5())
@test norm(sol.u - sol2.u) ≈ 0
