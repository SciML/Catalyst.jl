function testnetwork(rn, rn2)
    rxids = 1:numreactions(rn)
    rxids2 = 1:numreactions(rn2)
    @test params(rn) == params(rn2)
    @test paramsmap(rn) == paramsmap(rn2)
    @test species(rn) == species(rn2)
    @test speciesmap(rn) == speciesmap(rn2)
    @test dependents.(rn,rxids) == dependents.(rn2,rxids2)
    @test substrates.(rn,rxids) == substrates.(rn2,rxids2)
    @test products.(rn,rxids) == products.(rn2,rxids2)
    @test rateexpr.(rn,rxids) == rateexpr.(rn2,rxids2)
    @test oderatelawexpr.(rn,rxids) == oderatelawexpr.(rn2,rxids2)
    @test ssaratelawexpr.(rn,rxids) == ssaratelawexpr.(rn2,rxids2)
    @test ismassaction.(rn,rxids) == ismassaction.(rn2,rxids2)
end

# baseline network to test against
repressilator = @min_reaction_network begin
    hillr(P₃,α,K,n), ∅ --> m₁
    hillr(P₁,α,K,n), ∅ --> m₃
    (δ,γ), m₁ ↔ ∅
    (δ,γ), m₃ ↔ ∅
    β, m₁ --> m₁ + P₁
    10., m₃ --> m₃ + P₃
    2, P₁ --> ∅
    β*P₂, 2m₁ + 3P₁ --> P₂
end α K n δ γ β μ;
numrxs = numreactions(repressilator)

# empty network
rmin = @empty_reaction_network
foreach(param -> addparam!(rmin, param), params(repressilator))
foreach(spec -> addspecies!(rmin, spec), species(repressilator))
addreaction!(rmin, :(hillr(P₃,α,K,n)), :(∅ --> m₁) )
addreaction!(rmin, :(hillr(P₁,α,K,n)), :(∅ --> m₃) )
addreaction!(rmin, :((δ,γ)), :(m₁ ↔ ∅) )
addreaction!(rmin, :((δ,γ)), :(m₃ ↔ ∅) )
addreaction!(rmin, :β, :(m₁ --> m₁ + P₁) )
addreaction!(rmin, 10., :(m₃ --> m₃ + P₃) )
addreaction!(rmin, 2, :(P₁ --> ∅) )
addreaction!(rmin, :(β*P₂), :(2m₁ + 3P₁ --> P₂) )  
testnetwork(rmin, repressilator)

# empty network stoichiometry interface
rmin = @empty_reaction_network
foreach(param -> addparam!(rmin, param), params(repressilator))
foreach(spec -> addspecies!(rmin, spec), species(repressilator))
addreaction!(rmin, :(hillr(P₃,α,K,n)), (), (:m₁=>1,))
addreaction!(rmin, :(hillr(P₁,α,K,n)), (), (:m₃=>1,) )
addreaction!(rmin, :δ, (:m₁ => 1,), ())
addreaction!(rmin, :γ, (), (:m₁ => 1,))
addreaction!(rmin, :δ, (:m₃ => 1,), ())
addreaction!(rmin, :γ, (), (:m₃ => 1,))
addreaction!(rmin, :β, (:m₁=>1,), (:m₁=>1, :P₁=>1) )
addreaction!(rmin, 10., (:m₃=>1,), (:m₃=>1, :P₃=>1) )
addreaction!(rmin, 2, (:P₁=>1,),() )
addreaction!(rmin, :(β*P₂), (:m₁=>2,:P₁=>3), (:P₂=>1,) )  
testnetwork(rmin, repressilator)


# partial min_reaction_network
rmin2 = @min_reaction_network begin
    hillr(P₃,α,K,n), ∅ --> m₁
    hillr(P₁,α,K,n), ∅ --> m₃
    (δ,γ), m₁ ↔ ∅
    (δ,γ), m₃ ↔ ∅
    β, m₁ --> m₁ + P₁
    10., m₃ --> m₃ + P₃
end α K n δ γ β;
foreach(param -> addparam!(rmin2, param), params(repressilator))
foreach(spec -> addspecies!(rmin2, spec), species(repressilator))
addreaction!(rmin2, 2, :(P₁ --> ∅) )
addreaction!(rmin2, :(β*P₂), :(2m₁ + 3P₁ --> P₂) )  
testnetwork(rmin2, repressilator)

# partial min_reaction_network, stoichiometry interface
rmin2 = @min_reaction_network begin
    hillr(P₃,α,K,n), ∅ --> m₁
    hillr(P₁,α,K,n), ∅ --> m₃
    (δ,γ), m₁ ↔ ∅
    (δ,γ), m₃ ↔ ∅
    β, m₁ --> m₁ + P₁
    10., m₃ --> m₃ + P₃
end α K n δ γ β;
foreach(param -> addparam!(rmin2, param), params(repressilator))
foreach(spec -> addspecies!(rmin2, spec), species(repressilator))
addreaction!(rmin2, 2, (:P₁=>1,), ())
addreaction!(rmin2, :(β*P₂), (:m₁=>2,:P₁=>3), (:P₂=>1,) )  
testnetwork(rmin2, repressilator)


# simple example to test numerical agreement of solutions
fullrn = @reaction_network begin
    (k*X*(X-1)*(X-2),Y * X), X + 2X ⟺ Y+ X
end k
p = (1/6,)
u0 = [1000.,1000.]
tspan = (0.,5e-4)
foprob = ODEProblem(fullrn, u0, tspan, p)
fsol = solve(foprob,Tsit5(), abstol=1e-10, reltol=1e-10)

emptyrn = @empty_reaction_network
addspecies!(emptyrn, :X)
addspecies!(emptyrn, "Y")
addparam!(emptyrn, :k)
addreaction!(emptyrn, :((k*X*(X-1)*(X-2),Y * X)), :(X + 2X ⟺ Y+ X))
addodes!(emptyrn)
eoprob = ODEProblem(emptyrn, u0, tspan, p)
esol = solve(eoprob,Tsit5(), abstol=1e-10, reltol=1e-10)
@test norm(esol .- fsol, Inf) < 1e-7

minrn = @min_reaction_network begin
    X*Y, X + Y ⟾ X + 2X
end
addparam!(minrn, "k")
addreaction!(minrn, :(k*X*(X-1)*(X-2)), :(X + 2X ⟾ Y + X))
addodes!(minrn)
moprob = ODEProblem(minrn, u0, tspan, p)
msol = solve(moprob,Tsit5(), abstol=1e-10, reltol=1e-10)
@test norm(msol .- fsol, Inf) < 1e-7


# another simple example to test numerical agreement of solutions
# using stoichiometry interface
fullrn = @reaction_network begin
    (X^2/6,1.0), X ↔ Y + X
end
u0 = [1000.,1000.]
tspan = (0.,5e-4)
foprob = ODEProblem(fullrn, u0, tspan)
fsol = solve(foprob,Tsit5(), abstol=1e-10, reltol=1e-10)

emptyrn = @empty_reaction_network
addspecies!(emptyrn, :X)
addspecies!(emptyrn, "Y")
addreaction!(emptyrn, :(X^2/6), (:X=>1,), (:X=>1, :Y=>1))
addreaction!(emptyrn, 1.0, (:X=>1,:Y=>1), (:X=>1,))
addodes!(emptyrn)
eoprob = ODEProblem(emptyrn, u0, tspan)
esol = solve(eoprob,Tsit5(), abstol=1e-10, reltol=1e-10)
@test norm(esol .- fsol, Inf) < 1e-7
