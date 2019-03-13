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
rxids = 1:numrxs

function testnetwork(rn)
    rxidsmin = 1:numreactions(rn)
    @test params(rn) == params(repressilator)
    @test paramsmap(rn) == paramsmap(repressilator)
    @test species(rn) == species(repressilator)
    @test speciesmap(rn) == speciesmap(repressilator)
    @test dependents.(rn,rxidsmin) == dependents.(repressilator,rxids)
    @test substrates.(rn,rxidsmin) == substrates.(repressilator,rxids)
    @test products.(rn,rxidsmin) == products.(repressilator,rxids)
    @test rateexpr.(rn,rxidsmin) == rateexpr.(repressilator,rxids)
    @test oderatelawexpr.(rn,rxidsmin) == oderatelawexpr.(repressilator,rxids)
    @test ssaratelawexpr.(rn,rxidsmin) == ssaratelawexpr.(repressilator,rxids)
    @test ismassaction.(rn,rxidsmin) == ismassaction.(repressilator,rxids)
end

# empty network
rmin = @empty_reaction_network

# parameters
foreach(param -> addparam!(rmin, param), params(repressilator))

# species
foreach(spec -> addspecies!(rmin, spec), species(repressilator))

# reactions
addreaction!(rmin, :(hillr(P₃,α,K,n)), :(∅ --> m₁) )
addreaction!(rmin, :(hillr(P₁,α,K,n)), :(∅ --> m₃) )
addreaction!(rmin, :((δ,γ)), :(m₁ ↔ ∅) )
addreaction!(rmin, :((δ,γ)), :(m₃ ↔ ∅) )
addreaction!(rmin, :β, :(m₁ --> m₁ + P₁) )
addreaction!(rmin, 10., :(m₃ --> m₃ + P₃) )
addreaction!(rmin, 2, :(P₁ --> ∅) )
addreaction!(rmin, :(β*P₂), :(2m₁ + 3P₁ --> P₂) )  

testnetwork(rmin)

# partial min_reaction_network

rmin2 = @min_reaction_network begin
    hillr(P₃,α,K,n), ∅ --> m₁
    hillr(P₁,α,K,n), ∅ --> m₃
    (δ,γ), m₁ ↔ ∅
    (δ,γ), m₃ ↔ ∅
    β, m₁ --> m₁ + P₁
    10., m₃ --> m₃ + P₃
end α K n δ γ β;

# params, should only add missing ones
foreach(param -> addparam!(rmin2, param), params(repressilator))

# species, should only add missing ones
foreach(spec -> addspecies!(rmin2, spec), species(repressilator))

# reactions
addreaction!(rmin2, 2, :(P₁ --> ∅) )
addreaction!(rmin2, :(β*P₂), :(2m₁ + 3P₁ --> P₂) )  

testnetwork(rmin2)


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
