#! format: off

using Catalyst, ModelingToolkit

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
    @parameters k
    k, A --> 0
end
rntest(rn, :name)

name = :blah
rn = @reaction_network $name begin
    @parameters k
    k, A --> 0
end
rntest(rn, :blah)

rn = @reaction_network begin
    @parameters k
    k, A --> 0
end
rntest(rn, nameof(rn))

function makern(; name)
    @reaction_network $name begin
    @parameters k
        k, A --> 0
    end
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
    @parameters k k2 n
    @species A(t) B(t) C(t) D(t) H(t)
    k*D*hill(B,k2,B*D*H,n), 3*A  --> 2*C
end
@parameters k k2 n
@variables t A(t) B(t) C(t) D(t) H(t)
@test issetequal([A,B,C,D,H], species(rn))
@test issetequal([k,k2,n], parameters(rn))

# test interpolation within the DSL
@parameters α k k1 k2
@variables t A(t) B(t) C(t) D(t)
AA = A
AAA = A^2 + B
rn = @reaction_network rn begin
    @parameters k
    @species A(t) B(t) C(t) D(t)
    k*$AAA, C --> D
end
rn2 = ReactionSystem([Reaction(k*AAA, [C], [D])], t; name=:rn)
@test rn == rn2

rn = @reaction_network rn begin
    @parameters k
    @species A(t) C(t) D(t)
    k, $AA + C --> D
end
rn2 = ReactionSystem([Reaction(k, [AA,C], [D])], t; name=:rn)
@test rn == rn2

BB = B; A2 = A
rn = @reaction_network rn begin
    @parameters k1 k2
    (k1,k2), C + $A2 + $BB + $A2 <--> $BB + $BB
end
rn2 = ReactionSystem([Reaction(k1, [C, A, B], [B], [1,2,1],[2]),
                      Reaction(k2, [B], [C, A, B], [2], [1,2,1])],
                     t; name=:rn)
@test rn == rn2

kk1 = k^2*A
kk2 = k1+k2
rn = @reaction_network rn begin
    @parameters α k k1 k2
    α+$kk1*$kk2*$AA, 2*$AA + B --> $AA
end
rn2 = ReactionSystem([Reaction(α+kk1*kk2*AA, [A, B], [A], [2, 1], [1])], t; name=:rn)
@test rn == rn2

@testset "make_reaction_system can be called from another module" begin
    ex = quote
        (Ka, Depot --> Central)
        (CL / Vc, Central --> 0)
    end
    # Line number nodes aren't ignored so have to be manually removed
    Base.remove_linenums!(ex)
    @test eval(Catalyst.make_reaction_system(ex)) isa ReactionSystem
end

rx = @reaction k*h, A + 2*B --> 3*C + D
@parameters k h
@variables t A(t) B(t) C(t) D(t)
@test rx == Reaction(k*h,[A,B],[C,D],[1,2],[3,1])

ex = k*A^2 + B
V  = A
rx = @reaction b+$ex, 2*$V + C--> ∅
@parameters b
@test rx == Reaction(b+ex, [A,C], nothing, [2,1], nothing)

### test floating point stoichiometry work ###
@parameters k
@variables t B(t) C(t) D(t)
rx1 = Reaction(k,[B,C],[B,D], [2.5,1],[3.5, 2.5])
rx2 = Reaction(2*k, [B], [D], [1], [2.5])
rx3 = Reaction(2*k, [B], [D], [2.5], [2])
@named mixedsys = ReactionSystem([rx1,rx2,rx3],t,[B,C,D],[k])
osys = convert(ODESystem, mixedsys; combinatoric_ratelaws=false)
rn = @reaction_network mixedsys begin
    @parameters k
    k, 2.5*B + C --> 3.5*B + 2.5*D
    2*k, B --> 2.5*D
    2*k, 2.5*B --> 2*D
end
@test rn == mixedsys
