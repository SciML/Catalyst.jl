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

