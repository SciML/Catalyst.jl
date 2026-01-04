### Preparations ###

# Fetch packages.
using Catalyst, Test


### TransportReaction Creation Tests ###

# Tests TransportReaction with non-trivial rate.
let
    rs = @reaction_network begin
        @parameters dV dE [edgeparameter = true]
        (p, 1), 0 <--> X
    end
    @unpack dV, dE, X = rs

    tr = TransportReaction(dV * dE, X)
    @test isequal(tr.rate, dV * dE)
end
# Test reactions with constants in rate.
let
    t = default_t()
    @species X(t) Y(t)

    tr_1 = TransportReaction(1.5, X)
    tr_1_macro = @transport_reaction 1.5 X
    @test isequal(tr_1.rate, tr_1_macro.rate)
    @test isequal(tr_1.species, tr_1_macro.species)

    tr_2 = TransportReaction(π, Y)
    tr_2_macro = @transport_reaction π Y
    @test isequal(tr_2.rate, tr_2_macro.rate)
    @test isequal(tr_2.species, tr_2_macro.species)
end

### Spatial Reactions Getters Correctness ###

# Test case 1.
let
    tr_1 = @transport_reaction dX X
    tr_2 = @transport_reaction dY1 * dY2 Y

    # @test ModelingToolkit.getname.(species(tr_1)) == ModelingToolkit.getname.(spatial_species(tr_1)) == [:X] # species(::TransportReaction) currently not supported.
    # @test ModelingToolkit.getname.(species(tr_2)) == ModelingToolkit.getname.(spatial_species(tr_2)) == [:Y]
    @test ModelingToolkit.getname.(spatial_species(tr_1)) == [:X]
    @test ModelingToolkit.getname.(spatial_species(tr_2)) == [:Y]
    @test ModelingToolkit.getname.(parameters(tr_1)) == [:dX]
    @test ModelingToolkit.getname.(parameters(tr_2)) == [:dY1, :dY2]

    # @test issetequal(species(tr_1), [tr_1.species])
    # @test issetequal(species(tr_2), [tr_2.species])
    @test issetequal(spatial_species(tr_1), [tr_1.species])
    @test issetequal(spatial_species(tr_2), [tr_2.species])
end

# Test case 2.
let
    rs = @reaction_network begin
        @species X(t) Y(t)
        @parameters dX dY1 dY2
    end
    @unpack X, Y, dX, dY1, dY2 = rs
    tr_1 = TransportReaction(dX, X)
    tr_2 = TransportReaction(dY1 * dY2, Y)
    # @test isequal(species(tr_1), [X])
    # @test isequal(species(tr_1), [X])
    @test issetequal(spatial_species(tr_2), [Y])
    @test issetequal(spatial_species(tr_2), [Y])
    @test issetequal(parameters(tr_1), [dX])
    @test issetequal(parameters(tr_2), [dY1, dY2])
end

### Error Tests ###

# Tests that creation of TransportReaction with non-parameters in rate yield errors.
# Tests that errors are throw even when the rate is highly nested.
# Tests that error is thrown when non-permitted rate symbol is used.
let
    t = default_t()
    @species X(t) Y(t)
    @parameters D1 D2 D3
    @test_throws ErrorException TransportReaction(D1 + D2 * (D3 + Y), X)
    @test_throws ErrorException TransportReaction(Y, X)
    @test_throws Exception @eval @transport_reaction ∅ X
end

### Other Tests ###

# Test Interpolation
# Does not currently work. The 3 tr_macro_ lines generate errors.
let
    rs = @reaction_network begin
        @species X(t) Y(t) Z(t)
        @parameters dX dY1 dY2 dZ
    end
    @unpack X, Y, Z, dX, dY1, dY2, dZ = rs
    rate1 = dX
    rate2 = dY1 * dY2
    species3 = Z
    tr_1 = TransportReaction(dX, X)
    tr_2 = TransportReaction(dY1 * dY2, Y)
    tr_3 = TransportReaction(dZ, Z)
    tr_macro_1 = @transport_reaction $dX X
    tr_macro_2 = @transport_reaction $(rate2) Y
    @test_broken false
    # tr_macro_3 = @transport_reaction dZ $species3 # Currently does not work, something with meta programming.

    @test isequal(tr_1, tr_macro_1)
    @test isequal(tr_2, tr_macro_2)
    # @test isequal(tr_3, tr_macro_3)
end

# Checks that the `hash` functions work for `TransportReaction`s.
let
    tr1 = @transport_reaction D1 X
    tr2 = @transport_reaction D1 X
    tr3 = @transport_reaction D2 X
    hash(tr1, 0x0000000000000001) == hash(tr2, 0x0000000000000001)
    hash(tr2, 0x0000000000000001) != hash(tr3, 0x0000000000000001)
end
