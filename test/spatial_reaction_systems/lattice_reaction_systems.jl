### Preparations ###

# Fetch packages.
using Catalyst, Graphs, Test

# Fetch test networks.
include("../spatial_test_networks.jl")

# Pre declares a grid.
grids = [very_small_2d_cartesian_grid, very_small_2d_masked_grid, very_small_2d_graph_grid]

### Tests LatticeReactionSystem Getters Correctness ###

# Test case 1.
let 
    rs = @reaction_network begin
        (p, 1), 0 <--> X
    end
    tr = @transport_reaction d X   
    for grid in grids 
        lrs = LatticeReactionSystem(rs, [tr], grid)

        @unpack X, p = rs
        d = edge_parameters(lrs)[1]
        @test issetequal(species(lrs), [X])  
        @test issetequal(spatial_species(lrs), [X])   
        @test issetequal(parameters(lrs), [p, d])   
        @test issetequal(vertex_parameters(lrs), [p])  
        @test issetequal(edge_parameters(lrs), [d])      
    end
end

# Test case 2.
let 
    rs = @reaction_network begin
        @parameters p1 p2 [edgeparameter=true]
    end
    @unpack p1, p2 = rs

    @test !isedgeparameter(p1)
    @test isedgeparameter(p2)
end

# Test case 3.
let 
    rs = @reaction_network begin
        @parameters pX pY dX [edgeparameter=true] dY
        (pX, 1), 0 <--> X
        (pY, 1), 0 <--> Y
    end
    tr_1 = @transport_reaction dX X  
    tr_2 = @transport_reaction dY Y    
    for grid in grids
        lrs = LatticeReactionSystem(rs, [tr_1, tr_2], grid)

        @unpack X, Y, pX, pY, dX, dY = rs
        @test issetequal(species(lrs), [X, Y])
        @test issetequal(spatial_species(lrs), [X, Y])
        @test issetequal(parameters(lrs), [pX, pY, dX, dY])
        @test issetequal(vertex_parameters(lrs), [pX, pY, dY])
        @test issetequal(edge_parameters(lrs), [dX])
    end
end

# Test case 4.
let 
    rs = @reaction_network begin
        @parameters dX p
        (pX, 1), 0 <--> X
        (pY, 1), 0 <--> Y
    end
    tr_1 = @transport_reaction dX X   
    for grid in grids
        lrs = LatticeReactionSystem(rs, [tr_1], grid)

        @unpack dX, p, X, Y, pX, pY = rs
        @test issetequal(species(lrs), [X, Y])
        @test issetequal(spatial_species(lrs), [X])
        @test issetequal(parameters(lrs), [dX, p, pX, pY])
        @test issetequal(vertex_parameters(lrs), [dX, p, pX, pY])
        @test issetequal(edge_parameters(lrs), [])
    end
end

# Test case 5.
let 
    rs = @reaction_network begin
        @species W(t)
        @parameters pX pY dX [edgeparameter=true] dY
        (pX, 1), 0 <--> X
        (pY, 1), 0 <--> Y
        (pZ, 1), 0 <--> Z
        (pV, 1), 0 <--> V
    end
    @unpack dX, X, V = rs
    @parameters dV dW
    @variables t
    @species W(t)
    tr_1 = TransportReaction(dX, X)  
    tr_2 = @transport_reaction dY Y    
    tr_3 = @transport_reaction dZ Z 
    tr_4 = TransportReaction(dV, V)  
    tr_5 = TransportReaction(dW, W)   
    for grid in grids  
        lrs = LatticeReactionSystem(rs, [tr_1, tr_2, tr_3, tr_4, tr_5], grid)

        @unpack pX, pY, pZ, pV, dX, dY, X, Y, Z, V = rs
        dZ, dV, dW = edge_parameters(lrs)[2:end]
        @test issetequal(species(lrs), [W, X, Y, Z, V])
        @test issetequal(spatial_species(lrs), [X, Y, Z, V, W])
        @test issetequal(parameters(lrs), [pX, pY, dX, dY, pZ, pV, dZ, dV, dW])
        @test issetequal(vertex_parameters(lrs), [pX, pY, dY, pZ, pV])
        @test issetequal(edge_parameters(lrs), [dX, dZ, dV, dW])
    end
end

# Test case 6.
let 
    rs = @reaction_network customname begin
        (p, 1), 0 <--> X
    end
    tr = @transport_reaction d X    
    for grid in grids
        lrs = LatticeReactionSystem(rs, [tr], grid)

        @test nameof(lrs) == :customname
    end
end

### Tests Spatial Reactions Getters Correctness ###

# Test case 1.
let 
    tr_1 = @transport_reaction dX X    
    tr_2 = @transport_reaction dY1*dY2 Y   

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
    tr_2 = TransportReaction(dY1*dY2, Y)
    # @test isequal(species(tr_1), [X])
    # @test isequal(species(tr_1), [X])
    @test issetequal(spatial_species(tr_2), [Y])
    @test issetequal(spatial_species(tr_2), [Y])
    @test issetequal(parameters(tr_1), [dX])
    @test issetequal(parameters(tr_2), [dY1, dY2])
end

### Tests Spatial Reactions Generation ###

# Tests TransportReaction with non-trivial rate.
let 
    rs = @reaction_network begin
        @parameters dV dE [edgeparameter=true] 
        (p,1), 0 <--> X
    end
    @unpack dV, dE, X = rs
    
    tr = TransportReaction(dV*dE, X)
    @test isequal(tr.rate, dV*dE)
end

# Tests transport_reactions function for creating TransportReactions.
let 
    rs = @reaction_network begin
        @parameters d
        (p,1), 0 <--> X
    end
    @unpack d, X = rs
    trs = TransportReactions([(d, X), (d, X)])
    @test isequal(trs[1], trs[2])
end

# Test reactions with constants in rate.
let 
    @variables t
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

### Test Interpolation ###

# Does not currently work. The 3 tr_macro_ lines generate errors.
# Test case 1.
let
    rs = @reaction_network begin
        @species X(t) Y(t) Z(t)
        @parameters dX dY1 dY2 dZ
    end
    @unpack X, Y, Z, dX, dY1, dY2, dZ = rs
    rate1 = dX
    rate2 = dY1*dY2 
    species3 = Z
    tr_1 = TransportReaction(dX, X)
    tr_2 = TransportReaction(dY1*dY2, Y)
    tr_3 = TransportReaction(dZ, Z)
    tr_macro_1 = @transport_reaction $dX X
    tr_macro_2 = @transport_reaction $(rate2) Y
    # tr_macro_3 = @transport_reaction dZ $species3 # Currently does not work, something with meta programming.
    
    @test isequal(tr_1, tr_macro_1)
    @test isequal(tr_2, tr_macro_2) # Unsure why these fails, since for components equality hold: `isequal(tr_1.species, tr_macro_1.species)` and `isequal(tr_1.rate, tr_macro_1.rate)` are both true.
    # @test isequal(tr_3, tr_macro_3)
end

### Tests Error generation ###

# Test creation of TransportReaction with non-parameters in rate.
# Tests that it works even when rate is highly nested.
let 
    @variables t
    @species X(t) Y(t)
    @parameters D1 D2 D3
    @test_throws ErrorException TransportReaction(D1 + D2*(D3 + Y), X)
    @test_throws ErrorException TransportReaction(Y, X)
end

# Network where diffusion species is not declared in non-spatial network.
let 
    rs = @reaction_network begin
        (p, d), 0 <--> X
    end
    tr = @transport_reaction D Y    
    for grid in grids
        @test_throws ErrorException LatticeReactionSystem(rs, [tr], grid)
    end
end

# Network where the rate depend on a species
let 
    rs = @reaction_network begin
        @species Y(t)
        (p, d), 0 <--> X
    end
    tr = @transport_reaction D*Y X
    for grid in grids
        @test_throws ErrorException LatticeReactionSystem(rs, [tr], grid)
    end
end

# Network with edge parameter in non-spatial reaction rate.
let 
    rs = @reaction_network begin
        @parameters p [edgeparameter=true]
        (p, d), 0 <--> X
    end
    tr = @transport_reaction D X
    for grid in grids
        @test_throws ErrorException LatticeReactionSystem(rs, [tr], grid)
    end
end

# Network where metadata has been added in rs (which is not seen in transport reaction).
let 
    rs = @reaction_network begin
        @species X(t) [description="Species with added metadata"]
        (p, d), 0 <--> X
    end
    tr = @transport_reaction D X
    for grid in grids
        @test_throws ErrorException LatticeReactionSystem(rs, [tr], grid)

        rs = @reaction_network begin
            @parameters D [description="Parameter with added metadata"]
            (p, d), 0 <--> X
        end
        tr = @transport_reaction D X
        @test_throws ErrorException LatticeReactionSystem(rs, [tr], grid)
    end
end


### Tests Grid Vertex and Edge Number Computation ###

# Tests that the correct numbers are computed for num_edges.
let 
    # Function counting the values in an iterator by stepping through it.
    function iterator_count(iterator)
        count = 0
        foreach(e -> count+=1, iterator)
        return count
    end

    # Cartesian and masked grid (test diagonal edges as well).
    for lattice in [small_1d_cartesian_grid, small_2d_cartesian_grid, small_3d_cartesian_grid, 
                random_1d_masked_grid, random_2d_masked_grid, random_3d_masked_grid]
        lrs1 = LatticeReactionSystem(SIR_system, SIR_srs_1, lattice)
        lrs2 = LatticeReactionSystem(SIR_system, SIR_srs_1, lattice; diagonal_connections=true)
        @test lrs1.num_edges == iterator_count(lrs1.edge_iterator)    
        @test lrs2.num_edges == iterator_count(lrs2.edge_iterator)    
    end

    # Graph grids (cannot test diagonal connections).
    for lattice in [small_2d_graph_grid, small_3d_graph_grid, undirected_cycle, small_directed_cycle, unconnected_graph]
        lrs1 = LatticeReactionSystem(SIR_system, SIR_srs_1, lattice)
        @test lrs1.num_edges == iterator_count(lrs1.edge_iterator)    
    end
end

