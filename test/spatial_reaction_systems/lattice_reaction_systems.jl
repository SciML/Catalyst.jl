### Preparations ###

# Fetch packages.
using Catalyst, Graphs, Test

# Pre declares a grid.
grid = Graphs.grid([2, 2])

### Tests LatticeReactionSystem Getters Correctness ###

# Test case 1.
let 
    rs = @reaction_network begin
        (p, 1), 0 <--> X
    end
    tr = @transport_reaction d X    
    lrs = LatticeReactionSystem(rs, tr, grid)

    @test ModelingToolkit.getname.(species(lrs)) == [:X]   
    @test ModelingToolkit.getname.(spatial_species(lrs)) == [:X]   
    @test ModelingToolkit.getname.(parameters(lrs)) == [:p, :d]   
    @test ModelingToolkit.getname.(vertex_parameters(lrs)) == [:p]  
    @test ModelingToolkit.getname.(edge_parameters(lrs)) == [:d]      
end

# Test case 2.
let 
    rs = @reaction_network begin
        @parameters pX pY dX [edgeparameter=true] dY
        (pX, 1), 0 <--> X
        (pY, 1), 0 <--> Y
    end
    tr_1 = @transport_reaction dX X  
    tr_2 = @transport_reaction dY Y    
    lrs = LatticeReactionSystem(rs, [tr_1, tr_2], grid)

    @test ModelingToolkit.getname.(species(lrs)) == [:X, :Y]   
    @test ModelingToolkit.getname.(spatial_species(lrs)) == [:X, :Y]   
    @test ModelingToolkit.getname.(parameters(lrs)) == [:pX, :pY, :dX, :dY]   
    @test ModelingToolkit.getname.(vertex_parameters(lrs)) == [:pX, :pY, :dY]  
    @test ModelingToolkit.getname.(edge_parameters(lrs)) == [:dX]      
end

# Test case 3.
let 
    rs = @reaction_network begin
        @parameters dX p
        (pX, 1), 0 <--> X
        (pY, 1), 0 <--> Y
    end
    tr_1 = @transport_reaction dX X   
    lrs = LatticeReactionSystem(rs, tr_1, grid)

    @test ModelingToolkit.getname.(species(lrs)) == [:X, :Y]   
    @test ModelingToolkit.getname.(spatial_species(lrs)) == [:X]   
    @test ModelingToolkit.getname.(parameters(lrs)) == [:dX, :p, :pX, :pY]   
    @test ModelingToolkit.getname.(vertex_parameters(lrs)) == [:dX, :p, :pX, :pY]  
    @test ModelingToolkit.getname.(edge_parameters(lrs)) == []      
end

# Test case 4.
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
    lrs = LatticeReactionSystem(rs, [tr_1, tr_2, tr_3, tr_4, tr_5], grid)

    @test ModelingToolkit.getname.(species(lrs)) == [:W, :X, :Y, :Z, :V]   
    @test ModelingToolkit.getname.(spatial_species(lrs)) == [:X, :Y, :Z, :V, :W]   
    @test ModelingToolkit.getname.(parameters(lrs)) == [:pX, :pY, :dX, :dY, :pZ, :pV, :dZ, :dV, :dW]   
    @test ModelingToolkit.getname.(vertex_parameters(lrs)) == [:pX, :pY, :dY, :pZ, :pV]  
    @test ModelingToolkit.getname.(edge_parameters(lrs)) == [:dX, :dZ, :dV, :dW]      
end

### Tests Spatial Reactions Getters Correctness ###

# Test case 1.
let 
    tr_1 = @transport_reaction dX X    
    tr_2 = @transport_reaction dY1*dY2 Y     
    @test ModelingToolkit.getname.(species(tr_1)) == ModelingToolkit.getname.(spatial_species(tr_1)) == [:X]
    @test ModelingToolkit.getname.(species(tr_2)) == ModelingToolkit.getname.(spatial_species(tr_2)) == [:Y]
    @test ModelingToolkit.getname.(parameters(tr_1)) == [:dX]
    @test ModelingToolkit.getname.(parameters(tr_2)) == [:dY1, :dY2]
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
    @test isequal(species(tr_1), [X])
    @test isequal(species(tr_1), [X])
    @test isequal(spatial_species(tr_2), [Y])
    @test isequal(spatial_species(tr_2), [Y])
    @test isequal(parameters(tr_1), [dX])
    @test isequal(parameters(tr_2), [dY1, dY2])
end

### Test Interpolation ###

# Does not currently work. The 3 tr_macro_ lines generate errors.
# Test case 1.
# let
#     rs = @reaction_network begin
#         @species X(t) Y(t) Z(t)
#         @parameters dX dY1 dY2 dZ
#     end
#     @unpack X, Y, Z, dX, dY1, dY2, dZ = rs
#     rate1 = dX
#     rate2 = dY1*dY2 
#     species3 = Z
#     tr_1 = TransportReaction(dX, X)
#     tr_2 = TransportReaction(dY1*dY2, Y)
#     tr_2 = TransportReaction(dZ, Z)
#     tr_macro_1 = @transport_reaction $dX X
#     tr_macro_2 = @transport_reaction $(rate2) Y
#     tr_macro_3 = @transport_reaction dZ $species3
#     
#     @teest isequal(tr_1, tr_macro_1)
#     @teest isequal(tr_2, tr_macro_2)
#     @teest isequal(tr_3, tr_macro_3)
# end

### Tests Error generation ###

# Network where diffusion species is not declared in non-spatial network.
let 
    rs = @reaction_network begin
        (p, d), 0 <--> X
    end
    tr = @transport_reaction D Y    
    @test_throws ErrorException LatticeReactionSystem(rs, tr, grid)
end

# Network where the rate depend on a species
let 
    rs = @reaction_network begin
        @species Y(t)
        (p, d), 0 <--> X
    end
    tr = @transport_reaction D*Y X
    @test_throws ErrorException LatticeReactionSystem(rs, tr, grid)
end

# Network with edge parameter in non-spatial reaction rate.
let 
    rs = @reaction_network begin
        @parameters p [edgeparameter=true]
        (p, d), 0 <--> X
    end
    tr = @transport_reaction D X
    @test_throws ErrorException LatticeReactionSystem(rs, tr, grid)
end

# Network where metadata has been added in rs (which is not seen in transport reaction).
let 
    rs = @reaction_network begin
        @species X(t) [description="Species with added metadata"]
        (p, d), 0 <--> X
    end
    tr = @transport_reaction D X
    @test_throws ErrorException LatticeReactionSystem(rs, tr, grid)

    rs = @reaction_network begin
        @parameters D [description="Parameter with added metadata"]
        (p, d), 0 <--> X
    end
    tr = @transport_reaction D X
    @test_throws ErrorException LatticeReactionSystem(rs, tr, grid)
end

