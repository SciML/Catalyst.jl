using Catalyst, ModelingToolkit

# Test creating networks with/without options.
@reaction_network begin
    (k1,k2), A <--> B
end
@reaction_network begin
@parameters k1 k2
    (k1,k2), A <--> B
end
@reaction_network begin
@parameters k1 k2
@species A B
    (k1,k2), A <--> B
end
@reaction_network begin
@species A B
    (k1,k2), A <--> B
end

@reaction_network begin
@parameters begin
    k1 
    k2
end
    (k1,k2), A <--> B
end
@reaction_network begin
@species begin
    A
    B
end
    (k1,k2), A <--> B
end
@reaction_network begin
@parameters begin
    k1 
    k2
end
@species begin
    A
    B
end
    (k1,k2), A <--> B
end

@reaction_network name begin
    (k1,k2), A <--> B
end
@reaction_network name begin
@parameters k1 k2
    (k1,k2), A <--> B
end
@reaction_network name begin
@species A B
    (k1,k2), A <--> B
end
@reaction_network name begin
@parameters k1 k2
@species A B
    (k1,k2), A <--> B
end

@reaction_network name begin
    (k1,k2), A <--> B
end
@reaction_network name begin
    (k1,k2), A <--> B
@parameters k1 k2
end
@reaction_network name begin
    (k1,k2), A <--> B
@species A B
end
@reaction_network name begin
    (k1,k2), A <--> B
@parameters k1 k2
@species A B
end


@reaction_network name begin
@parameters begin
    k1 
    k2
end
    (k1,k2), A <--> B
end
@reaction_network name begin
@species begin
    A
    B
end
    (k1,k2), A <--> B
end
@reaction_network name begin
@parameters begin
    k1 
    k2
end
@species begin
    A
    B
end
    (k1,k2), A <--> B
end

# Checks that some created networks are identical.
rn1 = @reaction_network name begin
    (k1,k2), A <--> B
end
rn2 = @reaction_network name begin
@parameters k1 k2
    (k1,k2), A <--> B
end
rn3 = @reaction_network name begin
@species A B
    (k1,k2), A <--> B
end
rn4 = @reaction_network name begin
@parameters k1 k2
@species A B
    (k1,k2), A <--> B
end


