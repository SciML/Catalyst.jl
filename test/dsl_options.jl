using Catalyst, ModelingToolkit

# Test creating networks with/without options.
@reaction_network begin (k1, k2), A <--> B end
@reaction_network begin
    @parameters k1 k2
    (k1, k2), A <--> B
end
@reaction_network begin
    @parameters k1 k2
    @species A B
    (k1, k2), A <--> B
end
@reaction_network begin
    @species A B
    (k1, k2), A <--> B
end

@reaction_network begin
    @parameters begin
        k1
        k2
    end
    (k1, k2), A <--> B
end
@reaction_network begin
    @species begin
        A
        B
    end
    (k1, k2), A <--> B
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
    (k1, k2), A <--> B
end

@reaction_network name begin (k1, k2), A <--> B end
@reaction_network name begin
    @parameters k1 k2
    (k1, k2), A <--> B
end
@reaction_network name begin
    @species A B
    (k1, k2), A <--> B
end
@reaction_network name begin
    @parameters k1 k2
    @species A B
    (k1, k2), A <--> B
end

@reaction_network name begin (k1, k2), A <--> B end
@reaction_network name begin
    (k1, k2), A <--> B
    @parameters k1 k2
end
@reaction_network name begin
    (k1, k2), A <--> B
    @species A B
end
@reaction_network name begin
    (k1, k2), A <--> B
    @parameters k1 k2
    @species A B
end

@reaction_network name begin
    @parameters begin
        k1
        k2
    end
    (k1, k2), A <--> B
end
@reaction_network name begin
    @species begin
        A
        B
    end
    (k1, k2), A <--> B
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
    (k1, k2), A <--> B
end

# Checks that some created networks are identical.
rn1 = @reaction_network name begin (k1, k2), A <--> B end
rn2 = @reaction_network name begin
    @parameters k1 k2
    (k1, k2), A <--> B
end
rn3 = @reaction_network name begin
    @species A B
    (k1, k2), A <--> B
end
rn4 = @reaction_network name begin
    @parameters k1 k2
    @species A B
    (k1, k2), A <--> B
end
#@test isequal(species(rn1),species(rn2))
#@test isequal(species(rn2),species(rn3))
#@test isequal(species(rn3),species(rn4))
@test isequal(parameters(rn1), parameters(rn2))
@test isequal(parameters(rn2), parameters(rn3))
@test isequal(parameters(rn3), parameters(rn4))

parameters(rn1)

rn5 = @reaction_network name begin (k1, k2), A <--> B end
rn6 = @reaction_network name begin
    @parameters k2 k1
    @species B A
    (k1, k2), A <--> B
end
#@test !isequal(species(rn5),species(rn6)) # This does not work, as ReactionSystem does not accept ordering of species, but does that itself.
@test !isequal(parameters(rn5), parameters(rn6))


# Checks that the rights things are put in vectors
rn7 = @reaction_network name begin
    @parameters p d1 d2
    @species A B
        p, 0 --> A
        1, A --> B
        (d1,d2), (A,B) --> 0
    end
rn8 = @reaction_network name begin
    p, 0 --> A
    1, A --> B
    (d1,d2), (A,B) --> 0
end
@test isequal(parameters(rn7), parameters(rn8))

@parameters p d1 d2
@variables t A(t) B(t)
@test isequal(parameters(rn8)[1],p)
@test isequal(parameters(rn8)[2],d1)
@test isequal(parameters(rn8)[3],d2)
@test isequal(species(rn8)[1],A)
@test isequal(species(rn8)[2],B)