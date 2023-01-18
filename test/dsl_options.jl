#! format: off

using Catalyst, ModelingToolkit, OrdinaryDiffEq

### Test creating networks with/without options. ###
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

rn5 = @reaction_network name begin (k1, k2), A <--> B end
rn6 = @reaction_network name begin
    @parameters k2 k1
    @species B A
    (k1, k2), A <--> B
end
#@test !isequal(species(rn5),species(rn6)) # This does not work, as ReactionSystem does not accept ordering of species, but does that itself.
@test !isequal(parameters(rn5), parameters(rn6))


### Checks that the rights things are put in vectors. ###
rn7 = @reaction_network name begin
    @parameters p d1 d2
    @species A B
    p, 0 --> A
    1, A --> B
    (d1, d2), (A, B) --> 0
end
rn8 = @reaction_network name begin
    p, 0 --> A
    1, A --> B
    (d1, d2), (A, B) --> 0
end
@test isequal(parameters(rn7), parameters(rn8))

@parameters p d1 d2
@variables t A(t) B(t)
@test isequal(parameters(rn8)[1], p)
@test isequal(parameters(rn8)[2], d1)
@test isequal(parameters(rn8)[3], d2)
@test isequal(species(rn8)[1], A)
@test isequal(species(rn8)[2], B)


### Tests that order is preserved when set ###
rn9 = @reaction_network name begin
    @species X1 X2 X3 X4
    k4, 0 --> X4
    k3, 0 --> X3
    k2, 0 --> X2
    k1, 0 --> X1
end
@test isequal(Symbol.(getfield.(species(rn9),:f)),[:X1, :X2, :X3, :X4]) # Should use other notation to get syms, but can't find it.

rn10 = @reaction_network name begin
    @parameters k1 k2 k3 k4
    k4, 0 --> X4
    k3, 0 --> X3
    k2, 0 --> X2
    k1, 0 --> X1
end
@test isequal(Symbol.(parameters(rn10)),[:k1, :k2, :k3, :k4])

rn11 = @reaction_network name begin
    @species X1 X2 X3 X4 Y1 Y2 Y3 Y4
    @parameters k1 k2 k3 k4 l1 l2 l3 l4
    k4*Y3+l4, 0 --> X4 + Y4
    k3*Y1, 0 --> X3
    k2+l2+l1, Y2 --> X2
    k1, 0 --> X1
end
@test isequal(Symbol.(getfield.(species(rn11),:f)),[:X1, :X2, :X3, :X4, :Y1, :Y2, :Y3, :Y4])
@test isequal(Symbol.(parameters(rn11)),[:k1, :k2, :k3, :k4, :l1, :l2, :l3, :l4])


#### Tests that defaults work. ###
rn12 = @reaction_network name begin
    @parameters p=1.0 d1 d2=5
    @species A B=4
    p, 0 --> A
    1, A --> B
    (d1, d2), (A, B) --> 0
end

rn13 = @reaction_network name begin
@parameters p1=1.0 p2=2.0 k1=4.0 k2=5.0 v=8.0 K=9.0 n=3 d=10.0
@species X=4.0 Y=3.0 X2Y=2.0 Z=1.0
    (p1,p2), 0 --> (X,Y)
    (k1,k2), 2X + Y --> X2Y
    hill(X2Y,v,K,n), 0 --> Z
    d, (X,Y,X2Y,Z) --> 0
end
u0_13 = []
p_13 = []

rn14 = @reaction_network name begin
@parameters p1=1.0 p2 k1=4.0 k2 v=8.0 K n=3 d
@species X=4.0 Y X2Y Z=1.0
    (p1,p2), 0 --> (X,Y)
    (k1,k2), 2X + Y --> X2Y
    hill(X2Y,v,K,n), 0 --> Z
    d, (X,Y,X2Y,Z) --> 0
end
u0_14 = [:p2=>2.0, :k2=>5.0, :K=>9.0, :d=>10.0]
p_14 = [:Y=>3.0, :X2Y=>2.0]

rn15 = @reaction_network name begin
@parameters p1 p2 k1 k2 v K n d
@species X Y X2Y Z
    (p1,p2), 0 --> (X,Y)
    (k1,k2), 2X + Y --> X2Y
    hill(X2Y,v,K,n), 0 --> Z
    d, (X,Y,X2Y,Z) --> 0
end
u0_15 = [:p1=>1.0, :p2=>2.0, :k1=>4.0, :k2=>5.0, :v=>8.0, :K=>9.0, :n=>3, :d=>10.0]
p_15 = [:X=>4.0, :Y=>3.0, :X2Y=>2.0, :Z=>1.0]

uEnd_13 = solve(ODEProblem(rn10,u0_13,(0.0,10.0),p_13),Rosenbrock23()).u[end]
uEnd_14 = solve(ODEProblem(rn11,u0_14,(0.0,10.0),p_14),Rosenbrock23()).u[end]
uEnd_15 = solve(ODEProblem(rn12,u0_15,(0.0,10.0),p_15),Rosenbrock23()).u[end]

@test isapprox(uEnd_13, uEnd_14, rtol = 1e3 * eps())
@test isapprox(uEnd_14, uEnd_15, rtol = 1e3 * eps())
