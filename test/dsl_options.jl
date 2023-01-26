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
    @species A(t) B(t)
    (k1, k2), A <--> B
end
@reaction_network begin
    @species A(t) B(t)
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
        A(t)
        B(t)
    end
    (k1, k2), A <--> B
end
@reaction_network begin
    @parameters begin
        k1
        k2
    end
    @species begin
        A(t)
        B(t)
    end
    (k1, k2), A <--> B
end

@reaction_network name begin (k1, k2), A <--> B end
@reaction_network name begin
    @parameters k1 k2
    (k1, k2), A <--> B
end
@reaction_network name begin
    @species A(t) B(t)
    (k1, k2), A <--> B
end
@reaction_network name begin
    @parameters k1 k2
    @species A(t) B(t)
    (k1, k2), A <--> B
end

@reaction_network name begin (k1, k2), A <--> B end
@reaction_network name begin
    (k1, k2), A <--> B
    @parameters k1 k2
end
@reaction_network name begin
    (k1, k2), A <--> B
    @species A(t) B(t)
end
@reaction_network name begin
    (k1, k2), A <--> B
    @parameters k1 k2
    @species A(t) B(t)
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
        A(t)
        B(t)
    end
    (k1, k2), A <--> B
end
@reaction_network name begin
    @parameters begin
        k1
        k2
    end
    @species begin
        A(t)
        B(t)
    end
    (k1, k2), A <--> B
end


### Tests that when either @species or @parameters is given, the other is infered properly. ###
@variables t

rn1 = @reaction_network begin
    k*X, A + B --> 0
end
@test issetequal(species(rn1), @species A(t) B(t))
@test issetequal(parameters(rn1), @parameters k X) 

rn2 = @reaction_network begin
    @species A(t) B(t) X(t)
    k*X, A + B --> 0
end
@test issetequal(species(rn2), @species A(t) B(t) X(t))
@test issetequal(parameters(rn2), @parameters k) 

rn3 = @reaction_network begin
    @parameters k
    k*X, A + B --> 0
end
@test issetequal(species(rn3), @species A(t) B(t))
@test issetequal(parameters(rn3), @parameters k X) 

rn4 = @reaction_network begin
    @species A(t) B(t) X(t)
    @parameters k
    k*X, A + B --> 0
end
@test issetequal(species(rn4), @species A(t) B(t) X(t))
@test issetequal(parameters(rn4), @parameters k) 

rn5 = @reaction_network begin
    @parameters k B [isconstantspecies=true]
    k*X, A + B --> 0
end
@test issetequal(species(rn5), @species A(t))
@test issetequal(parameters(rn5), @parameters k B X) 


### Tests that when some species or parameters are left out, the others are set properly. ###
@variables t

rn6 = @reaction_network begin
    @species A(t) 
    k*X, A + B --> 0
end
@test issetequal(species(rn6), @species A(t) B(t))
@test issetequal(parameters(rn6), @parameters k X) 

rn7 = @reaction_network begin
    @species A(t) X(t) 
    k*X, A + B --> 0
end
@test issetequal(species(rn7), @species A(t) X(t) B(t))
@test issetequal(parameters(rn7), @parameters k) 

rn7 = @reaction_network begin
    @parameters B [isconstantspecies=true]
    k*X, A + B --> 0
end
@test issetequal(species(rn7), @species A(t))
@test issetequal(parameters(rn7), @parameters B k X) 

rn8 = @reaction_network begin
    @parameters B [isconstantspecies=true] k 
    k*X, A + B --> 0
end
@test issetequal(species(rn8), @species A(t))
@test issetequal(parameters(rn8), @parameters B k X) 

rn9 = @reaction_network begin
    @parameters k1 X1 
    @species A1(t) B1(t)
    k1*X1, A1 + B1 --> 0
    k2*X2, A2 + B2 --> 0
end
@test issetequal(species(rn9), @species A1(t) B1(t) A2(t) B2(t))
@test issetequal(parameters(rn9), @parameters k1 X1 k2 X2) 

rn10 = @reaction_network begin
    @parameters k1 X2 B2 [isconstantspecies=true] 
    @species A1(t) X1(t)
    k1*X1, A1 + B1 --> 0
    k2*X2, A2 + B2 --> 0
end
@test issetequal(species(rn10), @species A1(t) X1(t) B1(t) A2(t))
@test issetequal(parameters(rn10), @parameters k1 X2 B2 k2) 

rn11 = @reaction_network begin
    @parameters k1 k2
    @species X1(t)
    k1*X1, A1 + B1 --> 0
    k2*X2, A2 + B2 --> 0
end
@test issetequal(species(rn11), @species X1(t) A1(t) A2(t) B1(t) B2(t))
@test issetequal(parameters(rn11), @parameters k1 k2 X2) 


### Checks that some created networks are identical. ###
rn12 = @reaction_network name begin (k1, k2), A <--> B end
rn13 = @reaction_network name begin
    @parameters k1 k2
    (k1, k2), A <--> B
end
rn14 = @reaction_network name begin
    @species A(t) B(t)
    (k1, k2), A <--> B
end
rn15 = @reaction_network name begin
    @parameters k1 k2
    @species A(t) B(t)
    (k1, k2), A <--> B
end
@test isequal(species(rn12),species(rn13))
@test isequal(species(rn13),species(rn14))
@test isequal(species(rn14),species(rn15))
@test isequal(parameters(rn12), parameters(rn13))
@test isequal(parameters(rn13), parameters(rn14))
@test isequal(parameters(rn14), parameters(rn15))

rn16 = @reaction_network name begin (k1, k2), A <--> B end
rn17 = @reaction_network name begin
    @parameters k2 k1
    @species B(t) A(t)
    (k1, k2), A <--> B
end
@test !isequal(species(rn16),species(rn17))
@test !isequal(parameters(rn16), parameters(rn17))


### Checks that the rights things are put in vectors. ###
rn18 = @reaction_network name begin
    @parameters p d1 d2
    @species A(t) B(t)
    p, 0 --> A
    1, A --> B
    (d1, d2), (A, B) --> 0
end
rn19 = @reaction_network name begin
    p, 0 --> A
    1, A --> B
    (d1, d2), (A, B) --> 0
end
@test isequal(parameters(rn18), parameters(rn19))

@parameters p d1 d2
@variables t A(t) B(t)
@test isequal(parameters(rn18)[1], p)
@test isequal(parameters(rn18)[2], d1)
@test isequal(parameters(rn18)[3], d2)
@test isequal(species(rn18)[1], A)
@test isequal(species(rn18)[2], B)

rn20 = @reaction_network name begin
    @species X(t)
    @parameters S
    mm(X,v,K), 0 --> Y
    (k1,k2), 2Y <--> Y2
    d*Y, S*(Y2+Y) --> 0
end
rn21 = @reaction_network name begin
    @species X(t) Y(t) Y2(t)
    @parameters v K k1 k2 d S
    mm(X,v,K), 0 --> Y
    (k1,k2), 2Y <--> Y2
    d*Y, S*(Y2+Y) --> 0
end
rn22 = @reaction_network name begin
    @species X(t) Y2(t)
    @parameters d k1
    mm(X,v,K), 0 --> Y
    (k1,k2), 2Y <--> Y2
    d*Y, S*(Y2+Y) --> 0
end

@parameters v K k1 k2 d S
@variables t X(t) Y(t) Y2(t)
@test length(parameters(rn20))==length(parameters(rn21))==length(parameters(rn22))==6
@test length(species(rn20))==length(species(rn21))==length(species(rn22))==3
@test issetequal(parameters(rn20),parameters(rn21))
@test issetequal(parameters(rn21),parameters(rn22))
@test issetequal(parameters(rn22),[v K k1 k2 d S])
@test issetequal(species(rn20), species(rn21))
@test issetequal(species(rn21), species(rn22))
@test issetequal(species(rn22), [X Y Y2])

### Tests that order is preserved when set. ###
rn23 = @reaction_network name begin
    @species X1(t) X2(t) X3(t) X4(t)
    k4, 0 --> X4
    k3, 0 --> X3
    k2, 0 --> X2
    k1, 0 --> X1
end
@test isequal(map(Symbol ∘ ModelingToolkit.operation, states(rn23)),[:X1, :X2, :X3, :X4])

rn24 = @reaction_network name begin
    @parameters k1 k2 k3 k4
    k4, 0 --> X4
    k3, 0 --> X3
    k2, 0 --> X2
    k1, 0 --> X1
end
@test isequal(map(Symbol, parameters(rn24)),[:k1, :k2, :k3, :k4])

rn25 = @reaction_network name begin
    @species X1(t) X2(t) X3(t) X4(t) Y1(t) Y2(t) Y3(t) Y4(t)
    @parameters k1 k2 k3 k4 l1 l2 l3 l4
    k4*Y3+l4, 0 --> X4 + Y4
    k3*Y1, 0 --> X3
    k2+l2+l1, Y2 --> X2
    k1, 0 --> X1
end
@test isequal(map(Symbol ∘ ModelingToolkit.operation, states(rn25)),[:X1, :X2, :X3, :X4, :Y1, :Y2, :Y3, :Y4])
@test isequal(map(Symbol, parameters(rn25)),[:k1, :k2, :k3, :k4, :l1, :l2, :l3, :l4])


#### Tests that defaults work. ###
rn26 = @reaction_network name begin
    @parameters p=1.0 d1 d2=5
    @species A(t) B(t)=4
    p, 0 --> A
    1, A --> B
    (d1, d2), (A, B) --> 0
end

rn27 = @reaction_network name begin
@parameters p1=1.0 p2=2.0 k1=4.0 k2=5.0 v=8.0 K=9.0 n=3 d=10.0
@species X(t)=4.0 Y(t)=3.0 X2Y(t)=2.0 Z(t)=1.0
    (p1,p2), 0 --> (X,Y)
    (k1,k2), 2X + Y --> X2Y
    hill(X2Y,v,K,n), 0 --> Z
    d, (X,Y,X2Y,Z) --> 0
end
u0_27 = []
p_27 = []

rn28 = @reaction_network name begin
@parameters p1=1.0 p2 k1=4.0 k2 v=8.0 K n=3 d
@species X(t)=4.0 Y(t) X2Y(t) Z(t)=1.0
    (p1,p2), 0 --> (X,Y)
    (k1,k2), 2X + Y --> X2Y
    hill(X2Y,v,K,n), 0 --> Z
    d, (X,Y,X2Y,Z) --> 0
end
u0_28 = [:p2=>2.0, :k2=>5.0, :K=>9.0, :d=>10.0]
p_28 = [:Y=>3.0, :X2Y=>2.0]

rn29 = @reaction_network name begin
@parameters p1 p2 k1 k2 v K n d
@species X(t) Y(t) X2Y(t) Z(t)
    (p1,p2), 0 --> (X,Y)
    (k1,k2), 2X + Y --> X2Y
    hill(X2Y,v,K,n), 0 --> Z
    d, (X,Y,X2Y,Z) --> 0
end
u0_29 = [:p1=>1.0, :p2=>2.0, :k1=>4.0, :k2=>5.0, :v=>8.0, :K=>9.0, :n=>3, :d=>10.0]
p_29 = [:X=>4.0, :Y=>3.0, :X2Y=>2.0, :Z=>1.0]

uEnd_27 = solve(ODEProblem(rn27,u0_27,(0.0,10.0),p_27),Rosenbrock23()).u[end]
uEnd_28 = solve(ODEProblem(rn28,u0_28,(0.0,10.0),p_28),Rosenbrock23()).u[end]
uEnd_29 = solve(ODEProblem(rn29,u0_29,(0.0,10.0),p_29),Rosenbrock23()).u[end]

@test isapprox(uEnd_27, uEnd_28, rtol = 1e3 * eps())
@test isapprox(uEnd_28, uEnd_29, rtol = 1e3 * eps())
