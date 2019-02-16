#Tests the jacobian
network1 = @reaction_network begin
    (2.0,1.0),   ∅ ↔ X
    (3.0,1.0),   ∅ ↔ Y
    (5.0,2.0),   X + Y ↔ XY
end
J = zeros(3,3)
network1.jac(J,[1.,1.,1.],[],0.)
@test J == [-6. -5. 2.; -5. -6. 2.; 5. 5. -2.]

network2 = @min_reaction_network begin
    (2.0,1.0),   ∅ ↔ X
    (3.0,1.0),   ∅ ↔ Y
    (5.0,2.0),   X + Y ↔ XY
end
addodes!(network2)
J = zeros(3,3)
network2.jac(J,[1.,1.,1.],[],0.)
@test J == [-6. -5. 2.; -5. -6. 2.; 5. 5. -2.]

network3 = @reaction_network begin
    (p1,1.0),   ∅ ↔ X
    (p2,1.0),   ∅ ↔ Y
    (p3*X,1.0), X + Y ↔ XY
end p1 p2 p3
J = zeros(3,3)
for i = 1:100
    u = 10*rand(3)
    p = 10*rand(3)
    network3.jac(J,u,p,0.)
    @test maximum(J .- [-1-2*p[3]*u[2]*u[1] -p[3]*u[1]*u[1] 1.; -2*p[3]*u[2]*u[1] -1-p[3]*u[1]*u[1] 1; 2*p[3]*u[2]*u[1] p[3]*u[1]*u[1] -1.])<0.0001
end

#Tests the parameter jacobian.
network1 = @reaction_network begin
    (p1,1.0),   ∅ ↔ X
    (1.0,p2),   X ↔ Y
    (p3,p4),    Y ↔ Z
end p1 p2 p3 p4
J = zeros(3,4)
network1.paramjac(J,[2.,3.,4.],[2.,3.,4.,5.],0.)
@test J == [1. 3. 0. 0.; -0. -3. -3. 4.; 0. 0. 3. -4.]

network2 = @reaction_network begin
    (1.,p1*p2),   ∅ ↔ X
    (p3^3+p1,p2),   2X ↔ Y
end p1 p2 p3
J = zeros(2,3)
for i = 1:100
    u = 10*rand(2)
    p = 10*rand(3)
    network2.paramjac(J,u,p,0.)
    @test maximum(J .- [-u[1]*p[2]-u[1]^2 -u[1]*p[1]+2*u[2] -6*p[3]^2*u[1]^2/2; u[1]^2/2 -u[2] 3*p[3]^2*u[1]^2/2])<0.0001
end
