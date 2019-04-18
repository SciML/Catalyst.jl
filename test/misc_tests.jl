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

network4 = @reaction_network begin
    k1, 2A --> B
    k2, B --> 2A
    k3, A + B --> C
    k4, C --> A + B
    k5, 3C --> 3A
    k6, 0 --> 2B
    hill(A,k7,k8,2), ∅ --> B 
end k1 k2 k3 k4 k5 k6 k7 k8
p = rand(8)
u = rand(3)
t = 0.
function n4f(du, u, p, t)
    A = u[1]; B = u[2]; C = u[3];
    k1 = p[1]; k2 = p[2]; k3 = p[3]; k4 = p[4]; k5 = p[5]; k6 = p[6]; k7 = p[7]; k8 = p[8]
    du[1] = -k1*A^2 + 2*k2*B - k3*A*B + k4*C + k5*C^3/2
    du[2] = .5*k1*A^2 -k2*B -k3*A*B + k4*C + 2*k6 + k7*A^2/(k8^2+A^2)
    du[3] = k3*A*B-k4*C-k5*C^3/2
end
function n4Jac(J, u, p, t)
    A = u[1]; B = u[2]; C = u[3];
    k1 = p[1]; k2 = p[2]; k3 = p[3]; k4 = p[4]; k5 = p[5]; k6 = p[6]; k7 = p[7]; k8 = p[8]
    J[1,1] = -2*k1*A - k3*B 
    J[1,2] = 2*k2 - k3*A
    J[1,3] = k4 + 3*k5*C^2 / 2
    J[2,1] = k1*A -k3*B + 2*k7*k8^2*A/(k8^2+A^2)^2
    J[2,2] = -k2-k3*A
    J[2,3] = k4
    J[3,1] = k3*B
    J[3,2] = k3*A
    J[3,3] = -k4 - 3*k5*C^2/2
end
du  = zeros(3); network4.f(du,u,p,t);
dua = zeros(3); n4f(dua,u,p,t);
J   = zeros(3,3); network4.jac(J,u,p,t);
Ja  = zeros(3,3); n4Jac(Ja,u,p,t);
@test maximum(abs.(du .- dua)) < 10*eps(Float64)
@test maximum(abs.(J .- Ja)) < 10*eps(Float64)

# Test the sparse jacobian
network5 = @min_reaction_network begin
    k1, 2A --> B
    k2, B --> 2A
    k3, A + B --> C
    k4, C --> A + B
    k5, 3C --> 3A
    k6, 0 --> 2B
    hill(A,k7,k8,2), ∅ --> B 
end k1 k2 k3 k4 k5 k6 k7 k8
addodes!(network5, sparse_jac=true)
Js = sparse(ones(3,3))
network5.jac(Js, u, p, t)
@test maximum(abs.(Matrix(Js) .- Ja)) < 10*eps(Float64)

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
