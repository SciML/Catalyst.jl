### Fetch required packages and reaction networks ###
using DiffEqBiological, Random, Statistics, StochasticDiffEq, Test
include("test_networks.jl")

### Compares to the manually calcualted function ###
identical_networks = Vector{Pair}()

function real_f_1(du,u,p,t)
    X1,X2,X3 = u
    p,k1,k2,k3,d = p
    du[1] = 2*p - k1*X1
    du[2] = k1*X1 - k2*X2 - k3*X2
    du[3] = k2*X2 + k3*X2 - d*X3
end
function real_g_1(du,u,p,t)
    p,k1,k2,k3,d = p
    X1,X2,X3 = u
    du[1,1] = 2*sqrt(p)
    du[1,2] = -sqrt(k1*X1)
    du[1,3] = 0
    du[1,4] = 0
    du[1,5] = 0
    du[2,1] = 0
    du[2,2] = sqrt(k1*X1)
    du[2,3] = -sqrt(k2*X2)
    du[2,4] = -sqrt(k3*X2)
    du[2,5] = 0
    du[3,1] = 0
    du[3,2] = 0
    du[3,3] = sqrt(k2*X2)
    du[3,4] = sqrt(k3*X2)
    du[3,5] = -sqrt(d*X3)
end
push!(identical_networks, reaction_networks_standard[8] => (real_f_1,real_g_1))

function real_f_2(du,u,p,t)
    X1, = u
    v,K,n,d = p
    du[1] = v/10+hill(X1,v,K,n) - d*X1
end
function real_g_2(du,u,p,t)
    X1, = u
    v,K,n,d = p
    du[1,1] = sqrt(v/10+hill(X1,v,K,n))
    du[1,2] = -sqrt(d*X1)
end
push!(identical_networks, reaction_networks_hill[6] => (real_f_2,real_g_2))

function real_f_3(du,u,p,t)
    X1,X2,X3,X4,X5,X6,X7 = u
    k1,k2,k3,k4,k5,k6 = p
    du[1] = -k1*X1*X2 + k2*X3
    du[2] = -k1*X1*X2 + k2*X3
    du[3] = k1*X1*X2 - k2*X3 - k3*X3*X4 + k4*X5
    du[4] = -k3*X3*X4 + k4*X5
    du[5] = k3*X3*X4 - k4*X5 - k5*X5*X6 + k6*X7
    du[6] = -k5*X5*X6 + k6*X7
    du[6] = k5*X5*X6 - k6*X7
end
function real_g_3(du,u,p,t)
    X1,X2,X3 = u
    k1,k2,k3,k4,k5,k6 = p
    du = zeros(7,6)
    du[1,1] = -sqrt(k1*X1*X2)
    du[1,2] = sqrt(k2*X3)
    du[2,1] = -sqrt(k1*X1*X2)
    du[2,2] = sqrt(k2*X3)
    du[3,1] = sqrt(k1*X1*X2)
    du[3,2] = -sqrt(k2*X3)
    du[3,3] = -sqrt(k3*X3*X4)
    du[3,4] = sqrt(k4*X5)
    du[4,3] = -sqrt(k3*X3*X4)
    du[4,4] = sqrt(k4*X5)
    du[5,3] = sqrt(k3*X3*X4)
    du[5,4] = -sqrt(k4*X5)
    du[5,5] = -sqrt(k5*X5*X6)
    du[5,6] = sqrt(k6*X7)
    du[6,5] = -sqrt(k5*X5*X6)
    du[6,6] = sqrt(k6*X7)
    du[7,5] = sqrt(k5*X5*X6)
    du[7,6] = -sqrt(k6*X7)
end
push!(identical_networks, reaction_networks_constraint[9] => (real_f_3,real_g_3))

for (i,networks) in enumerate(identical_networks)
    for factor in [1e-1, 1e0, 1e1], repeat in 1:5
        u0 = 100. .+ factor*rand(length(networks[1].states))
        p = 0.01 .+ factor*rand(length(networks[1].ps))
        (i==2) && (u0[1] += 1000.)
        (i==3) ? (p[2:2:6] .*= 100.) : (p[1] += 500.)
        prob1 = SDEProblem(networks[1],u0,(0.,100.),p)
        sol1 = solve(prob1,ImplicitEM(),saveat=0.01,maxiters=1e7)
        prob2 = SDEProblem(networks[2][1],networks[2][2],u0,(0.,100.),p)
        sol2 = solve(prob2,ImplicitEM(),saveat=0.01,maxiters=1e7)
        for i = 1:length(u0)
            vals1 = getindex.(sol1.u[1000:end],i);
            vals2 = getindex.(sol1.u[1000:end],i);
            @test 0.8 < mean(vals1)/mean(vals2) < 1.25
            @test 0.8 < std(vals1)/std(vals2) < 1.25
        end
    end
end


### Tries to create a large number of problem, ensuring there are no errors (cannot solve as solution likely to go into negatives). ###
for reaction_network in reaction_networks_all
    for factor in [1e-2, 1e-1, 1e0, 1e1]
        u0 = factor*rand(length(reaction_network.states))
        p = factor*rand(length(reaction_network.ps))
        prob = SDEProblem(reaction_network,u0,(0.,1.),p)
    end
end
