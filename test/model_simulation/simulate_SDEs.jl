### Fetch Packages and Reaction Networks ###

# Fetch packages.
using Catalyst, Random, Statistics, StochasticDiffEq, Test
using ModelingToolkit: get_postprocess_fbody

# Sets rnd number.
using StableRNGs
rng = StableRNG(12345)

# Fetch test networks.
include("../test_networks.jl")

### Compares to Manually Calcualted Function ###

let
    identical_networks = Vector{Pair}()

    function real_f_1(du, u, p, t)
        X1, X2, X3 = u
        p, k1, k2, k3, d = p
        du[1] = 2 * p - k1 * X1
        du[2] = k1 * X1 - k2 * X2 - k3 * X2
        du[3] = k2 * X2 + k3 * X2 - d * X3
    end
    function real_g_1(du, u, p, t)
        X1, X2, X3 = u
        p, k1, k2, k3, d = p
        du[1, 1] = 2 * sqrt(p)
        du[1, 2] = -sqrt(k1 * X1)
        du[1, 3] = 0
        du[1, 4] = 0
        du[1, 5] = 0
        du[2, 1] = 0
        du[2, 2] = sqrt(k1 * X1)
        du[2, 3] = -sqrt(k2 * X2)
        du[2, 4] = -sqrt(k3 * X2)
        du[2, 5] = 0
        du[3, 1] = 0
        du[3, 2] = 0
        du[3, 3] = sqrt(k2 * X2)
        du[3, 4] = sqrt(k3 * X2)
        du[3, 5] = -sqrt(d * X3)
    end
    push!(identical_networks,
          reaction_networks_standard[8] => (real_f_1, real_g_1, zeros(3, 5)))

    function real_f_2(du, u, p, t)
        X1, = u
        v, K, n, d = p
        du[1] = v / 10 + v * X1^n / (X1^n + K^n) - d * X1
    end
    function real_g_2(du, u, p, t)
        X1, = u
        v, K, n, d = p
        du[1, 1] = sqrt(v / 10 + v * X1^n / (X1^n + K^n))
        du[1, 2] = -sqrt(d * X1)
    end
    push!(identical_networks,
          reaction_networks_hill[6] => (real_f_2, real_g_2, zeros(1, 2)))

    function real_f_3(du, u, p, t)
        X1, X2, X3, X4, X5, X6, X7 = u
        k1, k2, k3, k4, k5, k6 = p
        du[1] = -k1 * X1 * X2 + k2 * X3
        du[2] = -k1 * X1 * X2 + k2 * X3
        du[3] = k1 * X1 * X2 - k2 * X3 - k3 * X3 * X4 + k4 * X5
        du[4] = -k3 * X3 * X4 + k4 * X5
        du[5] = k3 * X3 * X4 - k4 * X5 - k5 * X5 * X6 + k6 * X7
        du[6] = -k5 * X5 * X6 + k6 * X7
        du[7] = k5 * X5 * X6 - k6 * X7
    end
    function real_g_3(du, u, p, t)
        X1, X2, X3, X4, X5, X6, X7 = u
        k1, k2, k3, k4, k5, k6 = p
        fill!(du, 0)
        du[1, 1] = -sqrt(k1 * X1 * X2)
        du[1, 2] = sqrt(k2 * X3)
        du[2, 1] = -sqrt(k1 * X1 * X2)
        du[2, 2] = sqrt(k2 * X3)
        du[3, 1] = sqrt(k1 * X1 * X2)
        du[3, 2] = -sqrt(k2 * X3)
        du[3, 3] = -sqrt(k3 * X3 * X4)
        du[3, 4] = sqrt(k4 * X5)
        du[4, 3] = -sqrt(k3 * X3 * X4)
        du[4, 4] = sqrt(k4 * X5)
        du[5, 3] = sqrt(k3 * X3 * X4)
        du[5, 4] = -sqrt(k4 * X5)
        du[5, 5] = -sqrt(k5 * X5 * X6)
        du[5, 6] = sqrt(k6 * X7)
        du[6, 5] = -sqrt(k5 * X5 * X6)
        du[6, 6] = sqrt(k6 * X7)
        du[7, 5] = sqrt(k5 * X5 * X6)
        du[7, 6] = -sqrt(k6 * X7)
    end
    push!(identical_networks,
          reaction_networks_constraint[9] => (real_f_3, real_g_3, zeros(7, 6)))

    for (i, networks) in enumerate(identical_networks)
        for factor in [1e-1, 1e0, 1e1], repeat in 1:3
            u0 = 100.0 .+ factor * rand(rng, length(unknowns(networks[1])))
            p = 0.01 .+ factor * rand(rng, length(parameters(networks[1])))
            (i == 2) && (u0[1] += 1000.0)
            (i == 3) ? (p[2:2:6] .*= 1000.0; u0 .+= 1000) : (p[1] += 500.0)
            prob1 = SDEProblem(networks[1], u0, (0.0, 100.0), p)
            prob2 = SDEProblem(networks[2][1], networks[2][2], u0, (0.0, 100.0), p,
                               noise_rate_prototype = networks[2][3])
            du1 = similar(u0)
            du2 = similar(u0)
            prob1.f.f(du1, u0, p, 0.0)
            prob2.f.f(du2, u0, p, 0.0)
            @test all(isapprox.(du1, du2))
            g1 = zeros(numspecies(networks[1]), numreactions(networks[1]))
            g2 = copy(g1)
            prob1.f.g(g1, u0, p, 0.0)
            prob2.f.g(g2, u0, p, 0.0)
            @test all(isapprox.(g1, g2))
        end
    end
end

### Checks Simulations Don't Error ###

# Tries to create a large number of problem, ensuring there are no errors (cannot solve as solution likely to go into negatives). 
let
    for reaction_network in reaction_networks_all
        for factor in [1e-2, 1e-1, 1e0, 1e1]
            u0 = factor * rand(rng, length(states(reaction_network)))
            p = factor * rand(rng, length(parameters(reaction_network)))
            prob = SDEProblem(reaction_network, u0, (0.0, 1.0), p)
        end
    end
end

### Noise Scaling ###

# Tests with multiple noise scaling parameters directly in the macro.
let 
    noise_scaling_network_1 = @reaction_network begin 
        @parameters η1 η2
        (k1, k2), X1 ↔ X2, ([noise_scaling=η1],[noise_scaling=η2]) 
    end
    u0 = [:X1 => 1000.0, :X2 => 3000.0]
    sol_1_1 = solve(SDEProblem(noise_scaling_network_1, u0, (0.0, 1000.0), [:k1 => 2.0, :k2 => 0.66, :η1 => 2.0, :η2 => 2.0]), ImplicitEM(); saveat=1.0)
    sol_1_2 = solve(SDEProblem(noise_scaling_network_1, u0, (0.0, 1000.0), [:k1 => 2.0, :k2 => 0.66, :η1 => 2.0, :η2 => 0.2]), ImplicitEM(); saveat=1.0)
    sol_1_3 = solve(SDEProblem(noise_scaling_network_1, u0, (0.0, 1000.0), [:k1 => 2.0, :k2 => 0.66, :η1 => 0.2, :η2 => 0.2]), ImplicitEM(); saveat=1.0)
    @test var(sol_1_1[1,:]) > var(sol_1_2[1,:]) > var(sol_1_3[1,:]) 

    noise_scaling_network_2 = @reaction_network begin 
        @parameters η[1:2]
        (k1, k2), X1 ↔ X2, ([noise_scaling=η[1]],[noise_scaling=η[2]])  
    end
    @unpack k1, k2, η = noise_scaling_network_2
    sol_2_1 = solve(SDEProblem(noise_scaling_network_2, u0, (0.0, 1000.0), [k1 => 2.0, k2 => 0.66, η[1] => 2.0, η[2] => 2.0]), ImplicitEM(); saveat=1.0)
    sol_2_2 = solve(SDEProblem(noise_scaling_network_2, u0, (0.0, 1000.0), [k1 => 2.0, k2 => 0.66, η[1] => 2.0, η[2] => 0.2]), ImplicitEM(); saveat=1.0)
    sol_2_3 = solve(SDEProblem(noise_scaling_network_2, u0, (0.0, 1000.0), [k1 => 2.0, k2 => 0.66, η[1] => 0.2, η[2] => 0.2]), ImplicitEM(); saveat=1.0)
    @test var(sol_2_1[1,:]) > var(sol_2_2[1,:]) > var(sol_2_3[1,:]) 
end

# Tests using default values for noise scaling.
# Tests when reaction system is created programmatically.
# Tests @noise_scaling_parameters macro.
let
    η_stored = :η
    @variables t
    @species X1(t) X2(t)
    p_syms = @parameters $(η_stored) k1 k2

    r1 = Reaction(k1,[X1],[X2],[1],[1]; metadata = [:noise_scaling => η_stored])
    r2 = Reaction(k2,[X2],[X1],[1],[1]; metadata = [:noise_scaling => η_stored])
    @named noise_scaling_network = ReactionSystem([r1, r2], t, [X1, X2], [k1, k2, p_syms[1]])

    u0 = [:X1 => 1100.0, :X2 => 3900.0]
    p = [:k1 => 2.0, :k2 => 0.5, :η=>0.0]
    @test_broken SDEProblem(noise_scaling_network, u0, (0.0, 1000.0), p).ps[:η] == 0.0 # Broken due to SII/MTK stuff.
end

# Complicated test with many combinations of options.
# Tests the noise_scaling_parameters getter.
let
    noise_scaling_network = @reaction_network begin 
        @parameters k1 par1 [description="Parameter par1"] par2 η1 η2=0.0 [description="Parameter η2"] η3=1.0 η4
        (p, d), 0 ↔ X1, ([noise_scaling=η1],[noise_scaling=η2]) 
        (k1, k2), X1 ↔ X2, ([noise_scaling=η3],[noise_scaling=η4])  
    end
    u0 = [:X1 => 500.0, :X2 => 500.0]
    p = [:p => 20.0, :d => 0.1, :η1 => 0.0, :η3 => 0.0, :η4 => 0.0, :k1 => 2.0, :k2 => 2.0, :par1 => 1000.0, :par2 => 1000.0]
    
    @test getdescription(parameters(noise_scaling_network)[2]) == "Parameter par1"
    @test getdescription(parameters(noise_scaling_network)[5]) == "Parameter η2"
    
    sprob = SDEProblem(noise_scaling_network, u0, (0.0, 1000.0), p)
    @test sprob.ps[:η1] == sprob.ps[:η2] == sprob.ps[:η3] == sprob.ps[:η4] == 0.0
end

# Tests default_noise_scaling_option.
# Tests using sys. indexing for setting noise scaling parameter values.
# tests for default value used in noise scaling parameter.
let
    noise_scaling_network = @reaction_network begin
        @parameters η1 η2=0.1
        @default_noise_scaling η1
        (p,d), 0 <--> X1, ([noise_scaling=η2],[noise_scaling=η2])
        (p,d), 0 <--> X2, ([description="Y"],[description="Y"])
        (p,d), 0 <--> X3
    end
    noise_scaling_network = complete(noise_scaling_network)
    
    u0 = [:X1 => 1000.0, :X2 => 1000.0, :X3 => 1000.0]
    ps = [noise_scaling_network.p => 1000.0, noise_scaling_network.d => 1.0, noise_scaling_network.η1 => 1.0]
    sol = solve(SDEProblem(noise_scaling_network, u0, (0.0, 1000.0), ps), ImplicitEM(); saveat=1.0)
    @test var(sol[:X1]) < var(sol[:X2]) 
    @test var(sol[:X1]) < var(sol[:X3]) 
end

# Tests  using complicated noise scaling expressions
let
    noise_scaling_network = @reaction_network begin
        @parameters η1 η2 η3 η4
        @species N1(t) N2(t)=0.5
        @variables N3(t)
        @default_noise_scaling η1 + N1 + 5.0
        p, 0 --> X1 
        p, 0 --> X2, [noise_scaling=0.33η2^1.2 + N2]
        p, 0 --> X3, [noise_scaling=N3*η3] 
        p, 0 --> X4, [noise_scaling=exp(-η4) - 0.008]
        p, 0 --> X5, [noise_scaling=0.0]
        d, (X1, X2, X3, X4, X5) --> 0, ([], [noise_scaling=0.33η2^1.2 + N2], [noise_scaling=N3*η3], [noise_scaling=exp(-η4) - 0.008], [noise_scaling=0.0])
    end

    u0 = [:X1 => 1000.0, :X2 => 1000.0, :X3 => 1000.0, :X4 => 1000.0, :X5 => 1000.0, :N1 => 3.0, :N3 => 0.33]
    ps = [:p => 1000.0, :d => 1.0, :η1 => 1.0, :η2 => 1.4, :η3 => 0.33, :η4 => 4.0]
    sol = solve(SDEProblem(noise_scaling_network, u0, (0.0, 1000.0), ps), ImplicitEM(); saveat=1.0, adaptive=false, dt=0.1)
    @test var(sol[:X1]) > var(sol[:X2]) > var(sol[:X3]) > var(sol[:X4]) > var(sol[:X5])  
end


### Checks Simulations Don't Error ###

#Tries to create a large number of problem, ensuring there are no errors (cannot solve as solution likely to go into negatives). 
let
    for reaction_network in reaction_networks_all
        for factor in [1e-2, 1e-1, 1e0, 1e1]
            u0 = factor * rand(rng, length(unknowns(reaction_network)))
            p = factor * rand(rng, length(parameters(reaction_network)))
            prob = SDEProblem(reaction_network, u0, (0.0, 1.0), p)
        end
    end
end

### Other Tests ###

# No parameter test.
let
    no_param_network = @reaction_network begin (1.2, 5), X1 ↔ X2 end
    for factor in [1e3, 1e4]
        u0 = factor * (1.0 .+ rand(rng, length(unknowns(no_param_network))))
        prob = SDEProblem(no_param_network, u0, (0.0, 1000.0))
        sol = solve(prob, ImplicitEM())
        vals1 = getindex.(sol.u[1:end], 1)
        vals2 = getindex.(sol.u[1:end], 2)
        @test mean(vals1) > mean(vals2)
    end
end
