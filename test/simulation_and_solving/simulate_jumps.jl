### Prepares Tests ###

# Fetch packages.
using Catalyst, JumpProcesses, Statistics, Test

# Sets stable rng number.
using StableRNGs
rng = StableRNG(12345)
seed = rand(rng, 1:100)

# Fetch test functions and networks.
include("../test_functions.jl")
include("../test_networks.jl")


### Basic Tests ###

# Compares jump simulations generated through hCatalyst, and manually created systems.
let
    # Declares vectors to store information about each network in. Initial conditions are approximately
    # at the steady state.
    catalyst_networks = []
    manual_networks = []
    u0_syms = []
    ps_syms = []
    u0s = []
    ps = []
    sps = []

    # Manual declaration of the first network (and its simulation settings).
    rate_1_1(u, p, t) = p[1]
    rate_1_2(u, p, t) = p[2] * u[1]
    rate_1_3(u, p, t) = p[3] * u[2]
    rate_1_4(u, p, t) = p[4] * u[2]
    rate_1_5(u, p, t) = p[5] * u[3]
    rate_1_6(u, p, t) = p[6] * u[3]
    rate_1_7(u, p, t) = p[7] * u[4]
    rate_1_8(u, p, t) = p[8] * u[4]
    affect_1_1!(integrator) = (integrator.u[1] += 1)
    affect_1_2!(integrator) = (integrator.u[1] -= 1; integrator.u[2] += 1)
    affect_1_3!(integrator) = (integrator.u[2] -= 1; integrator.u[1] += 1)
    affect_1_4!(integrator) = (integrator.u[2] -= 1; integrator.u[3] += 1)
    affect_1_5!(integrator) = (integrator.u[3] -= 1; integrator.u[2] += 1)
    affect_1_6!(integrator) = (integrator.u[3] -= 1; integrator.u[4] += 1)
    affect_1_7!(integrator) = (integrator.u[4] -= 1; integrator.u[3] += 1)
    affect_1_8!(integrator) = (integrator.u[4] -= 1)
    jump_1_1 = ConstantRateJump(rate_1_1, affect_1_1!)
    jump_1_2 = ConstantRateJump(rate_1_2, affect_1_2!)
    jump_1_3 = ConstantRateJump(rate_1_3, affect_1_3!)
    jump_1_4 = ConstantRateJump(rate_1_4, affect_1_4!)
    jump_1_5 = ConstantRateJump(rate_1_5, affect_1_5!)
    jump_1_6 = ConstantRateJump(rate_1_6, affect_1_6!)
    jump_1_7 = ConstantRateJump(rate_1_7, affect_1_7!)
    jump_1_8 = ConstantRateJump(rate_1_8, affect_1_8!)
    jumps_1 = (jump_1_1, jump_1_2, jump_1_3, jump_1_4, jump_1_5, jump_1_6, jump_1_7,
                jump_1_8)
    push!(catalyst_networks, reaction_networks_standard[5])
    push!(manual_networks, jumps_1)
    push!(u0_syms, [:X1, :X2, :X3, :X4])
    push!(ps_syms, [:p, :k1, :k2, :k3, :k4, :k5, :k6, :d])
    push!(u0s, [:X1 => 40, :X2 => 30, :X3 => 20, :X4 => 10])
    push!(ps, [:p => 10.0, :k1 => 1.0, :k2 => 1.0, :k3 => 1.0, :k4 => 1.0, :k5 => 1.0, :k6 => 1.0, :d => 1.0])
    push!(sps, :X4)

    # Manual declaration of the second network (and its simulation settings).
    rate_2_1(u, p, t) = p[1] / 10 + p[1] * (u[1]^p[3]) / (u[1]^p[3] + p[2]^p[3])
    rate_2_2(u, p, t) = p[4] * u[1] * u[2]
    rate_2_3(u, p, t) = p[5] * u[3]
    rate_2_4(u, p, t) = p[6] * u[3]
    rate_2_5(u, p, t) = p[7] * u[1]
    rate_2_6(u, p, t) = p[7] * u[2]
    rate_2_7(u, p, t) = p[7] * u[3]
    affect_2_1!(integrator) = (integrator.u[1] += 1; integrator.u[2] += 1)
    function affect_2_2!(integrator)
        (integrator.u[1] -= 1; integrator.u[2] -= 1; integrator.u[3] += 1)
    end
    function affect_2_3!(integrator)
        (integrator.u[1] += 1; integrator.u[2] += 1; integrator.u[3] -= 1)
    end
    affect_2_4!(integrator) = (integrator.u[3] -= 1; integrator.u[1] += 1)
    affect_2_5!(integrator) = (integrator.u[1] -= 1)
    affect_2_6!(integrator) = (integrator.u[2] -= 1)
    affect_2_7!(integrator) = (integrator.u[3] -= 1)
    jump_2_1 = ConstantRateJump(rate_2_1, affect_2_1!)
    jump_2_2 = ConstantRateJump(rate_2_2, affect_2_2!)
    jump_2_3 = ConstantRateJump(rate_2_3, affect_2_3!)
    jump_2_4 = ConstantRateJump(rate_2_4, affect_2_4!)
    jump_2_5 = ConstantRateJump(rate_2_5, affect_2_5!)
    jump_2_6 = ConstantRateJump(rate_2_6, affect_2_6!)
    jump_2_7 = ConstantRateJump(rate_2_7, affect_2_7!)
    jumps_2 = (jump_2_1, jump_2_2, jump_2_3, jump_2_4, jump_2_5, jump_2_6, jump_2_7)
    push!(catalyst_networks, reaction_networks_hill[7])
    push!(manual_networks, jumps_2)
    push!(u0_syms, [:X1, :X2, :X3])
    push!(ps_syms, [:v, :K, :n, :k1, :k2, :k3, :d])
    push!(u0s, [:X1 => 10, :X2 => 10, :X3 => 10])
    push!(ps, [:v => 20.0, :K => 2.0, :n => 2, :k1 => 0.1, :k2 => 0.1, :k3 => 0.1, :d => 1.0])
    push!(sps, :X3)

    # Manual declaration of the third network (and its simulation settings).
    rate_3_1(u, p, t) = p[1] * binomial(u[1], 1)
    rate_3_2(u, p, t) = p[2] * binomial(u[2], 2)
    rate_3_3(u, p, t) = p[3] * binomial(u[2], 2)
    rate_3_4(u, p, t) = p[4] * binomial(u[3], 3)
    rate_3_5(u, p, t) = p[5] * binomial(u[3], 3)
    rate_3_6(u, p, t) = p[6] * binomial(u[4], 4)
    affect_3_1!(integrator) = (integrator.u[1] -= 1; integrator.u[2] += 2)
    affect_3_2!(integrator) = (integrator.u[2] -= 2; integrator.u[1] += 1)
    affect_3_3!(integrator) = (integrator.u[2] -= 2; integrator.u[3] += 3)
    affect_3_4!(integrator) = (integrator.u[3] -= 3; integrator.u[2] += 2)
    affect_3_5!(integrator) = (integrator.u[3] -= 3; integrator.u[4] += 4)
    affect_3_6!(integrator) = (integrator.u[4] -= 4; integrator.u[3] += 3)
    jump_3_1 = ConstantRateJump(rate_3_1, affect_3_1!)
    jump_3_2 = ConstantRateJump(rate_3_2, affect_3_2!)
    jump_3_3 = ConstantRateJump(rate_3_3, affect_3_3!)
    jump_3_4 = ConstantRateJump(rate_3_4, affect_3_4!)
    jump_3_5 = ConstantRateJump(rate_3_5, affect_3_5!)
    jump_3_6 = ConstantRateJump(rate_3_6, affect_3_6!)
    jumps_3 = (jump_3_1, jump_3_2, jump_3_3, jump_3_4, jump_3_5, jump_3_6)
    push!(catalyst_networks, reaction_networks_conserved[5])
    push!(manual_networks, jumps_3)
    push!(u0_syms, [:X1, :X2, :X3, :X4])
    push!(ps_syms, [:k1, :k2, :k3, :k4, :k5, :k6])
    push!(u0s, [:X1 => 15, :X2 => 5, :X3 => 5, :X4 => 5])
    push!(ps, [:k1 => 0.1, :k2 => 0.1, :k3 => 0.1, :k4 => 0.1, :k5 => 0.1, :k6 => 0.1])
    push!(sps, :X4)

    # Loops through all cases, checks that identical simulations are generated with/without Catalyst.
    for (rn_catalyst, rn_manual, u0_sym, ps_sym, u0_1, ps_1, sp) in 
            zip(catalyst_networks, manual_networks, u0_syms, ps_syms, u0s, ps, sps)

        # Simulates the Catalyst-created model.
        dprob_1 = DiscreteProblem(rn_catalyst, u0_1, (0.0, 10000.0), ps_1)
        jprob_1 = JumpProblem(rn_catalyst, dprob_1, Direct(); rng)
        sol1 = solve(jprob_1, SSAStepper(); seed, saveat = 1.0)
        
        # Simulates the manually written model
        u0_2 = map_to_vec(u0_1, u0_sym)
        ps_2 = map_to_vec(ps_1, ps_sym)
        dprob_2 = DiscreteProblem(u0_2, (0.0, 10000.0), ps_2)
        jprob_2 = JumpProblem(dprob_2, Direct(), rn_manual...; rng)
        sol2 = solve(jprob_2, SSAStepper(); seed, saveat = 1.0)
        
        # Checks that the means are similar (the test have been check that it holds across a large 
        # number of simulates, even without seed).
        @test mean(sol1[sp]) ≈ mean(sol2[findfirst(u0_sym .== sp),:]) rtol = 1e-1
    end
end

# Checks that simulations for a large number of potential systems are completed (and don't error).
let
    for rn in reaction_networks_all
        u0 = rnd_u0_Int64(rn, rng)
        ps = rnd_ps(rn, rng)
        dprob = DiscreteProblem(rn, u0, (0.0, 1.0), ps)
        jprob = JumpProblem(rn, dprob, Direct(); rng)
        @test SciMLBase.successful_retcode(solve(jprob, SSAStepper()))
    end
end

### Other Tests ###

# Tests simulating a network without parameters.
let
    no_param_network = @reaction_network begin 
        (1.2, 5), X1 ↔ X2 
    end
    u0 = rnd_u0_Int64(no_param_network, rng)
    dprob = DiscreteProblem(no_param_network, u0, (0.0, 1000.0))
    jprob = JumpProblem(no_param_network, dprob, Direct(); rng)
    sol = solve(jprob, SSAStepper())
    @test mean(sol[:X1]) > mean(sol[:X2])
end
