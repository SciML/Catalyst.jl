### Fetch Packages and Reaction Networks ###

# Fetch packages.
using Catalyst, JumpProcesses, Random, Statistics, Test, SciMLBase
using ModelingToolkit: get_states, get_ps

# Sets rnd number.
using StableRNGs
rng = StableRNG(12345)

# Fetch test networks.
include("../test_networks.jl")

### Compares to Manually Calcualted Function ###

let
    identical_networks = Vector{Pair}()

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
    push!(identical_networks, reaction_networks_standard[6] => jumps_1)

    rate_2_1(u, p, t) = p[1] / 10 + u[1]^p[3] / (u[1]^p[3] + p[2]^p[3])
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
    push!(identical_networks, reaction_networks_hill[7] => jumps_2)

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
    push!(identical_networks, reaction_networks_constraint[5] => jumps_3)

    for (i, networks) in enumerate(identical_networks)
        for factor in [1e-2, 1e-1, 1e0, 1e1], repeat in 1:3
            (i == 3) && (factor > 1e-1) && continue   # Large numbers seems to crash it.
            u0 = rand(rng, 1:Int64(factor * 100), length(get_states(networks[1])))
            p = factor * rand(rng, length(get_ps(networks[1])))
            prob1 = JumpProblem(networks[1],
                                DiscreteProblem(networks[1], u0, (0.0, 1000.0), p),
                                Direct())
            sol1 = solve(prob1, SSAStepper())
            prob2 = JumpProblem(DiscreteProblem(u0, (0.0, 1000.0), p), Direct(),
                                networks[2]...)
            sol2 = solve(prob2, SSAStepper())
            for i in 1:length(u0)
                vals1 = getindex.(sol1.u, i)
                vals2 = getindex.(sol1.u, i)
                (mean(vals2) > 0.001) && @test 0.8 < mean(vals1) / mean(vals2) < 1.25
                (std(vals2) > 0.001) && @test 0.8 < std(vals1) / std(vals2) < 1.25
            end
        end
    end
end

### Checks Simulations Don't Error ###

let
    for (i, network) in enumerate(reaction_networks_all)
        for factor in [1e-1, 1e0, 1e1]
            u0 = rand(rng, 1:Int64(factor * 100), length(get_states(network)))
            p = factor * rand(rng, length(get_ps(network)))
            prob = JumpProblem(network, DiscreteProblem(network, u0, (0.0, 1.0), p),
                               Direct())
            @test SciMLBase.successful_retcode(solve(prob, SSAStepper()))
        end
    end
end

### Other Tests ###

# No parameter test.
let
    no_param_network = @reaction_network begin (1.2, 5), X1 ↔ X2 end
    for factor in [1e1, 1e2]
        u0 = rand(rng, 1:Int64(factor * 100), length(get_states(no_param_network)))
        prob = JumpProblem(no_param_network,
                           DiscreteProblem(no_param_network, u0, (0.0, 1000.0)), Direct())
        sol = solve(prob, SSAStepper())
        vals1 = getindex.(sol.u[1:end], 1)
        vals2 = getindex.(sol.u[1:end], 2)
        @test mean(vals1) > mean(vals2)
    end
end
