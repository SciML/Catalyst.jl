### Preparations ###

# Fetch packages.
using JumpProcesses
using Random, Statistics, SparseArrays, Test

# Fetch test networks.
include("../spatial_test_networks.jl")

### Correctness Tests ###

# Tests that there are no errors during runs.
let
    for grid in [small_2d_grid, short_path, small_directed_cycle]
        for srs in [Vector{TransportReaction}(), SIR_srs_1, SIR_srs_2]
            lrs = LatticeReactionSystem(SIR_system, srs, grid)
            u0_1 = [:S => 999, :I => 1, :R => 0]
            u0_2 = [:S => round.(Int64, 500.0 .+ 500.0 * rand_v_vals(lrs.lattice)), :I => 1, :R => 0, ]
            u0_3 = [:S => 950, :I => round.(Int64, 50 * rand_v_vals(lrs.lattice)), :R => round.(Int64, 50 * rand_v_vals(lrs.lattice))]
            u0_4 = [:S => round.(500.0 .+ 500.0 * rand_v_vals(lrs.lattice)), :I => round.(50 * rand_v_vals(lrs.lattice)), :R => round.(50 * rand_v_vals(lrs.lattice))]
            u0_5 = make_u0_matrix(u0_3, vertices(lrs.lattice), map(s -> Symbol(s.f), species(lrs.rs)))
            for u0 in [u0_1, u0_2, u0_3, u0_4, u0_5]
                p1 = [:α => 0.1 / 1000, :β => 0.01]
                p2 = [:α => 0.1 / 1000, :β => 0.02 * rand_v_vals(lrs.lattice)]
                p3 = [
                    :α => 0.1 / 2000 * rand_v_vals(lrs.lattice),
                    :β => 0.02 * rand_v_vals(lrs.lattice),
                ]
                p4 = make_u0_matrix(p1, vertices(lrs.lattice), Symbol.(parameters(lrs.rs)))
                for pV in [p1] #, p2, p3, p4] # Removed until spatial non-diffusion parameters are supported.
                    pE_1 = map(sp -> sp => 0.01, ModelingToolkit.getname.(edge_parameters(lrs)))
                    pE_2 = map(sp -> sp => 0.01, ModelingToolkit.getname.(edge_parameters(lrs)))
                    pE_3 = map(sp -> sp => rand_e_vals(lrs.lattice, 0.01), ModelingToolkit.getname.(edge_parameters(lrs)))
                    pE_4 = make_u0_matrix(pE_3, edges(lrs.lattice), ModelingToolkit.getname.(edge_parameters(lrs)))
                    for pE in [pE_1, pE_2, pE_3, pE_4]
                        println("\n\n\nHere: 1")
                        println(pE)
                        println(typeof(pE))
                        println((pV, pE))
                        dprob = DiscreteProblem(lrs, u0, (0.0, 100.0), (pV, pE))
                        println(dprob.p)
                        println(typeof(dprob.p))
                        jprob = JumpProblem(lrs, dprob, NSM())
                        @test SciMLBase.successful_retcode(solve(jprob, SSAStepper()))
                    end
                end
            end
        end
    end
end

### Input Handling Tests ###

# Tests that the correct hopping rates and initial conditions are generated.
# In this base case, hoppping rates should be on the form D_{s,i,j}.
let
    # Prepares the system.
    lrs = LatticeReactionSystem(SIR_system, SIR_srs_2, small_2d_grid)

    # Prepares various u0 input types.
    u0_1 = [:I => 2.0, :S => 1.0, :R => 3.0]
    u0_2 = [:I => fill(2., nv(small_2d_grid)), :S => 1.0, :R => 3.0]
    u0_3 = [1.0, 2.0, 3.0]
    u0_4 = [1.0, fill(2., nv(small_2d_grid)), 3.0]
    u0_5 = permutedims(hcat(fill(1., nv(small_2d_grid)), fill(2., nv(small_2d_grid)), fill(3., nv(small_2d_grid))))

    # Prepare various (compartment) parameter input types.
    pV_1 = [:β => 0.2, :α => 0.1]
    pV_2 = [:β => fill(0.2, nv(small_2d_grid)), :α => 1.0]
    pV_3 = [0.1, 0.2]
    pV_4 = [0.1, fill(0.2, nv(small_2d_grid))]
    pV_5 = permutedims(hcat(fill(0.1, nv(small_2d_grid)), fill(0.2, nv(small_2d_grid))))

    # Prepare various (diffusion) parameter input types.
    pE_1 = [:dI => 0.02, :dS => 0.01, :dR => 0.03]
    pE_2 = [:dI => 0.02, :dS => fill(0.01, ne(small_2d_grid)), :dR => 0.03]
    pE_3 = [0.01, 0.02, 0.03]
    pE_4 = [fill(0.01, ne(small_2d_grid)), 0.02, 0.03]
    pE_5 = permutedims(hcat(fill(0.01, ne(small_2d_grid)), fill(0.02, ne(small_2d_grid)), fill(0.03, ne(small_2d_grid))))

    # Checks hopping rates and u0 are correct.
    true_u0 = [fill(1.0, 1, 25); fill(2.0, 1, 25); fill(3.0, 1, 25)]
    true_hopping_rates = cumsum.([fill(dval, length(v)) for dval in [0.01,0.02,0.03], v in small_2d_grid.fadjlist])
    true_maj_scaled_rates = [0.1, 0.2]
    true_maj_reactant_stoch = [[1 => 1, 2 => 1], [2 => 1]]
    true_maj_net_stoch = [[1 => -1, 2 => 1], [2 => -1, 3 => 1]]
    for u0 in [u0_1, u0_2, u0_3, u0_4, u0_5]
        # Provides parameters as a tupple.
        for pV in [pV_1, pV_3], pE in [pE_1, pE_2, pE_3, pE_4, pE_5]
            dprob = DiscreteProblem(lrs, u0, (0.0, 100.0), (pV,pE))
            jprob = JumpProblem(lrs, dprob, NSM())
            @test jprob.prob.u0 == true_u0
            @test jprob.discrete_jump_aggregation.hop_rates.hop_const_cumulative_sums == true_hopping_rates
            @test jprob.massaction_jump.scaled_rates == true_maj_scaled_rates
            @test jprob.massaction_jump.reactant_stoch  == true_maj_reactant_stoch
            @test jprob.massaction_jump.net_stoch == true_maj_net_stoch
        end
        # Provides parameters as a combined vector.
        for pV in [pV_1], pE in [pE_1, pE_2]
            dprob = DiscreteProblem(lrs, u0, (0.0, 100.0), [pE; pV])
            jprob = JumpProblem(lrs, dprob, NSM())
            @test jprob.prob.u0 == true_u0
            @test jprob.discrete_jump_aggregation.hop_rates.hop_const_cumulative_sums == true_hopping_rates
            @test jprob.massaction_jump.scaled_rates == true_maj_scaled_rates
            @test jprob.massaction_jump.reactant_stoch  == true_maj_reactant_stoch
            @test jprob.massaction_jump.net_stoch == true_maj_net_stoch
        end
    end
end

### ABC Model Test (from JumpProcesses) ###
let 
    # Preparations (stuff used in JumpProcesses examples ported over here, could be written directly into code).
    Nsims = 100
    reltol = 0.05
    non_spatial_mean = [65.7395, 65.7395, 434.2605] #mean of 10,000 simulations
    dim = 1
    linear_size = 5
    num_nodes = linear_size^dim
    dims = Tuple(repeat([linear_size], dim))
    domain_size = 1.0 #μ-meter
    mesh_size = domain_size / linear_size
    rates = [0.1 / mesh_size, 1.0]
    diffusivity = 1.0
    num_species = 3

    # Make model.
    rn = @reaction_network begin
        (kB,kD), A + B <--> C
    end
    tr_1 = @transport_reaction D A
    tr_2 = @transport_reaction D B
    tr_3 = @transport_reaction D C
    lattice = Graphs.grid(dims)
    lrs = LatticeReactionSystem(rn, [tr_1, tr_2, tr_3], lattice)

    # Set simulation parameters and create problems
    u0 = [:A => [0,0,500,0,0], :B => [0,0,500,0,0], :C => 0]
    tspan = (0.0, 10.0)
    pV = [:kB => rates[1], :kD => rates[2]]
    pE = [:D => diffusivity]
    dprob = DiscreteProblem(lrs, u0, tspan, (pV, pE))
    jump_problems = [JumpProblem(lrs, dprob, alg(); save_positions = (false, false)) for alg in [NSM, DirectCRDirect]] # NRM doesn't work. Might need Cartesian grid.

    # Tests.
    function get_mean_end_state(jump_prob, Nsims)
        end_state = zeros(size(jump_prob.prob.u0))
        for i in 1:Nsims
            sol = solve(jump_prob, SSAStepper())
            end_state .+= sol.u[end]
        end
        end_state / Nsims
    end
    for jprob in jump_problems
        solution = solve(jprob, SSAStepper())
        mean_end_state = get_mean_end_state(jprob, Nsims)
        mean_end_state = reshape(mean_end_state, num_species, num_nodes)
        diff = sum(mean_end_state, dims = 2) - non_spatial_mean
        for (i, d) in enumerate(diff)
            @test abs(d) < reltol * non_spatial_mean[i]
        end
    end
end