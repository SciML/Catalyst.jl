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
        for srs in [Vector{DiffusionReaction}(), SIR_srs_1, SIR_srs_2]
            lrs = LatticeReactionSystem(SIR_system, srs, grid)
            u0_1 = make_values_int([:S => 999.0, :I => 1.0, :R => 0.0])
            u0_2 = make_values_int([
                                       :S => 500.0 .+ 500.0 * rand_v_vals(lrs.lattice),
                                       :I => 1.0,
                                       :R => 0.0,
                                   ])
            u0_3 = make_values_int([
                                       :S => 950.0,
                                       :I => 50 * rand_v_vals(lrs.lattice),
                                       :R => 50 * rand_v_vals(lrs.lattice),
                                   ])
            u0_4 = make_values_int([
                                       :S => 500.0 .+ 500.0 * rand_v_vals(lrs.lattice),
                                       :I => 50 * rand_v_vals(lrs.lattice),
                                       :R => 50 * rand_v_vals(lrs.lattice),
                                   ])
            u0_5 = make_values_int(make_u0_matrix(u0_3, vertices(lrs.lattice),
                                                  map(s -> Symbol(s.f), species(lrs.rs))))
            for u0 in [u0_1, u0_2, u0_3, u0_4, u0_5]
                p1 = [:α => 0.1 / 1000, :β => 0.01]
                p2 = [:α => 0.1 / 1000, :β => 0.02 * rand_v_vals(lrs.lattice)]
                p3 = [
                    :α => 0.1 / 2000 * rand_v_vals(lrs.lattice),
                    :β => 0.02 * rand_v_vals(lrs.lattice),
                ]
                p4 = make_u0_matrix(p1, vertices(lrs.lattice), Symbol.(parameters(lrs.rs)))
                for pV in [p1] #, p2, p3, p4] # Removed until spatial non-diffusion parameters are supported.
                    pE_1 = map(sp -> sp => 0.01,
                               ModelingToolkit.getname.(diffusion_parameters(lrs)))
                    pE_2 = map(sp -> sp => 0.01,
                               ModelingToolkit.getname.(diffusion_parameters(lrs)))
                    pE_3 = map(sp -> sp => rand_e_vals(lrs.lattice, 0.01),
                               ModelingToolkit.getname.(diffusion_parameters(lrs)))
                    pE_4 = make_u0_matrix(pE_3, edges(lrs.lattice),
                                          ModelingToolkit.getname.(diffusion_parameters(lrs)))
                    for pE in [pE_1, pE_2, pE_3, pE_4]
                        dprob = DiscreteProblem(lrs, u0, (0.0, 100.0), (pV, pE))
                        jprob = JumpProblem(lrs, dprob, NSM())
                        solve(jprob, SSAStepper())
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
    pC_1 = [:β => 0.2, :α => 0.1]
    pC_2 = [:β => fill(0.2, nv(small_2d_grid)), :α => 1.0]
    pC_3 = [0.1, 0.2]
    pC_4 = [0.1, fill(0.2, nv(small_2d_grid))]
    pC_5 = permutedims(hcat(fill(0.1, nv(small_2d_grid)), fill(0.2, nv(small_2d_grid))))

    # Prepare various (diffusion) parameter input types.
    pD_1 = [:dI => 0.02, :dS => 0.01, :dR => 0.03]
    pD_2 = [:dI => 0.02, :dS => fill(0.01, ne(small_2d_grid)), :dR => 0.03]
    pD_3 = [0.01, 0.02, 0.03]
    pD_4 = [fill(0.01, ne(small_2d_grid)), 0.02, 0.03]
    pD_5 = permutedims(hcat(fill(0.01, ne(small_2d_grid)), fill(0.02, ne(small_2d_grid)), fill(0.03, ne(small_2d_grid))))

    # Checks hopping rates and u0 are correct.
    true_u0 = [fill(1.0, 1, 25); fill(2.0, 1, 25); fill(3.0, 1, 25)]
    true_hopping_rates = cumsum.([fill(dval, length(v)) for dval in [0.01,0.02,0.03], v in small_2d_grid.fadjlist])
    for u0 in [u0_1, u0_2, u0_3, u0_4, u0_5]
        # Provides parameters as a tupple.
        for pC in [pC_1, pC_3], pD in [pD_1, pD_2, pD_3, pD_4, pD_5]
            dprob = DiscreteProblem(lrs, u0, (0.0, 100.0), (pC,pD))
            jprob = JumpProblem(lrs, dprob, NSM())
            @test jprob.prob.u0 == true_u0
            @test jprob.discrete_jump_aggregation.hop_rates.hop_const_cumulative_sums == true_hopping_rates
        end
        # Provides parameters as a combined vector.
        for pC in [pC_1], pD in [pD_1, pD_2]
            dprob = DiscreteProblem(lrs, u0, (0.0, 100.0), [pD; pC])
            jprob = JumpProblem(lrs, dprob, NSM())
            @test jprob.prob.u0 == true_u0
            @test jprob.discrete_jump_aggregation.hop_rates.hop_const_cumulative_sums == true_hopping_rates
        end
    end
end

### Hoping Rates Computation Tests ###

# Currently not in use, but will be added as more cases are enabled.
if false
    # Tests hopping rates of the form D_{s,i,j}
    let
        # Prepares the system.
        lrs = LatticeReactionSystem(SIR_system, SIR_srs_2, small_2d_grid)
        u0 = [:I => 2.0, :S => 1.0, :R => 3.0]
        pC = [:β => 0.2, :α => 0.1]

        # Prepare various (diffusion) parameter input types.
        pD_1 = [:dI => 0.02, :dS => 0.01, :dR => 0.03]
        pD_2 = [:dI => 0.02, :dS => fill(0.01, ne(lrs.lattice)), :dR => 0.03]
        pD_3 = [0.01, 0.02, 0.03]
        pD_4 = [fill(0.01, ne(lrs.lattice)), 0.02, 0.03]
        pD_5 = permutedims(hcat(fill(0.01, ne(lrs.lattice)), fill(0.02, ne(lrs.lattice)), fill(0.03, ne(lrs.lattice))))

        # Checks hopping rates are correct.
        true_hopping_rates = []
        # Provides parameters as a tupple.
        for pD in [pD_1, pD_2, pD_3, pD_4, pD_5]
            dprob = DiscreteProblem(lrs, u0, (0.0, 100.0), (pC,pD))
            jprob = JumpProblem(lrs, dprob, NSM())
            @test jprob.discrete_jump_aggregation.hop_rates.hop_const_cumulative_sums == true_hopping_rates
        end
        # Provides parameters as a combined vector.
        for pD in [pD_1, pD_2]
            dprob = DiscreteProblem(lrs, u0, (0.0, 100.0), [pD; pC])
            jprob = JumpProblem(lrs, dprob, NSM())
            @test jprob.discrete_jump_aggregation.hop_rates.hop_const_cumulative_sums == true_hopping_rates
        end
    end

    # Tests hopping rates of the form D_{s}
    let
        # Prepares the system.
        lrs = LatticeReactionSystem(SIR_system, SIR_srs_2, small_2d_grid)
        u0 = [:I => 2.0, :S => 1.0, :R => 3.0]
        pC = [:β => 0.2, :α => 0.1]

        # Prepare various (diffusion) parameter input types.
        pD_1 = [:dI => 0.02, :dS => 0.01, :dR => 0.03]
        pD_3 = [0.01, 0.02, 0.03]

        # Checks hopping rates are correct.
        true_hopping_rates = []
        # Provides parameters as a tupple.
        for pD in [pD_1, pD_3]
            dprob = DiscreteProblem(lrs, u0, (0.0, 100.0), (pC,pD))
            jprob = JumpProblem(lrs, dprob, NSM())
            @test jprob.discrete_jump_aggregation.hop_rates.hop_const_cumulative_sums == true_hopping_rates
        end
        # Provides parameters as a combined vector.
        for pD in [pD_1]
            dprob = DiscreteProblem(lrs, u0, (0.0, 100.0), [pD; pC])
            jprob = JumpProblem(lrs, dprob, NSM())
            @test jprob.discrete_jump_aggregation.hop_rates.hop_const_cumulative_sums == true_hopping_rates
        end
    end

    # Tests hopping rates of the form D_{s,i}
    let
        # Prepares special system and diffusion reactions (depending on compartments).
        SIR_system_special = @reaction_network begin
            @parameters comp_dS comp_dI
            α, S + I --> 2I
            β, I --> R
        end
        @unpack comp_dS, comp_dI = SIR_system_special
        SIR_srs_special = [DiffusionReaction(comp_dS, :S), DiffusionReaction(comp_dI, :I)]

        # Prepares the system.
        lrs = LatticeReactionSystem(SIR_system_special, SIR_srs_special, small_2d_grid)
        u0 = [:I => 2.0, :S => 1.0, :R => 3.0]
        pD = []

        # Prepare various (compartment) parameter input types.
        pC_1 = [:β => 0.2, :α => 0.1, :comp_dI => 0.02, :comp_dS => 0.01]
        pC_2 = [:β => fill(0.2, nv(small_2d_grid)), :α => 1.0, :comp_dI => fill(0.02, nv(small_2d_grid)), :comp_dS => 0.01]
        pC_3 = [:β => 0.2, :α => 1.0, :comp_dI => fill(0.02, nv(small_2d_grid)), :comp_dS => 0.01]
        pC_4 = [0.1, 0.2, 0.01, 0.02]
        pC_5 = [0.1, 0.2, 0.01, fill(0.02, nv(small_2d_grid))]
        pC_6 = [0.1, fill(0.2, nv(small_2d_grid)), 0.01, fill(0.02, nv(small_2d_grid))]
        pC_7 = permutedims(hcat(fill(0.1, nv(small_2d_grid)), fill(0.2, nv(small_2d_grid)), fill(0.01, nv(small_2d_grid)), fill(0.02, nv(small_2d_grid))))

        # Checks hopping rates are correct.
        true_hopping_rates = []
        # Provides parameters as a tupple.
        for pC in [pC_1, pC_2, pC_4, pC_5]
            dprob = DiscreteProblem(lrs, u0, (0.0, 100.0), (pC,pD))
            jprob = JumpProblem(lrs, dprob, NSM())
            @test jprob.discrete_jump_aggregation.hop_rates.hop_const_cumulative_sums == true_hopping_rates
        end
        # Provides parameters as a combined vector.
        for pC in [pC_1, pC_2]
            dprob = DiscreteProblem(lrs, u0, (0.0, 100.0), [pD; pC])
            jprob = JumpProblem(lrs, dprob, NSM())
            @test jprob.discrete_jump_aggregation.hop_rates.hop_const_cumulative_sums == true_hopping_rates
        end
    end

    # Tests hopping rates of the form D_{s}*L_{i,j}
    # Since dS and dI are constant across the system, these can both be diffusionparameter or not (bot hcases tested).
    let
        SIR_system_special = @reaction_network begin
            @parameters dS [diffusionparameter = true] dI connection_d [diffusionparameter = true] 
            α, S + I --> 2I
            β, I --> R
        end
        @unpack dS, dI = SIR_system_special
        SIR_srs_special = [DiffusionReaction(dS*connection_d, :S), DiffusionReaction(dI*connection_d, :I)]

        # Prepares the system.
        lrs = LatticeReactionSystem(SIR_system_special, SIR_srs_special, small_2d_grid)
        u0 = [:I => 2.0, :S => 1.0, :R => 3.0]
        pC = [:β => 0.2, :α => 0.1, :dI => 0.01]

        # Prepare various (compartment) parameter input types.
        pC_1 = [:β => 0.2, :α => 0.1, ]
        pC_2 = [:β => fill(0.2, nv(small_2d_grid)), :α => 1.0]
        pC_3 = [0.1, 0.2]
        pC_4 = [0.1, fill(0.2, nv(small_2d_grid))]
        pC_5 = permutedims(hcat(fill(0.1, nv(small_2d_grid)), fill(0.2, nv(small_2d_grid))))

        # Prepare various (diffusion) parameter input types.
        pD_1 = [:connection_d => 2.0, :dS => 0.005]
        pD_2 = [:connection_d => fill(0.2, ne(small_2d_grid)), :dS => 0.005]
        pD_3 = [0.005, 2.0]
        pD_4 = [0.005, fill(0.2, ne(small_2d_grid))]

        # Checks hopping rates are correct.
        true_hopping_rates = []
        # Provides parameters as a tupple.
        for pD in [pD_1, pD_2, pD_3, pD_4]
            dprob = DiscreteProblem(lrs, u0, (0.0, 100.0), (pC,pD))
            jprob = JumpProblem(lrs, dprob, NSM())
            @test jprob.discrete_jump_aggregation.hop_rates.hop_const_cumulative_sums == true_hopping_rates
        end
        # Provides parameters as a combined vector.
        for pD in [pD_1, pD_2]
            dprob = DiscreteProblem(lrs, u0, (0.0, 100.0), [pD; pC])
            jprob = JumpProblem(lrs, dprob, NSM())
            @test jprob.discrete_jump_aggregation.hop_rates.hop_const_cumulative_sums == true_hopping_rates
        end
    end

    # Tests hopping rates of the form D_{s,i}*L_{i,j}
    let
        SIR_system_special = @reaction_network begin
            @parameters comp_dS comp_dI connection_d [diffusionparameter = true] 
            α, S + I --> 2I
            β, I --> R
        end
        @unpack comp_dS, comp_dI, connection_d = SIR_system_special
        SIR_srs_special = [DiffusionReaction(comp_dS*connection_d, :S), DiffusionReaction(comp_dI*connection_d, :I)]

        # Prepares the system.
        lrs = LatticeReactionSystem(SIR_system_special, SIR_srs_special, small_2d_grid)
        u0 = [:I => 2.0, :S => 1.0, :R => 3.0]

        # Prepare various (compartment) parameter input types.
        pC_1 = [:β => 0.2, :α => 0.1, :comp_dI => 0.02, :comp_dS => 0.01]
        pC_2 = [:β => fill(0.2, nv(small_2d_grid)), :α => 1.0, :comp_dI => fill(0.02, nv(small_2d_grid)), :comp_dS => 0.01]
        pC_3 = [:β => 0.2, :α => 1.0, :comp_dI => fill(0.02, nv(small_2d_grid)), :comp_dS => 0.01]
        pC_4 = [0.1, 0.2, 0.01, 0.02]
        pC_5 = [0.1, 0.2, 0.01, fill(0.02, nv(small_2d_grid))]
        pC_6 = [0.1, fill(0.2, nv(small_2d_grid)), 0.01, fill(0.02, nv(small_2d_grid))]
        pC_7 = permutedims(hcat(fill(0.1, nv(small_2d_grid)), fill(0.2, nv(small_2d_grid)), fill(0.01, nv(small_2d_grid)), fill(0.02, nv(small_2d_grid))))

        # Prepare various (diffusion) parameter input types.
        pD_1 = [:connection_d => 2.0]
        pD_2 = [:connection_d => fill(0.2, ne(small_2d_grid))]
        pD_3 = [2.0]
        pD_4 = [fill(0.2, ne(small_2d_grid))]
        pD_5 = [fill(0.2, 1, ne(small_2d_grid))]

        # Checks hopping rates are correct.
        true_hopping_rates = []
        # Provides parameters as a tupple.
        for pC in [pC_1, pC_2, pC_4, pC_5], pD in [pD_1, pD_2, pD_3, pD_4, pD_5]
            dprob = DiscreteProblem(lrs, u0, (0.0, 100.0), (pC,pD))
            jprob = JumpProblem(lrs, dprob, NSM())
            @test jprob.discrete_jump_aggregation.hop_rates.hop_const_cumulative_sums == true_hopping_rates
        end
        # Provides parameters as a combined vector.
        for pC in [pC_1, pC_2], pD in [pD_1, pD_2]
            dprob = DiscreteProblem(lrs, u0, (0.0, 100.0), [pD; pC])
            jprob = JumpProblem(lrs, dprob, NSM())
            @test jprob.discrete_jump_aggregation.hop_rates.hop_const_cumulative_sums == true_hopping_rates
        end
    end
end