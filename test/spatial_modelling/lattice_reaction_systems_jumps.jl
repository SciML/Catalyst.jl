### Preparations ###

# Fetch packages.
using JumpProcesses, Statistics, SparseArrays, Test

# Fetch test networks.
include("../spatial_test_networks.jl")


### General Tests ###

# Tests that there are no errors during runs for a variety of input forms.
let
    for grid in [small_2d_graph_grid, small_2d_cartesian_grid, small_2d_masked_grid]
        for srs in [Vector{TransportReaction}(), SIR_srs_1, SIR_srs_2]
            lrs = LatticeReactionSystem(SIR_system, srs, grid)
            u0_1 = [:S => 999, :I => 1, :R => 0]
            u0_2 = [:S => round.(Int64, 500 .+ 500 * rand_v_vals(lrs)), :I => 1, :R => 0]
            u0_3 = [
                :S => round.(Int64, 500 .+ 500 * rand_v_vals(lrs)),
                :I => round.(Int64, 50 * rand_v_vals(lrs)),
                :R => round.(Int64, 50 * rand_v_vals(lrs)),
            ]
            for u0 in [u0_1, u0_2, u0_3]
                pV_1 = [:α => 0.1 / 1000, :β => 0.01]
                pV_2 = [:α => 0.1 / 1000, :β => 0.02 * rand_v_vals(lrs)]
                pV_3 = [
                    :α => 0.1 / 2000 * rand_v_vals(lrs),
                    :β => 0.02 * rand_v_vals(lrs),
                ]
                for pV in [pV_1, pV_2, pV_3]
                    pE_1 = [sp => 0.01 for sp in spatial_param_syms(lrs)]
                    pE_2 = [sp => rand_e_vals(lrs) / 50.0 for sp in spatial_param_syms(lrs)]
                    for pE in [pE_1, pE_2]
                        isempty(spatial_param_syms(lrs)) && (pE = Vector{Pair{Symbol, Float64}}())
                        dprob = DiscreteProblem(lrs, u0, (0.0, 1.0), [pV; pE])
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
# In this base case, hopping rates should be on the form D_{s,i,j}.
let
    # Prepares the system.
    grid = small_2d_graph_grid
    lrs = LatticeReactionSystem(SIR_system, SIR_srs_2, grid)

    # Prepares various u0 input types.
    u0_1 = [:I => 2.0, :S => 1.0, :R => 3.0]
    u0_2 = [:I => fill(2.0, nv(grid)), :S => 1.0, :R => 3.0]

    # Prepare various (compartment) parameter input types.
    pV_1 = [:β => 0.2, :α => 0.1]
    pV_2 = [:β => fill(0.2, nv(grid)), :α => 1.0]

    # Prepare various (diffusion) parameter input types.
    pE_1 = [:dI => 0.02, :dS => 0.01, :dR => 0.03]
    dS_vals = spzeros(num_verts(lrs), num_verts(lrs))
    foreach(e -> (dS_vals[e[1], e[2]] = 0.01), edge_iterator(lrs))
    pE_2 = [:dI => 0.02, :dS => dS_vals, :dR => 0.03]

    # Checks hopping rates and u0 are correct.
    true_u0 = [fill(1.0, 1, 25); fill(2.0, 1, 25); fill(3.0, 1, 25)]
    true_hopping_rates = cumsum.([fill(dval, length(v)) for dval in [0.01, 0.02, 0.03], v in grid.fadjlist])
    true_maj_scaled_rates = [0.1, 0.2]
    true_maj_reactant_stoch = [[1 => 1, 2 => 1], [2 => 1]]
    true_maj_net_stoch = [[1 => -1, 2 => 1], [2 => -1, 3 => 1]]
    for u0 in [u0_1, u0_2]
        for pV in [pV_1, pV_2], pE in [pE_1, pE_2]
            dprob = DiscreteProblem(lrs, u0, (0.0, 100.0), [pE; pV])
            jprob = JumpProblem(lrs, dprob, NSM())
            @test jprob.prob.u0 == true_u0
            @test jprob.discrete_jump_aggregation.hop_rates.hop_const_cumulative_sums == true_hopping_rates
            @test jprob.massaction_jump.reactant_stoch == true_maj_reactant_stoch
            @test all(issetequal(ns1, ns2) for (ns1, ns2) in zip(jprob.massaction_jump.net_stoch, true_maj_net_stoch))
        end
    end
end


### SpatialMassActionJump Testing ###

# Checks that the correct structures are produced.
let
    # Network for reference:
    # A, ∅ → X
    # 1, 2X + Y → 3X
    # B, X → Y
    # 1, X → ∅
    # srs = [@transport_reaction dX X]
    # Create LatticeReactionSystem
    lrs = LatticeReactionSystem(brusselator_system, brusselator_srs_1, small_3d_graph_grid)

    # Create JumpProblem
    u0 = [:X => 1, :Y => rand(1:10, num_verts(lrs))]
    tspan = (0.0, 100.0)
    ps = [:A => 1.0, :B => 5.0 .+ rand_v_vals(lrs), :dX => rand_e_vals(lrs)]
    dprob = DiscreteProblem(lrs, u0, tspan, ps)
    jprob = JumpProblem(lrs, dprob, NSM())

    # Checks internal structures.
    jprob.massaction_jump.uniform_rates == [1.0, 0.5, 10.0] # 0.5 is due to combinatoric /2! in (2X + Y).
    jprob.massaction_jump.spatial_rates[1, :] == ps[2][2]
    # Test when new SII functions are ready, or we implement them in Catalyst.
    # @test isequal(to_int(getfield.(reactions(reactionsystem(lrs)), :netstoich)), jprob.massaction_jump.net_stoch)
    # @test isequal(to_int(Pair.(getfield.(reactions(reactionsystem(lrs)), :substrates),getfield.(reactions(reactionsystem(lrs)), :substoich))), jprob.massaction_jump.net_stoch)

    # Checks that problems can be simulated.
    @test SciMLBase.successful_retcode(solve(jprob, SSAStepper()))
end

# Checks that heterogeneous vertex parameters work. Checks that birth-death system with different
# birth rates produce different means.
let
    # Create model.
    birth_death_network = @reaction_network begin
        (p, d), 0 <--> X
    end
    srs = [(@transport_reaction D X)]
    lrs = LatticeReactionSystem(birth_death_network, srs, very_small_2d_graph_grid)

    # Create JumpProblem.
    u0 = [:X => 1]
    tspan = (0.0, 100.0)
    ps = [:p => [0.1, 1.0, 10.0, 100.0], :d => 1.0, :D => 0.0]
    dprob = DiscreteProblem(lrs, u0, tspan, ps)
    jprob = JumpProblem(lrs, dprob, NSM())

    # Simulate model (a few repeats to ensure things don't succeed by change for uniform rates).
    # Check that higher p gives higher mean.
    for i in 1:5
        sol = solve(jprob, SSAStepper(); saveat = 1.0)
        @test mean(getindex.(sol.u, 1)) < mean(getindex.(sol.u, 2)) < mean(getindex.(sol.u, 3)) < mean(getindex.(sol.u, 4))
    end
end


### Tests taken from JumpProcesses ###

# ABC Model Test
let
    # Preparations (stuff used in JumpProcesses examples ported over here, could be written directly into code).
    Nsims = 100
    reltol = 0.05
    non_spatial_mean = [65.7395, 65.7395, 434.2605] # Mean of 10,000 simulations.
    dim = 1
    linear_size = 5
    num_nodes = linear_size^dim
    dims = Tuple(repeat([linear_size], dim))
    domain_size = 1.0 # μ-meter.
    mesh_size = domain_size / linear_size
    rates = [0.1 / mesh_size, 1.0]
    diffusivity = 1.0
    num_species = 3

    # Make model.
    rn = @reaction_network begin
        (kB, kD), A + B <--> C
    end
    tr_1 = @transport_reaction D A
    tr_2 = @transport_reaction D B
    tr_3 = @transport_reaction D C
    lattice = Graphs.grid(dims)
    lrs = LatticeReactionSystem(rn, [tr_1, tr_2, tr_3], lattice)

    # Set simulation parameters and create problems.
    u0 = [:A => [0, 0, 500, 0, 0], :B => [0, 0, 500, 0, 0], :C => 0]
    tspan = (0.0, 10.0)
    pV = [:kB => rates[1], :kD => rates[2]]
    pE = [:D => diffusivity]
    dprob = DiscreteProblem(lrs, u0, tspan, [pV; pE])
    # NRM could be added, but doesn't work. Might need Cartesian grid.
    jump_problems = [JumpProblem(lrs, dprob, alg(); save_positions = (false, false)) for alg in [NSM, DirectCRDirect]]

    # Run tests.
    function get_mean_end_state(jump_prob, Nsims)
        end_state = zeros(size(jump_prob.prob.u0))
        for i in 1:Nsims
            sol = solve(jump_prob, SSAStepper())
            end_state .+= sol.u[end]
        end
        return end_state / Nsims
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


### Other Tests ###

# Checks that providing a non-spatial `DiscreteProblem` to a `JumpProblem` gives an error.
let
    lrs = LatticeReactionSystem(binding_system, binding_srs, very_small_2d_masked_grid)
    dprob = DiscreteProblem(binding_system, binding_u0, (0.0, 10.0), binding_p[1:2])
    @test_throws ArgumentError JumpProblem(lrs, dprob, NSM())
end
