### Preparations ###

# Fetch packages.
using OrdinaryDiffEq
using Random, Statistics, SparseArrays, Test

# Sets rnd number.
using StableRNGs
rng = StableRNG(12345)

# Fetch test networks.
include("spatial_test_networks.jl")


### Test No Error During Runs ###
for grid in [small_2d_grid, short_path, small_directed_cycle]
    # Non-stiff case
    for srs in [Vector{DiffusionReaction}(), SIR_srs_1, SIR_srs_2]
        lrs = LatticeReactionSystem(SIR_system, srs, grid)
        u0_1 = [:S => 999.0, :I => 1.0, :R => 0.0]
        u0_2 = [:S => 500.0 .+ 500.0 * rand_v_vals(lrs.lattice), :I => 1.0, :R => 0.0]
        u0_3 = [
            :S => 950.0,
            :I => 50 * rand_v_vals(lrs.lattice),
            :R => 50 * rand_v_vals(lrs.lattice),
        ]
        u0_4 = [
            :S => 500.0 .+ 500.0 * rand_v_vals(lrs.lattice),
            :I => 50 * rand_v_vals(lrs.lattice),
            :R => 50 * rand_v_vals(lrs.lattice),
        ]
        u0_5 = make_u0_matrix(u0_3, vertices(lrs.lattice),
                              map(s -> Symbol(s.f), species(lrs.rs)))
        for u0 in [u0_1, u0_2, u0_3, u0_4, u0_5]
            p1 = [:α => 0.1 / 1000, :β => 0.01]
            p2 = [:α => 0.1 / 1000, :β => 0.02 * rand_v_vals(lrs.lattice)]
            p3 = [
                :α => 0.1 / 2000 * rand_v_vals(lrs.lattice),
                :β => 0.02 * rand_v_vals(lrs.lattice),
            ]
            p4 = make_u0_matrix(p1, vertices(lrs.lattice), Symbol.(parameters(lrs.rs)))
            for pV in [p1, p2, p3, p4]
                pE_1 = map(sp -> sp => 0.01, spatial_param_syms(lrs))
                pE_2 = map(sp -> sp => 0.01, spatial_param_syms(lrs))
                pE_3 = map(sp -> sp => rand_e_vals(lrs.lattice, 0.01),
                           spatial_param_syms(lrs))
                pE_4 = make_u0_matrix(pE_3, edges(lrs.lattice), spatial_param_syms(lrs))
                for pE in [pE_1, pE_2, pE_3, pE_4]
                    oprob = ODEProblem(lrs, u0, (0.0, 500.0), (pV, pE))
                    @test SciMLBase.successful_retcode(solve(oprob, Tsit5()))

                    oprob = ODEProblem(lrs, u0, (0.0, 10.0), (pV, pE); jac = false)
                    @test SciMLBase.successful_retcode(solve(oprob, Tsit5()))
                end
            end
        end
    end

    # Stiff case
    for srs in [Vector{DiffusionReaction}(), brusselator_srs_1, brusselator_srs_2]
        lrs = LatticeReactionSystem(brusselator_system, srs, grid)
        u0_1 = [:X => 1.0, :Y => 20.0]
        u0_2 = [:X => rand_v_vals(lrs.lattice, 10.0), :Y => 2.0]
        u0_3 = [:X => rand_v_vals(lrs.lattice, 20), :Y => rand_v_vals(lrs.lattice, 10)]
        u0_4 = make_u0_matrix(u0_3, vertices(lrs.lattice),
                              map(s -> Symbol(s.f), species(lrs.rs)))
        for u0 in [u0_1, u0_2, u0_3, u0_4]
            p1 = [:A => 1.0, :B => 4.0]
            p2 = [:A => 0.5 .+ rand_v_vals(lrs.lattice, 0.5), :B => 4.0]
            p3 = [
                :A => 0.5 .+ rand_v_vals(lrs.lattice, 0.5),
                :B => 4.0 .+ rand_v_vals(lrs.lattice, 1.0),
            ]
            p4 = make_u0_matrix(p2, vertices(lrs.lattice), Symbol.(parameters(lrs.rs)))
            for pV in [p1, p2, p3, p4]
                pE_1 = map(sp -> sp => 0.2, spatial_param_syms(lrs))
                pE_2 = map(sp -> sp => rand(), spatial_param_syms(lrs))
                pE_3 = map(sp -> sp => rand_e_vals(lrs.lattice, 0.2),
                           spatial_param_syms(lrs))
                pE_4 = make_u0_matrix(pE_3, edges(lrs.lattice), spatial_param_syms(lrs))
                for pE in [pE_1, pE_2, pE_3, pE_4]
                    oprob = ODEProblem(lrs, u0, (0.0, 10.0), (pV, pE))
                    @test SciMLBase.successful_retcode(solve(oprob, QNDF()))

                    oprob = ODEProblem(lrs, u0, (0.0, 10.0), (pV, pE); sparse = false)
                    @test SciMLBase.successful_retcode(solve(oprob, QNDF()))
                end
            end
        end
    end
end

### Tests Simulation Correctness ###

# Checks that non-spatial brusselator simulation is identical to all on an unconnected lattice.
let
    lrs = LatticeReactionSystem(brusselator_system, brusselator_srs_1, unconnected_graph)
    u0 = [:X => 2.0 + 2.0 * rand(), :Y => 10.0 + 10.0 * rand()]
    pV = brusselator_p
    pE = [:dX => 0.2]
    oprob_nonspatial = ODEProblem(brusselator_system, u0, (0.0, 100.0), pV)
    oprob_spatial = ODEProblem(lrs, u0, (0.0, 100.0), (pV, pE))
    sol_nonspatial = solve(oprob_nonspatial, QNDF(); abstol = 1e-12, reltol = 1e-12)
    sol_spatial = solve(oprob_spatial, QNDF(); abstol = 1e-12, reltol = 1e-12)

    for i in 1:nv(unconnected_graph)
        @test all(isapprox.(sol_nonspatial.u[end],
                            sol_spatial.u[end][((i - 1) * 2 + 1):((i - 1) * 2 + 2)]))
    end
end

# Checks that result becomes homogeneous on a connected lattice.
let
    lrs = LatticeReactionSystem(binding_system, binding_srs, undirected_cycle)
    u0 = [
        :X => 1.0 .+ rand_v_vals(lrs.lattice),
        :Y => 2.0 * rand_v_vals(lrs.lattice),
        :XY => 0.5,
    ]
    oprob = ODEProblem(lrs, u0, (0.0, 1000.0), binding_p; tstops = 0.1:0.1:1000.0)
    ss = solve(oprob, Tsit5()).u[end]

    @test all(isapprox.(ss[1:3:end], ss[1]))
    @test all(isapprox.(ss[2:3:end], ss[2]))
    @test all(isapprox.(ss[3:3:end], ss[3]))
end

### Tests Special Cases ###

# Creates network with various combiantions of Symbls and Nums in diffusion reactions.
let
    @parameters dS dI dR
    @variables t
    @species S(t) I(t) R(t)
    SIR_srs_numsym_1 = diffusion_reactions([(:dS, :S), (:dI, :I), (:dR, :R)])
    SIR_srs_numsym_2 = diffusion_reactions([(dS, :S), (dI, :I), (dR, :R)])
    SIR_srs_numsym_3 = diffusion_reactions([(:dS, S), (:dI, I), (:dR, R)])
    SIR_srs_numsym_4 = diffusion_reactions([(dS, S), (dI, I), (dR, R)])
    SIR_srs_numsym_5 = diffusion_reactions([(dS, :S), (:dI, I), (dR, :R)])
    SIR_srs_numsym_6 = diffusion_reactions([(:dS, :S), (:dI, I), (dR, R)])

    u0 = [:S => 990.0, :I => 20.0 * rand_v_vals(small_2d_grid), :R => 0.0]
    pV = SIR_p
    pE_1 = [:dS => 0.01, :dI => 0.01, :dR => 0.01]
    pE_2 = [dS => 0.01, dI => 0.01, dR => 0.01]
    pE_3 = [dS => 0.01, :dI => 0.01, :dR => 0.01]
    ss_explicit_base = solve(ODEProblem(LatticeReactionSystem(SIR_system, SIR_srs_numsym_1, small_2d_grid), u0, (0.0, 10.0), (pV, pE_1); jac = false), Tsit5()).u[end]
    ss_implicit_base = solve(ODEProblem(LatticeReactionSystem(SIR_system, SIR_srs_numsym_1, small_2d_grid), u0, (0.0, 10.0), (pV, pE_1); jac = true), Rosenbrock23()).u[end]

    for srs in [
            SIR_srs_numsym_1,
            SIR_srs_numsym_2,
            SIR_srs_numsym_3,
            SIR_srs_numsym_4,
            SIR_srs_numsym_5,
            SIR_srs_numsym_6,
        ], pE in [pE_1, pE_2, pE_3]
        lrs = LatticeReactionSystem(SIR_system, srs, small_2d_grid)
        ss_explicit = solve(ODEProblem(lrs, u0, (0.0, 10.0), (pV, pE); jac = false), Tsit5()).u[end]
        ss_implicit = solve(ODEProblem(lrs, u0, (0.0, 10.0), (pV, pE); jac = true), Rosenbrock23()).u[end]
        @test all(isapprox.(ss_explicit, ss_explicit_base))
        @test all(isapprox.(ss_implicit, ss_implicit_base))
    end
end

# Create network with vaious combinations of graph/di-graph and parameters.
let
    lrs_digraph = LatticeReactionSystem(SIR_system, SIR_srs_2, complete_digraph(3))
    lrs_graph = LatticeReactionSystem(SIR_system, SIR_srs_2, complete_graph(3))
    u0 = [:S => 990.0, :I => 20.0 * rand_v_vals(lrs_digraph.lattice), :R => 0.0]
    pV = SIR_p
    pE_digraph_1 = [:dS => [0.10, 0.10, 0.12, 0.12, 0.14, 0.14], :dI => 0.01, :dR => 0.01]
    pE_digraph_2 = [[0.10, 0.10, 0.12, 0.12, 0.14, 0.14], 0.01, 0.01]
    pE_digraph_3 = [0.10 0.10 0.12 0.12 0.14 0.14; 0.01 0.01 0.01 0.01 0.01 0.01;
                    0.01 0.01 0.01 0.01 0.01 0.01]
    pE_graph_1 = [:dS => [0.10, 0.12, 0.14], :dI => 0.01, :dR => 0.01]
    pE_graph_2 = [[0.10, 0.12, 0.14], 0.01, 0.01]
    pE_graph_3 = [0.10 0.12 0.14; 0.01 0.01 0.01; 0.01 0.01 0.01]
    oprob_digraph_1 = ODEProblem(lrs_digraph, u0, (0.0, 500.0), (pV, pE_digraph_1))
    oprob_digraph_2 = ODEProblem(lrs_digraph, u0, (0.0, 500.0), (pV, pE_digraph_2))
    oprob_digraph_3 = ODEProblem(lrs_digraph, u0, (0.0, 500.0), (pV, pE_digraph_3))
    oprob_graph_11 = ODEProblem(lrs_graph, u0, (0.0, 500.0), (pV, pE_digraph_1))
    oprob_graph_12 = ODEProblem(lrs_graph, u0, (0.0, 500.0), (pV, pE_graph_1))
    oprob_graph_21 = ODEProblem(lrs_graph, u0, (0.0, 500.0), (pV, pE_digraph_2))
    oprob_graph_22 = ODEProblem(lrs_graph, u0, (0.0, 500.0), (pV, pE_graph_2))
    oprob_graph_31 = ODEProblem(lrs_graph, u0, (0.0, 500.0), (pV, pE_digraph_3))
    oprob_graph_32 = ODEProblem(lrs_graph, u0, (0.0, 500.0), (pV, pE_graph_3))
    sim_end_digraph_1 = solve(oprob_digraph_1, Tsit5()).u[end]
    sim_end_digraph_2 = solve(oprob_digraph_2, Tsit5()).u[end]
    sim_end_digraph_3 = solve(oprob_digraph_3, Tsit5()).u[end]
    sim_end_graph_11 = solve(oprob_graph_11, Tsit5()).u[end]
    sim_end_graph_12 = solve(oprob_graph_12, Tsit5()).u[end]
    sim_end_graph_21 = solve(oprob_graph_21, Tsit5()).u[end]
    sim_end_graph_22 = solve(oprob_graph_22, Tsit5()).u[end]
    sim_end_graph_31 = solve(oprob_graph_31, Tsit5()).u[end]
    sim_end_graph_32 = solve(oprob_graph_32, Tsit5()).u[end]

    @test all(sim_end_digraph_1 .== sim_end_digraph_2 .== sim_end_digraph_3 .==
              sim_end_graph_11 .== sim_end_graph_12 .== sim_end_graph_21 .==
              sim_end_graph_22 .== sim_end_graph_31 .== sim_end_graph_32)
end

# Creates networks with empty species or parameters.
let
    binding_system_alt = @reaction_network begin
        @species X(t) Y(t) XY(t) Z(t) V(t) W(t)
        @parameters k1 k2 dX [diffusionparameter = true] dXY [diffusionparameter = true] dZ [
            diffusionparameter = true,
        ] dV [diffusionparameter = true] p1 p2
        (k1, k2), X + Y <--> XY
    end
    binding_srs_alt = [
        DiffusionReaction(:dX, :X),
        DiffusionReaction(:dXY, :XY),
        DiffusionReaction(:dZ, :Z),
        DiffusionReaction(:dV, :V),
    ]
    lrs_alt = LatticeReactionSystem(binding_system_alt, binding_srs_alt, small_2d_grid)
    u0_alt = [
        :X => 1.0,
        :Y => 2.0 * rand_v_vals(lrs_alt.lattice),
        :XY => 0.5,
        :Z => 2.0 * rand_v_vals(lrs_alt.lattice),
        :V => 0.5,
        :W => 1.0,
    ]
    p_alt = [
        :k1 => 2.0,
        :k2 => 0.1 .+ rand_v_vals(lrs_alt.lattice),
        :dX => 1.0 .+ rand_e_vals(lrs_alt.lattice),
        :dXY => 3.0,
        :dZ => rand_e_vals(lrs_alt.lattice),
        :dV => 0.2,
        :p1 => 1.0,
        :p2 => rand_v_vals(lrs_alt.lattice),
    ]
    oprob_alt = ODEProblem(lrs_alt, u0_alt, (0.0, 10.0), p_alt)
    ss_alt = solve(oprob_alt, Tsit5()).u[end]

    binding_srs_main = [DiffusionReaction(:dX, :X), DiffusionReaction(:dXY, :XY)]
    lrs = LatticeReactionSystem(binding_system, binding_srs_main, small_2d_grid)
    u0 = u0_alt[1:3]
    p = p_alt[1:4]
    oprob = ODEProblem(lrs, u0, (0.0, 10.0), p)
    ss = solve(oprob, Tsit5()).u[end]

    for i in 1:25
        @test isapprox(ss_alt[((i - 1) * 6 + 1):((i - 1) * 6 + 3)],
                       ss[((i - 1) * 3 + 1):((i - 1) * 3 + 3)]) < 1e-3
    end
end

# System with single spatial reaction.
let
    lrs_1 = LatticeReactionSystem(SIR_system, SIR_dif_S, small_2d_grid)
    lrs_2 = LatticeReactionSystem(SIR_system, [SIR_dif_S], small_2d_grid)
    u0 = [:S => 990.0, :I => 20.0 * rand_v_vals(lrs_1.lattice), :R => 0.0]
    pV = SIR_p
    pE = [:dS => 0.01]
    ss_1 = solve(ODEProblem(lrs_1, u0, (0.0, 500.0), (pV, pE)), Tsit5()).u[end]
    ss_2 = solve(ODEProblem(lrs_2, u0, (0.0, 500.0), (pV, pE)), Tsit5()).u[end]
    @test all(isequal.(ss_1, ss_2))
end

# Various ways to give parameters and initial conditions.
let
    lrs = LatticeReactionSystem(SIR_system, SIR_srs_2, very_small_2d_grid)
    u0_1 = [:S => 990.0, :I => [1.0, 3.0, 2.0, 5.0], :R => 0.0]
    u0_2 = [990.0, [1.0, 3.0, 2.0, 5.0], 0.0]
    u0_3 = [990.0 990.0 990.0 990.0; 1.0 3.0 2.0 5.0; 0.0 0.0 0.0 0.0]
    pV_1 = [:α => 0.1 / 1000, :β => [0.01, 0.02, 0.01, 0.03]]
    pV_2 = [0.1 / 1000, [0.01, 0.02, 0.01, 0.03]]
    pV_3 = [0.1/1000 0.1/1000 0.1/1000 0.1/1000; 0.01 0.02 0.01 0.03]
    pE_1 = [:dS => [0.01, 0.02, 0.03, 0.04], :dI => 0.01, :dR => 0.01]
    pE_2 = [[0.01, 0.02, 0.03, 0.04], :0.01, 0.01]
    pE_3 = [0.01 0.02 0.03 0.04; 0.01 0.01 0.01 0.01; 0.01 0.01 0.01 0.01]

    p1 = [
        :α => 0.1 / 1000,
        :β => [0.01, 0.02, 0.01, 0.03],
        :dS => [0.01, 0.02, 0.03, 0.04],
        :dI => 0.01,
        :dR => 0.01,
    ]
    ss_1_1 = solve(ODEProblem(lrs, u0_1, (0.0, 1.0), p1), Tsit5()).u[end]
    for u0 in [u0_1, u0_2, u0_3], pV in [pV_1, pV_2, pV_3], pE in [pE_1, pE_2, pE_3]
        ss = solve(ODEProblem(lrs, u0, (0.0, 1.0), (pV, pE)), Tsit5()).u[end]
        @test all(isequal.(ss, ss_1_1))
    end
end

# Checks that variosu combinations of jac and sparse gives the same result.
let
    lrs = LatticeReactionSystem(brusselator_system, brusselator_srs_1, small_2d_grid)
    u0 = [:X => rand_v_vals(lrs.lattice, 10), :Y => rand_v_vals(lrs.lattice, 10)]
    pV = brusselator_p
    pE = [:dX => 0.2]
    oprob = ODEProblem(lrs, u0, (0.0, 50.0), (pV, pE); jac = false, sparse = false)
    oprob_sparse = ODEProblem(lrs, u0, (0.0, 50.0), (pV, pE); jac = false, sparse = true)
    oprob_jac = ODEProblem(lrs, u0, (0.0, 50.0), (pV, pE); jac = true, sparse = false)
    oprob_sparse_jac = ODEProblem(lrs, u0, (0.0, 50.0), (pV, pE); jac = true, sparse = true)

    ss = solve(oprob, Rosenbrock23(); abstol = 1e-10, reltol = 1e-10).u[end]
    @test all(isapprox.(ss,
                        solve(oprob_sparse, Rosenbrock23(); abstol = 1e-10, reltol = 1e-10).u[end];
                        rtol = 0.0001))
    @test all(isapprox.(ss,
                        solve(oprob_jac, Rosenbrock23(); abstol = 1e-10, reltol = 1e-10).u[end];
                        rtol = 0.0001))
    @test all(isapprox.(ss,
                        solve(oprob_sparse_jac, Rosenbrock23(); abstol = 1e-10, reltol = 1e-10).u[end];
                        rtol = 0.0001))
end

# Splitting parameters by position
let
    lrs = LatticeReactionSystem(SIR_system, SIR_srs_2, small_2d_grid)
    u0 = [:S => 990.0, :I => 20.0 * rand_v_vals(lrs.lattice), :R => 0.0]
    p1 = ([:α => 0.1 / 1000, :β => 0.01], [:dS => 0.01, :dI => 0.01, :dR => 0.01])
    p2 = [:α => 0.1 / 1000, :β => 0.01, :dS => 0.01, :dI => 0.01, :dR => 0.01]
    oprob1 = ODEProblem(lrs, u0, (0.0, 500.0), p1; jac = false)
    oprob2 = ODEProblem(lrs, u0, (0.0, 500.0), p2; jac = false)

    @test all(isapprox.(solve(oprob1, Tsit5()).u[end], solve(oprob2, Tsit5()).u[end]))
end

### Compare to Hand-written Functions ###

# Compares the brusselator for a line of cells.
let
    function spatial_brusselator_f(du, u, p, t)
        # Non-spatial
        for i in 1:2:(length(u) - 1)
            du[i] = p[1] + 0.5 * (u[i]^2) * u[i + 1] - u[i] - p[2] * u[i]
            du[i + 1] = p[2] * u[i] - 0.5 * (u[i]^2) * u[i + 1]
        end

        # Spatial
        du[1] += p[3] * (u[3] - u[1])
        du[end - 1] += p[3] * (u[end - 3] - u[end - 1])
        for i in 3:2:(length(u) - 3)
            du[i] += p[3] * (u[i - 2] + u[i + 2] - 2u[i])
        end
    end
    function spatial_brusselator_jac(J, u, p, t)
        J .= 0
        # Non-spatial
        for i in 1:2:(length(u) - 1)
            J[i, i] = u[i] * u[i + 1] - 1 - p[2]
            J[i, i + 1] = 0.5 * (u[i]^2)
            J[i + 1, i] = p[2] - u[i] * u[i + 1]
            J[i + 1, i + 1] = -0.5 * (u[i]^2)
        end

        # Spatial
        J[1, 1] -= p[3]
        J[1, 3] += p[3]
        J[end - 1, end - 1] -= p[3]
        J[end - 1, end - 3] += p[3]
        for i in 3:2:(length(u) - 3)
            J[i, i] -= 2 * p[3]
            J[i, i - 2] += p[3]
            J[i, i + 2] += p[3]
        end
    end
    function spatial_brusselator_jac_sparse(J, u, p, t)
        # Spatial
        J.nzval .= 0.0
        J.nzval[7:6:(end - 9)] .= -2p[3]
        J.nzval[1] = -p[3]
        J.nzval[end - 3] = -p[3]
        J.nzval[3:3:(end - 4)] .= p[3]

        # Non-spatial
        for i in 1:1:Int64(lenth(u) / 2 - 1)
            j = 6(i - 1) + 1
            J.nzval[j] = u[i] * u[i + 1] - 1 - p[2]
            J.nzval[j + 1] = 0.5 * (u[i]^2)
            J.nzval[j + 3] = p[2] - u[i] * u[i + 1]
            J.nzval[j + 4] = -0.5 * (u[i]^2)
        end
        J.nzval[end - 3] = u[end - 1] * u[end] - 1 - p[end - 1]
        J.nzval[end - 2] = 0.5 * (u[end - 1]^2)
        J.nzval[end - 1] = p[2] - u[end - 1] * u[end]
        J.nzval[end] = -0.5 * (u[end - 1]^2)
    end
    function make_jac_prototype(u0)
        jac_prototype_pre = zeros(length(u0), length(u0))
        for i in 1:2:(length(u0) - 1)
            jac_prototype_pre[i, i] = 1
            jac_prototype_pre[i + 1, i] = 1
            jac_prototype_pre[i, i + 1] = 1
            jac_prototype_pre[i + 1, i + 1] = 1
        end
        for i in 3:2:(length(u0) - 1)
            jac_prototype_pre[i - 2, i] = 1
            jac_prototype_pre[i, i - 2] = 1
        end
        return sparse(jac_prototype_pre)
    end

    u0 = 2 * rand(10000)
    p = [1.0, 4.0, 0.1]
    tspan = (0.0, 100.0)

    ofun_hw_dense = ODEFunction(spatial_brusselator_f; jac = spatial_brusselator_jac)
    ofun_hw_sparse = ODEFunction(spatial_brusselator_f; jac = spatial_brusselator_jac,
                                 jac_prototype = make_jac_prototype(u0))

    lrs = LatticeReactionSystem(brusselator_system, brusselator_srs_1,
                                path_graph(Int64(length(u0) / 2)))
    u0V = [:X => u0[1:2:(end - 1)], :Y => u0[2:2:end]]
    pV = [:A => p[1], :B => p[2]]
    pE = [:dX => p[3]]
    ofun_aut_dense = ODEProblem(lrs, u0V, tspan, (pV, pE); jac = true, sparse = false).f
    ofun_aut_sparse = ODEProblem(lrs, u0V, tspan, (pV, pE); jac = true, sparse = true).f

    du_hw_dense = deepcopy(u0)
    du_hw_sparse = deepcopy(u0)
    du_aut_dense = deepcopy(u0)
    du_aut_sparse = deepcopy(u0)

    ofun_hw_dense(du_hw_dense, u0, p, 0.0)
    ofun_hw_sparse(du_hw_sparse, u0, p, 0.0)
    ofun_aut_dense(du_aut_dense, u0, p, 0.0)
    ofun_aut_sparse(du_aut_sparse, u0, p, 0.0)

    @test isapprox(du_hw_dense, du_aut_dense)
    @test isapprox(du_hw_sparse, du_aut_sparse)

    J_hw_dense = deepcopy(zeros(length(u0), length(u0)))
    J_hw_sparse = deepcopy(make_jac_prototype(u0))
    J_aut_dense = deepcopy(zeros(length(u0), length(u0)))
    J_aut_sparse = deepcopy(make_jac_prototype(u0))

    ofun_hw_dense.jac(J_hw_dense, u0, p, 0.0)
    ofun_hw_sparse.jac(J_hw_sparse, u0, p, 0.0)
    ofun_aut_dense.jac(J_aut_dense, u0, p, 0.0)
    ofun_aut_sparse.jac(J_aut_sparse, u0, p, 0.0)

    @test isapprox(J_hw_dense, J_aut_dense)
    @test isapprox(J_hw_sparse, J_aut_sparse)
end