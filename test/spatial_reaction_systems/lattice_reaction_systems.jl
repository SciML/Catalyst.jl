### Fetches Stuff ###

# Fetch packages.
using Catalyst, OrdinaryDiffEq, Random, Test
using BenchmarkTools, LinearAlgebra, Statistics
using Graphs

# Sets rnd number.
using StableRNGs
rng = StableRNG(12345)

### Helper Functions ###
rand_v_vals(grid) = rand(length(vertices(grid)))
rand_v_vals(grid, x::Number) = rand_v_vals(grid) * x
rand_e_vals(grid) = rand(length(edges(grid)))
rand_e_vals(grid, x::Number) = rand_e_vals(grid) * x

function make_u0_matrix(value_map, vals, symbols)
    (length(symbols) == 0) && (return zeros(0, length(vals)))
    d = Dict(value_map)
    return [(d[s] isa Vector) ? d[s][v] : d[s] for s in symbols, v in 1:length(vals)]
end

### Declares Models ###

# Small non-stiff network.
SIR_system = @reaction_network begin
    α, S + I --> 2I
    β, I --> R
end
SIR_p = [:α => 0.1 / 1000, :β => 0.01]
SIR_u0 = [:S => 999.0, :I => 1.0, :R => 0.0]

SIR_dif_S = DiffusionReaction(:dS, :S)
SIR_dif_I = DiffusionReaction(:dI, :I)
SIR_dif_R = DiffusionReaction(:dR, :R)
SIR_srs_1 = [SIR_dif_S]
SIR_srs_2 = [SIR_dif_S, SIR_dif_I, SIR_dif_R]

# Mid-sized non-stiff system.
CuH_Amination_system = @reaction_network begin
    10.0^kp1, CuoAc + Ligand --> CuoAcLigand
    10.0^kp2, CuoAcLigand + Silane --> CuHLigand + SilaneOAc
    10.0^k1, CuHLigand + Styrene --> AlkylCuLigand
    10.0^k_1, AlkylCuLigand --> CuHLigand + Styrene
    10.0^k2, AlkylCuLigand + Amine_E --> AlkylAmine + Cu_ELigand
    10.0^k_2, AlkylAmine + Cu_ELigand --> AlkylCuLigand + Amine_E
    10.0^k3, Cu_ELigand + Silane --> CuHLigand + E_Silane
    10.0^kam, CuHLigand + Amine_E --> Amine + Cu_ELigand
    10.0^kdc, CuHLigand + CuHLigand --> Decomposition
end
CuH_Amination_p = [
    :kp1 => 1.2,
    :kp2 => -0.72,
    :k1 => 0.57,
    :k_1 => -3.5,
    :k2 => -0.35,
    :k_2 => -0.77,
    :k3 => -0.025,
    :kam => -2.6,
    :kdc => -3.0,
]
CuH_Amination_u0 = [
    :CuoAc => 0.0065,
    :Ligand => 0.0072,
    :CuoAcLigand => 0.0,
    :Silane => 0.65,
    :CuHLigand => 0.0,
    :SilaneOAc => 0.0,
    :Styrene => 0.16,
    :AlkylCuLigand => 0.0,
    :Amine_E => 0.39,
    :AlkylAmine => 0.0,
    :Cu_ELigand => 0.0,
    :E_Silane => 0.0,
    :Amine => 0.0,
    :Decomposition => 0.0,
]

CuH_Amination_diff_1 = DiffusionReaction(:D1, :CuoAc)
CuH_Amination_diff_2 = DiffusionReaction(:D2, :Silane)
CuH_Amination_diff_3 = DiffusionReaction(:D3, :Cu_ELigand)
CuH_Amination_diff_4 = DiffusionReaction(:D4, :Amine)
CuH_Amination_diff_5 = DiffusionReaction(:D5, :CuHLigand)
CuH_Amination_srs_1 = [CuH_Amination_diff_1]
CuH_Amination_srs_2 = [
    CuH_Amination_diff_1,
    CuH_Amination_diff_2,
    CuH_Amination_diff_3,
    CuH_Amination_diff_4,
    CuH_Amination_diff_5,
]

# Small stiff system.
brusselator_system = @reaction_network begin
    A, ∅ → X
    1, 2X + Y → 3X
    B, X → Y
    1, X → ∅
end
brusselator_p = [:A => 1.0, :B => 4.0]

brusselator_dif_x = DiffusionReaction(:dX, :X)
brusselator_dif_y = DiffusionReaction(:dY, :Y)
brusselator_srs_1 = [brusselator_dif_x]
brusselator_srs_2 = [brusselator_dif_x, brusselator_dif_y]

# Mid-sized stiff system.
sigmaB_system = @reaction_network begin
    kDeg, (w, w2, w2v, v, w2v2, vP, σB, w2σB) ⟶ ∅
    kDeg, vPp ⟶ phos
    (kBw, kDw), 2w ⟷ w2
    (kB1, kD1), w2 + v ⟷ w2v
    (kB2, kD2), w2v + v ⟷ w2v2
    kK1, w2v ⟶ w2 + vP
    kK2, w2v2 ⟶ w2v + vP
    (kB3, kD3), w2 + σB ⟷ w2σB
    (kB4, kD4), w2σB + v ⟷ w2v + σB
    (kB5, kD5), vP + phos ⟷ vPp
    kP, vPp ⟶ v + phos
    v0 * ((1 + F * σB) / (K + σB)), ∅ ⟶ σB
    λW * v0 * ((1 + F * σB) / (K + σB)), ∅ ⟶ w
    λV * v0 * ((1 + F * σB) / (K + σB)), ∅ ⟶ v
end
sigmaB_p = [:kBw => 3600, :kDw => 18, :kB1 => 3600, :kB2 => 3600, :kB3 => 3600,
    :kB4 => 1800, :kB5 => 3600,
    :kD1 => 18, :kD2 => 18, :kD3 => 18, :kD4 => 1800, :kD5 => 18, :kK1 => 36, :kK2 => 12,
    :kP => 180, :kDeg => 0.7,
    :v0 => 0.4, :F => 30, :K => 0.2, :λW => 4, :λV => 4.5]
sigmaB_u0 = [
    :w => 1.0,
    :w2 => 1.0,
    :w2v => 1.0,
    :v => 1.0,
    :w2v2 => 1.0,
    :vP => 1.0,
    :σB => 1.0,
    :w2σB => 1.0,
    :vPp => 0.0,
    :phos => 0.4,
]

sigmaB_dif_σB = DiffusionReaction(:DσB, :σB)
sigmaB_dif_w = DiffusionReaction(:Dw, :w)
sigmaB_dif_v = DiffusionReaction(:Dv, :v)
sigmaB_srs_1 = [sigmaB_dif_σB]
sigmaB_srs_2 = [sigmaB_dif_σB, sigmaB_dif_w, sigmaB_dif_v]

### Declares Lattices ###

# Grids.
small_2d_grid = Graphs.grid([5, 5])
medium_2d_grid = Graphs.grid([20, 20])
large_2d_grid = Graphs.grid([100, 100])

small_3d_grid = Graphs.grid([5, 5, 5])
medium_3d_grid = Graphs.grid([20, 20, 20])
large_3d_grid = Graphs.grid([100, 100, 100])

# Paths.
short_path = path_graph(100)
long_path = path_graph(1000)

# Directed cycle.
small_directed_cycle = cycle_graph(100)
large_directed_cycle = cycle_graph(1000)

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
                pE_1 = map(sp -> sp => 0.01, lrs.spatial_params)
                pE_2 = map(sp -> sp => 0.01, lrs.spatial_params)
                pE_3 = map(sp -> sp => rand_e_vals(lrs.lattice, 0.01), lrs.spatial_params)
                pE_4 = make_u0_matrix(pE_3, edges(lrs.lattice), lrs.spatial_params)
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
                pE_1 = map(sp -> sp => 0.2, lrs.spatial_params)
                pE_2 = map(sp -> sp => rand(), lrs.spatial_params)
                pE_3 = map(sp -> sp => rand_e_vals(lrs.lattice, 0.2), lrs.spatial_params)
                pE_4 = make_u0_matrix(pE_3, edges(lrs.lattice), lrs.spatial_params)
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


### Tests Special Cases ###

# Create network with vaious combinations of graph/di-graph and parameters.
lrs_digraph = LatticeReactionSystem(SIR_system, SIR_srs_2, complete_digraph(3))
lrs_graph = LatticeReactionSystem(SIR_system, SIR_srs_2, complete_graph(3))
u0 = [:S => 990.0, :I => 20.0*rand_v_vals(lrs_digraph.lattice), :R => 0.0]
pV = SIR_p
pE_digraph_1 = [:dS => [0.10, 0.10, 0.12, 0.12, 0.14, 0.14], :dI => 0.01, :dR => 0.01]
pE_digraph_2 = [[0.10, 0.10, 0.12, 0.12, 0.14, 0.14], 0.01, 0.01]
pE_digraph_3 = [0.10 0.10 0.12 0.12 0.14 0.14; 0.01 0.01 0.01 0.01 0.01 0.01; 0.01 0.01 0.01 0.01 0.01 0.01]
pE_graph_1 = [:dS => [0.10, 0.12, 0.14], :dI => 0.01, :dR => 0.01]
pE_graph_2 = [[0.10, 0.12, 0.14], 0.01, 0.01]
pE_graph_3 = [0.10 0.12 0.14 ; 0.01 0.01 0.01; 0.01 0.01 0.01]
oprob_digraph_1 = ODEProblem(lrs_digraph, u0, (0.0,500.0), (pV,pE_digraph_1))
oprob_digraph_2 = ODEProblem(lrs_digraph, u0, (0.0,500.0), (pV,pE_digraph_2))
oprob_digraph_3 = ODEProblem(lrs_digraph, u0, (0.0,500.0), (pV,pE_digraph_3))
oprob_graph_11 = ODEProblem(lrs_graph, u0, (0.0,500.0), (pV,pE_digraph_1))
oprob_graph_12 = ODEProblem(lrs_graph, u0, (0.0,500.0), (pV,pE_graph_1))
oprob_graph_21 = ODEProblem(lrs_graph, u0, (0.0,500.0), (pV,pE_digraph_2))
oprob_graph_22 = ODEProblem(lrs_graph, u0, (0.0,500.0), (pV,pE_graph_2))
oprob_graph_31 = ODEProblem(lrs_graph, u0, (0.0,500.0), (pV,pE_digraph_3))
oprob_graph_32 = ODEProblem(lrs_graph, u0, (0.0,500.0), (pV,pE_graph_3))
sim_end_digraph_1 = solve(oprob_digraph_1, Tsit5()).u[end]
sim_end_digraph_2 = solve(oprob_digraph_2, Tsit5()).u[end]
sim_end_digraph_3 = solve(oprob_digraph_3, Tsit5()).u[end]
sim_end_graph_11 = solve(oprob_graph_11, Tsit5()).u[end]
sim_end_graph_12 = solve(oprob_graph_12, Tsit5()).u[end]
sim_end_graph_21 = solve(oprob_graph_21, Tsit5()).u[end]
sim_end_graph_22 = solve(oprob_graph_22, Tsit5()).u[end]
sim_end_graph_31 = solve(oprob_graph_31, Tsit5()).u[end]
sim_end_graph_32 = solve(oprob_graph_32, Tsit5()).u[end]

@test all(sim_end_digraph_1 .== sim_end_digraph_2 .== sim_end_digraph_3 .== sim_end_graph_11 .== sim_end_graph_12 .== sim_end_graph_21 .== sim_end_graph_22 .== sim_end_graph_31 .== sim_end_graph_32)

# Creates networks with empty species or parameters.
binding_system_1 = @reaction_network begin
    @species X(t) Y(t) XY(t) Z(t) V(t) W(t)
    @parameters k1 k2 dX dXY dZ dV p1 p2
    (k1, k2), X + Y <--> XY
end
binding_sr_1 = [DiffusionReaction(:dX, :X), DiffusionReaction(:dXY, :XY), DiffusionReaction(:dZ, :Z), DiffusionReaction(:dV, :V)]
binding_lrs_1 = LatticeReactionSystem(binding_system_1, binding_sr_1, Graphs.grid([5, 5]))
binding_u0_1 = [:X => 1.0, :Y => 2.0*rand_v_vals(binding_lrs_1.lattice), :XY => 0.5, :Z => 2.0*rand_v_vals(binding_lrs_1.lattice), :V => 0.5, :W => 1.0]
binding_p_1 = [:k1 => 2.0, :k2 => 0.1 .+ rand_v_vals(binding_lrs_1.lattice), :dX => 1.0 .+ rand_e_vals(binding_lrs_1.lattice), :dXY => 3.0, :dZ => rand_e_vals(binding_lrs_1.lattice), :dV => 0.2, :p1 => 1.0, :p2 => rand_v_vals(binding_lrs_1.lattice)]
binding_oprob_1 = ODEProblem(binding_lrs_1, binding_u0_1, (0.0,10.0), binding_p_1)
binding_sol_1 = solve(binding_oprob_1, Tsit5())

binding_system_2 = @reaction_network begin
    (k1, k2), X + Y <--> XY
end
binding_sr_2 = [DiffusionReaction(:dX, :X), DiffusionReaction(:dXY, :XY)]
binding_lrs_2 = LatticeReactionSystem(binding_system_2, binding_sr_2, Graphs.grid([5, 5]))
binding_u0_2 = binding_u0_1[1:3]
binding_p_2 = binding_p_1[1:4]
binding_oprob_2 = ODEProblem(binding_lrs_2, binding_u0_2, (0.0,10.0), binding_p_2)
binding_sol_2 = solve(binding_oprob_2, Tsit5())

for i = 1:25
    @test norm(binding_sol_1.u[end][(i-1)*6+1:(i-1)*6+3] .- binding_sol_2.u[end][(i-1)*3+1:(i-1)*3+3]) <1e-3
end


### Tests Runtimes ###
# Timings currently are from Torkel's computer.

# Small grid, small, non-stiff, system.
let
    lrs = LatticeReactionSystem(SIR_system, SIR_srs_2, small_2d_grid)
    u0 = [:S => 990.0, :I => 20.0 * rand_v_vals(lrs.lattice), :R => 0.0]
    pV = SIR_p
    pE = [:dS => 0.01, :dI => 0.01, :dR => 0.01]
    oprob = ODEProblem(lrs, u0, (0.0, 500.0), (pV, pE); jac = false)
    @test SciMLBase.successful_retcode(solve(oprob, Tsit5()))

    runtime_target = 0.00023
    runtime = minimum((@benchmark solve($oprob, Tsit5())).times) / 1000000000
    println("Small grid, small, non-stiff, system. Runtime: $(runtime), previous standard: $(runtime_target)")
    @test runtime < 1.2 * runtime_target
end

# Large grid, small, non-stiff, system.
let
    lrs = LatticeReactionSystem(SIR_system, SIR_srs_2, large_2d_grid)
    u0 = [:S => 990.0, :I => 20.0 * rand_v_vals(lrs.lattice), :R => 0.0]
    pV = SIR_p
    pE = [:dS => 0.01, :dI => 0.01, :dR => 0.01]
    oprob = ODEProblem(lrs, u0, (0.0, 500.0), (pV, pE); jac = false)
    @test SciMLBase.successful_retcode(solve(oprob, Tsit5()))

    runtime_target = 0.1
    runtime = minimum((@benchmark solve($oprob, Tsit5())).times) / 1000000000
    println("Large grid, small, non-stiff, system. Runtime: $(runtime), previous standard: $(runtime_target)")
    @test runtime < 1.2 * runtime_target
end

# Small grid, small, stiff, system.
let
    lrs = LatticeReactionSystem(brusselator_system, brusselator_srs_1, small_2d_grid)
    u0 = [:X => rand_v_vals(lrs.lattice, 10), :Y => rand_v_vals(lrs.lattice, 10)]
    pV = brusselator_p
    pE = [:dX => 0.2]
    oprob = ODEProblem(lrs, u0, (0.0, 100.0), (pV, pE))
    @test SciMLBase.successful_retcode(solve(oprob, QNDF()))

    runtime_target = 0.05
    runtime = minimum((@benchmark solve($oprob, QNDF())).times) / 1000000000
    println("Small grid, small, stiff, system. Runtime: $(runtime), previous standard: $(runtime_target)")
    @test runtime < 10 * runtime_target
end

# Large grid, small, stiff, system.
let
    lrs = LatticeReactionSystem(brusselator_system, brusselator_srs_1, large_2d_grid)
    u0 = [:X => rand_v_vals(lrs.lattice, 10), :Y => rand_v_vals(lrs.lattice, 10)]
    pV = brusselator_p
    pE = [:dX => 0.2]
    oprob = ODEProblem(lrs, u0, (0.0, 100.0), (pV, pE))
    @test SciMLBase.successful_retcode(solve(oprob, QNDF()))

    runtime_target = 140.0
    runtime = minimum((@benchmark solve($oprob, QNDF())).times) / 1000000000
    println("Large grid, small, stiff, system. Runtime: $(runtime), previous standard: $(runtime_target)")
    @test runtime < 10 * runtime_target
end

# Small grid, mid-sized, non-stiff, system.
let
    lrs = LatticeReactionSystem(CuH_Amination_system, CuH_Amination_srs_2, small_2d_grid)
    u0 = [
        :CuoAc => 0.005 .+ rand_v_vals(lrs.lattice, 0.005),
        :Ligand => 0.005 .+ rand_v_vals(lrs.lattice, 0.005),
        :CuoAcLigand => 0.0,
        :Silane => 0.5 .+ rand_v_vals(lrs.lattice, 0.5),
        :CuHLigand => 0.0,
        :SilaneOAc => 0.0,
        :Styrene => 0.16,
        :AlkylCuLigand => 0.0,
        :Amine_E => 0.39,
        :AlkylAmine => 0.0,
        :Cu_ELigand => 0.0,
        :E_Silane => 0.0,
        :Amine => 0.0,
        :Decomposition => 0.0,
    ]
    pV = CuH_Amination_p
    pE = [:D1 => 0.1, :D2 => 0.1, :D3 => 0.1, :D4 => 0.1, :D5 => 0.1]
    oprob = ODEProblem(lrs, u0, (0.0, 10.0), (pV, pE); jac = false)
    @test SciMLBase.successful_retcode(solve(oprob, Tsit5()))

    runtime_target = 0.00293
    runtime = minimum((@benchmark solve($oprob, Tsit5())).times) / 1000000000
    println("Small grid, mid-sized, non-stiff, system. Runtime: $(runtime), previous standard: $(runtime_target)")
    @test runtime < 10 * runtime_target
end

# Large grid, mid-sized, non-stiff, system.
let
    lrs = LatticeReactionSystem(CuH_Amination_system, CuH_Amination_srs_2, large_2d_grid)
    u0 = [
        :CuoAc => 0.005 .+ rand_v_vals(lrs.lattice, 0.005),
        :Ligand => 0.005 .+ rand_v_vals(lrs.lattice, 0.005),
        :CuoAcLigand => 0.0,
        :Silane => 0.5 .+ rand_v_vals(lrs.lattice, 0.5),
        :CuHLigand => 0.0,
        :SilaneOAc => 0.0,
        :Styrene => 0.16,
        :AlkylCuLigand => 0.0,
        :Amine_E => 0.39,
        :AlkylAmine => 0.0,
        :Cu_ELigand => 0.0,
        :E_Silane => 0.0,
        :Amine => 0.0,
        :Decomposition => 0.0,
    ]
    pV = CuH_Amination_p
    pE = [:D1 => 0.1, :D2 => 0.1, :D3 => 0.1, :D4 => 0.1, :D5 => 0.1]
    oprob = ODEProblem(lrs, u0, (0.0, 10.0), (pV, pE); jac = false)
    @test SciMLBase.successful_retcode(solve(oprob, Tsit5()))

    runtime_target = 1.257
    runtime = minimum((@benchmark solve($oprob, Tsit5())).times) / 1000000000
    println("Large grid, mid-sized, non-stiff, system. Runtime: $(runtime), previous standard: $(runtime_target)")
    @test runtime < 10 * runtime_target
end

# Small grid, mid-sized, stiff, system.
let
    lrs = LatticeReactionSystem(sigmaB_system, sigmaB_srs_2, small_2d_grid)
    u0 = [
        :w => 0.5 .+ rand_v_vals(lrs.lattice, 0.5),
        :w2 => 0.5 .+ rand_v_vals(lrs.lattice, 0.5),
        :w2v => 0.5 .+ rand_v_vals(lrs.lattice, 0.5),
        :v => 0.5 .+ rand_v_vals(lrs.lattice, 0.5),
        :w2v2 => 0.5 .+ rand_v_vals(lrs.lattice, 0.5),
        :vP => 0.5 .+ rand_v_vals(lrs.lattice, 0.5),
        :σB => 0.5 .+ rand_v_vals(lrs.lattice, 0.5),
        :w2σB => 0.5 .+ rand_v_vals(lrs.lattice, 0.5),
        :vPp => 0.0,
        :phos => 0.4,
    ]
    pV = sigmaB_p
    pE = [:DσB => 0.1, :Dw => 0.1, :Dv => 0.1]
    oprob = ODEProblem(lrs, u0, (0.0, 10.0), (pV, pE))
    @test SciMLBase.successful_retcode(solve(oprob, QNDF()))

    runtime_target = 0.023
    runtime = minimum((@benchmark solve($oprob, QNDF())).times) / 1000000000
    println("Small grid, mid-sized, stiff, system. Runtime: $(runtime), previous standard: $(runtime_target)")
    @test runtime < 10 * runtime_target
end

# Large grid, mid-sized, stiff, system.
let
    lrs = LatticeReactionSystem(sigmaB_system, sigmaB_srs_2, large_2d_grid)
    u0 = [
        :w => 0.5 .+ rand_v_vals(lrs.lattice, 0.5),
        :w2 => 0.5 .+ rand_v_vals(lrs.lattice, 0.5),
        :w2v => 0.5 .+ rand_v_vals(lrs.lattice, 0.5),
        :v => 0.5 .+ rand_v_vals(lrs.lattice, 0.5),
        :w2v2 => 0.5 .+ rand_v_vals(lrs.lattice, 0.5),
        :vP => 0.5 .+ rand_v_vals(lrs.lattice, 0.5),
        :σB => 0.5 .+ rand_v_vals(lrs.lattice, 0.5),
        :w2σB => 0.5 .+ rand_v_vals(lrs.lattice, 0.5),
        :vPp => 0.0,
        :phos => 0.4,
    ]
    pV = sigmaB_p
    pE = [:DσB => 0.1, :Dw => 0.1, :Dv => 0.1]
    oprob = ODEProblem(lrs, u0, (0.0, 10.0), (pV, pE))
    @test SciMLBase.successful_retcode(solve(oprob, QNDF()))

    runtime_target = 111.0
    runtime = minimum((@benchmark solve($oprob, QNDF())).times) / 1000000000
    println("Large grid, mid-sized, stiff, system. Runtime: $(runtime), previous standard: $(runtime_target)")
    @test runtime < 10 * runtime_target
end
