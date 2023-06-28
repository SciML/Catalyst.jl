### Fetches Stuff ###

# Fetch packages.
using Catalyst, OrdinaryDiffEq, Random, Test
using BenchmarkTools, Statistics
using Graphs

# Sets rnd number.
using StableRNGs
rng = StableRNG(12345)


### Helper Functions ###
rand_v_vals(grid) = rand(length(vertices(grid)))
rand_v_vals(grid, x::Number) = rand_v_vals(grid)*x
rand_e_vals(grid) = rand(length(edges(grid)))
rand_e_vals(grid, x::Number) = rand_e_vals(grid)*x

function make_u0_matrix(value_map, vals, symbols)
    (length(symbols)==0) && (return zeros(0, length(vals)))
    d = Dict(value_map)
    return [(d[s] isa Vector) ? d[s][v] : d[s] for s in symbols, v in 1:length(vals)]
end


### Declares Models ###

# Small non-stiff network.
binding_system = @reaction_network begin
    (kB,kD), X +Y <--> XY
end
binding_p = [:kB => 2.0, :kD => 0.5]

binding_dif_x = DiffusionReaction(:dX, :X)
binding_dif_y = DiffusionReaction(:dY, :Y)
binding_dif_xy = DiffusionReaction(:dXY, :XY)
binding_osr_xy1 = OnewaySpatialReaction(:d_ord_1, [:X, :Y], [:X, :Y], [1, 1], [1, 1])
binding_osr_xy2 = OnewaySpatialReaction(:d_ord_2, [:X, :Y], [:XY], [1, 1], [1])
binding_osr_xy3 = OnewaySpatialReaction(:d_ord_3, [:X, :Y], [:X, :Y], [2, 2], [2, 2])
binding_sr_1 = SpatialReaction(:d_sr_1, ([:X], [:Y]), ([:Y], [:X]), ([1], [1]), ([1], [1]))
binding_sr_2 = SpatialReaction(:d_sr_2, ([:X, :Y], [:XY]), ([:XY], []), ([1, 1], [2]), ([1], Vector{Int64}()))
binding_srs_1 = [binding_dif_x]
binding_srs_2 = [binding_dif_x, binding_dif_y, binding_dif_xy]
binding_srs_3 = [binding_sr_1, binding_sr_2]
binding_srs_4 = [binding_dif_x, binding_dif_y, binding_dif_xy, binding_osr_xy1, binding_osr_xy2, binding_osr_xy3, binding_sr_1, binding_sr_2]

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
CuH_Amination_p = [:kp1 => 1.2, :kp2 => -0.72, :k1 => 0.57, :k_1 => -3.5, :k2 => -0.35, :k_2 => -0.77, :k3 => -0.025, :kam => -2.6, :kdc => -3.0]
CuH_Amination_u0 = [:CuoAc => 0.0065, :Ligand => 0.0072, :CuoAcLigand => 0.0, :Silane => 0.65, :CuHLigand => 0.0, :SilaneOAc => 0.0, :Styrene => 0.16, :AlkylCuLigand => 0.0, :Amine_E => 0.39, :AlkylAmine => 0.0, :Cu_ELigand => 0.0, :E_Silane => 0.0, :Amine => 0.0, :Decomposition => 0.0]

CuH_Amination_diff_1 = DiffusionReaction(:D1, :CuoAc)
CuH_Amination_diff_2 = DiffusionReaction(:D2, :Silane)
CuH_Amination_diff_3 = DiffusionReaction(:D3, :Cu_ELigand)
CuH_Amination_diff_4 = DiffusionReaction(:D4, :Amine)
CuH_Amination_diff_5 = DiffusionReaction(:D5, :CuHLigand)
CuH_Amination_srs_1 = [CuH_Amination_diff_1]
CuH_Amination_srs_2 = [CuH_Amination_diff_1, CuH_Amination_diff_2, CuH_Amination_diff_3, CuH_Amination_diff_4, CuH_Amination_diff_5]

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
binding_osr_x2 = OnewaySpatialReaction(:d_ord_1, [:X], [:X], [2], [2])
brusselator_dif_sr = SpatialReaction(:D, ([:X], [:Y]), ([:Y], [:X]), ([1], [2]), ([1], [1]))
brusselator_srs_1 = [brusselator_dif_x]
brusselator_srs_2 = [brusselator_dif_x, brusselator_dif_y]
brusselator_srs_3 = [binding_osr_x2]
brusselator_srs_4 = [brusselator_dif_x, brusselator_dif_sr]

# Mid-sized stiff system.
sigmaB_system = @reaction_network begin
    kDeg,       (w,w2,w2v,v,w2v2,vP,σB,w2σB) ⟶ ∅
    kDeg,       vPp ⟶ phos
    (kBw,kDw),  2w ⟷ w2
    (kB1,kD1),  w2 + v ⟷ w2v
    (kB2,kD2),  w2v + v ⟷ w2v2
    kK1,        w2v ⟶ w2 + vP
    kK2,        w2v2 ⟶ w2v + vP
    (kB3,kD3),  w2 + σB ⟷ w2σB
    (kB4,kD4),  w2σB + v ⟷ w2v + σB
    (kB5,kD5),  vP + phos ⟷ vPp
    kP,         vPp ⟶ v + phos
    v0*((1+F*σB)/(K+σB)),     ∅ ⟶ σB
    λW*v0*((1+F*σB)/(K+σB)),  ∅ ⟶ w
    λV*v0*((1+F*σB)/(K+σB)),  ∅ ⟶ v
end
sigmaB_p = [:kBw => 3600, :kDw => 18, :kB1 => 3600, :kB2 => 3600, :kB3 => 3600, :kB4 => 1800, :kB5 => 3600, 
            :kD1 => 18, :kD2 => 18, :kD3 => 18, :kD4 => 1800, :kD5 => 18, :kK1 => 36, :kK2 => 12, :kP => 180, :kDeg => 0.7, 
            :v0 => 0.4, :F => 30, :K => 0.2, :λW => 4, :λV => 4.5]
sigmaB_u0 = [:w => 1.0, :w2 => 1.0, :w2v => 1.0, :v => 1.0, :w2v2 => 1.0, :vP => 1.0, :σB => 1.0, :w2σB => 1.0, :vPp => 0.0, :phos => 0.4]   

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
    for srs in [Vector{DiffusionReaction}(), binding_srs_1, binding_srs_2]  
        lrs = LatticeReactionSystem(binding_system, srs, grid)
        u0_1 = [:X => 1.0, :Y => 2.0, :XY => 0.0]
        u0_2 = [:X => rand_v_vals(lrs.lattice), :Y => 2.0, :XY => 0.0]
        u0_3 = [:X => 1.0, :Y => rand_v_vals(lrs.lattice), :XY => rand_v_vals(lrs.lattice)]
        u0_4 = [:X => rand_v_vals(lrs.lattice), :Y => rand_v_vals(lrs.lattice), :XY => rand_v_vals(lrs.lattice,3)]
        u0_5 = make_u0_matrix(u0_3, vertices(lrs.lattice), map(s -> Symbol(s.f), species(lrs.rs)))
        for u0 in [u0_1, u0_2, u0_3, u0_4, u0_5]
            p1 = [:kB => 2.0, :kD => 0.5]
            p2 = [:kB => 2.0, :kD => rand_v_vals(lrs.lattice)]
            p3 = [:kB => rand_v_vals(lrs.lattice), :kD => rand_v_vals(lrs.lattice)]
            p4 = make_u0_matrix(p1, vertices(lrs.lattice), Symbol.(parameters(lrs.rs)))
            for pV in [p1, p2, p3, p4]                
                pE_1 = map(sp -> sp => 0.2, lrs.spatial_params)
                pE_2 = map(sp -> sp => rand(), lrs.spatial_params)
                pE_3 = map(sp -> sp => rand_e_vals(lrs.lattice, 0.2), lrs.spatial_params)
                pE_4 = make_u0_matrix(pE_3, edges(lrs.lattice), lrs.spatial_params)
                for pE in [pE_1, pE_2, pE_3, pE_4]
                    oprob = ODEProblem(lrs, u0, (0.0,10.0), (pV, pE))
                    @test SciMLBase.successful_retcode(solve(oprob, Tsit5()))

                    oprob = ODEProblem(lrs, u0, (0.0,10.0), (pV, pE); jac=false)
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
        u0_4 = make_u0_matrix(u0_3, vertices(lrs.lattice), map(s -> Symbol(s.f), species(lrs.rs)))
        for u0 in [u0_1, u0_2, u0_3, u0_4]
            p1 = [:A => 1.0, :B => 4.0]
            p2 = [:A => 0.5 .+ rand_v_vals(lrs.lattice, 0.5), :B => 4.0]
            p3 = [:A => 0.5 .+ rand_v_vals(lrs.lattice, 0.5), :B => 4.0 .+ rand_v_vals(lrs.lattice, 1.0)]
            p4 = make_u0_matrix(p2, vertices(lrs.lattice), Symbol.(parameters(lrs.rs)))
            for pV in [p1, p2, p3, p4]
                pE_1 = map(sp -> sp => 0.2, lrs.spatial_params)
                pE_2 = map(sp -> sp => rand(), lrs.spatial_params)
                pE_3 = map(sp -> sp => rand_e_vals(lrs.lattice, 0.2), lrs.spatial_params)
                pE_4 = make_u0_matrix(pE_3, edges(lrs.lattice), lrs.spatial_params)
                for pE in [pE_1, pE_2, pE_3, pE_4]
                    oprob = ODEProblem(lrs, u0, (0.0,10.0), (pV, pE))
                    @test SciMLBase.successful_retcode(solve(oprob, QNDF()))

                    oprob = ODEProblem(lrs, u0, (0.0,10.0), (pV, pE); sparse=false)
                    @test SciMLBase.successful_retcode(solve(oprob, QNDF()))
                end
            end
        end
    end
end


### Tests Runtimes ###

# Timings currently are from Torkel's computer.


# Small grid, small, non-stiff, system.
let 
    lrs = LatticeReactionSystem(binding_system, binding_srs_2, small_2d_grid)
    u0 = [:X => rand_v_vals(lrs.lattice), :Y => rand_v_vals(lrs.lattice), :XY => rand_v_vals(lrs.lattice)]
    pV = binding_p
    pE = [:dX => 0.1, :dY => 0.2]
    oprob = ODEProblem(lrs, u0, (0.0,10.0), (pV,pE); jac=false)
    @test SciMLBase.successful_retcode(solve(oprob, Tsit5()))
    
    runtime_target = 0.001
    runtime = minimum((@benchmark solve($oprob, Tsit5())).times)/1000000000
    println("Small grid, small, non-stiff, system. Runtime: $(runtime), previous standard: $(runtime_target)")
    @test runtime < 1.2*runtime_target
end


# Large grid, small, non-stiff, system.
let 
    lrs = LatticeReactionSystem(binding_system, binding_srs_2, large_2d_grid)
    u0 = [:X => rand_v_vals(lrs.lattice), :Y => rand_v_vals(lrs.lattice), :XY => rand_v_vals(lrs.lattice)]
    pV = binding_p
    pE = [:dX => 0.1, :dY => 0.2]
    oprob = ODEProblem(lrs, u0, (0.0,10.0), (pV,pE); jac=false)
    @test SciMLBase.successful_retcode(solve(oprob, Tsit5()))
    
    runtime_target = 0.5
    runtime = minimum((@benchmark solve($oprob, Tsit5())).times)/1000000000
    println("Large grid, small, non-stiff, system. Runtime: $(runtime), previous standard: $(runtime_target)")
    @test runtime < 1.2*runtime_target
end

# Small grid, small, stiff, system.
let 
    lrs = LatticeReactionSystem(brusselator_system, brusselator_srs_1, small_2d_grid)
    u0 = [:X => rand_v_vals(lrs.lattice,10), :Y => rand_v_vals(lrs.lattice,10)]
    pV = brusselator_p
    pE = [:dX => 0.2,]
    oprob = ODEProblem(lrs, u0, (0.0,100.0), (pV,pE))
    @test SciMLBase.successful_retcode(solve(oprob, QNDF()))
    
    runtime_target = 0.1
    runtime = minimum((@benchmark solve($oprob, QNDF())).times)/1000000000
    println("Small grid, small, stiff, system. Runtime: $(runtime), previous standard: $(runtime_target)")
    @test runtime < 1.2*runtime_target
end

# Large grid, small, stiff, system.
let 
    lrs = LatticeReactionSystem(brusselator_system, brusselator_srs_1, large_2d_grid)
    u0 = [:X => rand_v_vals(lrs.lattice,10), :Y => rand_v_vals(lrs.lattice,10)]
    pV = brusselator_p
    pE = [:dX => 0.2,]
    oprob = ODEProblem(lrs, u0, (0.0,100.0), (pV,pE))
    @test SciMLBase.successful_retcode(solve(oprob, QNDF()))
    
    runtime_target = 200.0
    runtime = minimum((@benchmark solve($oprob, QNDF())).times)/1000000000
    println("Large grid, small, stiff, system. Runtime: $(runtime), previous standard: $(runtime_target)")
    @test runtime < 1.2*runtime_target
end


# Small grid, mid-sized, non-stiff, system.
let 
    lrs = LatticeReactionSystem(CuH_Amination_system, CuH_Amination_srs_2, small_2d_grid)
    u0 = [:CuoAc => 0.005 .+ rand_v_vals(lrs.lattice, 0.005), :Ligand => 0.005 .+ rand_v_vals(lrs.lattice, 0.005), :CuoAcLigand => 0.0, :Silane => 0.5 .+ rand_v_vals(lrs.lattice, 0.5), :CuHLigand => 0.0, :SilaneOAc => 0.0, :Styrene => 0.16, :AlkylCuLigand => 0.0, :Amine_E => 0.39, :AlkylAmine => 0.0, :Cu_ELigand => 0.0, :E_Silane => 0.0, :Amine => 0.0, :Decomposition => 0.0]
    pV = CuH_Amination_p
    pE = [:D1 => 0.1, :D2 => 0.1, :D3 => 0.1, :D4 => 0.1, :D5 => 0.1]
    oprob = ODEProblem(lrs, u0, (0.0,10.0), (pV,pE); jac=false)
    @test SciMLBase.successful_retcode(solve(oprob, Tsit5()))
    
    runtime_target = 0.005
    runtime = minimum((@benchmark solve($oprob, Tsit5())).times)/1000000000
    println("Small grid, mid-sized, non-stiff, system. Runtime: $(runtime), previous standard: $(runtime_target)")
    @test runtime < 1.2*runtime_target
end

# Large grid, mid-sized, non-stiff, system.
let 
    lrs = LatticeReactionSystem(CuH_Amination_system, CuH_Amination_srs_2, large_2d_grid)
    u0 = [:CuoAc => 0.005 .+ rand_v_vals(lrs.lattice, 0.005), :Ligand => 0.005 .+ rand_v_vals(lrs.lattice, 0.005), :CuoAcLigand => 0.0, :Silane => 0.5 .+ rand_v_vals(lrs.lattice, 0.5), :CuHLigand => 0.0, :SilaneOAc => 0.0, :Styrene => 0.16, :AlkylCuLigand => 0.0, :Amine_E => 0.39, :AlkylAmine => 0.0, :Cu_ELigand => 0.0, :E_Silane => 0.0, :Amine => 0.0, :Decomposition => 0.0]
    pV = CuH_Amination_p
    pE = [:D1 => 0.1, :D2 => 0.1, :D3 => 0.1, :D4 => 0.1, :D5 => 0.1]
    oprob = ODEProblem(lrs, u0, (0.0,10.0), (pV,pE); jac=false)
    @test SciMLBase.successful_retcode(solve(oprob, Tsit5()))
    
    runtime_target = 2
    runtime = minimum((@benchmark solve($oprob, Tsit5())).times)/1000000000
    println("Large grid, mid-sized, non-stiff, system. Runtime: $(runtime), previous standard: $(runtime_target)")
    @test runtime < 1.2*runtime_target
end

# Small grid, mid-sized, stiff, system.
let 
    lrs = LatticeReactionSystem(sigmaB_system, sigmaB_srs_2, small_2d_grid)
    u0 = [:w => 0.5 .+ rand_v_vals(lrs.lattice, 0.5), :w2 => 0.5 .+ rand_v_vals(lrs.lattice, 0.5), :w2v => 0.5 .+ rand_v_vals(lrs.lattice, 0.5), :v => 0.5 .+ rand_v_vals(lrs.lattice, 0.5), :w2v2 => 0.5 .+ rand_v_vals(lrs.lattice, 0.5), :vP => 0.5 .+ rand_v_vals(lrs.lattice, 0.5), :σB => 0.5 .+ rand_v_vals(lrs.lattice, 0.5), :w2σB => 0.5 .+ rand_v_vals(lrs.lattice, 0.5), :vPp => 0.0, :phos => 0.4]
    pV = sigmaB_p
    pE = [:DσB => 0.1, :Dw => 0.1, :Dv => 0.1]
    oprob = ODEProblem(lrs, u0, (0.0,10.0), (pV,pE))
    @test SciMLBase.successful_retcode(solve(oprob, QNDF()))
    
    runtime_target = 0.025
    runtime = minimum((@benchmark solve($oprob, QNDF())).times)/1000000000
    println("Small grid, mid-sized, non-stiff, system. Runtime: $(runtime), previous standard: $(runtime_target)")
    @test runtime < 1.2*runtime_target
end

# Large grid, mid-sized, stiff, system.
let 
    lrs = LatticeReactionSystem(sigmaB_system, sigmaB_srs_2, large_2d_grid)
    u0 = [:w => 0.5 .+ rand_v_vals(lrs.lattice, 0.5), :w2 => 0.5 .+ rand_v_vals(lrs.lattice, 0.5), :w2v => 0.5 .+ rand_v_vals(lrs.lattice, 0.5), :v => 0.5 .+ rand_v_vals(lrs.lattice, 0.5), :w2v2 => 0.5 .+ rand_v_vals(lrs.lattice, 0.5), :vP => 0.5 .+ rand_v_vals(lrs.lattice, 0.5), :σB => 0.5 .+ rand_v_vals(lrs.lattice, 0.5), :w2σB => 0.5 .+ rand_v_vals(lrs.lattice, 0.5), :vPp => 0.0, :phos => 0.4]
    pV = sigmaB_p
    pE = [:DσB => 0.1, :Dw => 0.1, :Dv => 0.1]
    oprob = ODEProblem(lrs, u0, (0.0,10.0), (pV,pE))
    @test SciMLBase.successful_retcode(solve(oprob, QNDF()))
    
    runtime_target = 150.
    runtime = minimum((@benchmark solve($oprob, QNDF())).times)/1000000000
    println("Large grid, mid-sized, non-stiff, system. Runtime: $(runtime), previous standard: $(runtime_target)")
    @test runtime < 1.2*runtime_target
end