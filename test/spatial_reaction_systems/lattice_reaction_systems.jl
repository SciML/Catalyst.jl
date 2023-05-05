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
    # Stiff case
    for srs in [Vector{SpatialReaction}(), brusselator_srs_1, brusselator_srs_2, brusselator_srs_3, brusselator_srs_4] 
        lrs = LatticeReactionSystem(brusselator_system, srs, grid)
        u0_1 = [:X => 1.0, :Y => 20.0]
        u0_2 = [:X => rand_v_vals(grid, 10.0), :Y => 2.0]
        u0_3 = [:X => rand_v_vals(grid, 20), :Y => rand_v_vals(grid, 10)]
        u0_4 = make_u0_matrix(u0_3, vertices(grid), map(s -> Symbol(s.f), species(lrs.rs)))
        for u0 in [u0_1, u0_2, u0_3, u0_4]
            p1 = [:A => 1.0, :B => 4.0]
            p2 = [:A => 0.5 .+ rand_v_vals(grid, 0.5), :B => 4.0]
            p3 = [:A => 0.5 .+ rand_v_vals(grid, 0.5), :B => 4.0 .+ rand_v_vals(grid, 1.0)]
            p4 = make_u0_matrix(p2, vertices(grid), Symbol.(parameters(lrs.rs)))
            for pV in [p1, p2, p3, p4]
                pE_1 = map(sp -> sp => 0.2, lrs.spatial_params)
                pE_2 = map(sp -> sp => rand(), lrs.spatial_params)
                pE_3 = map(sp -> sp => rand_e_vals(grid, 0.2), lrs.spatial_params)
                pE_4 = make_u0_matrix(pE_3, edges(grid), lrs.spatial_params)
                for pE in [pE_1, pE_2, pE_3, pE_4]
                    oprob = ODEProblem(lrs, u0, (0.0,10.0), (pV, pE))
                    @test SciMLBase.successful_retcode(solve(oprob, QNDF()))

                    oprob = ODEProblem(lrs, u0, (0.0,10.0), (pV, pE); sparse=false)
                    @test SciMLBase.successful_retcode(solve(oprob, QNDF()))
                end
            end
        end
    end

    # Non-stiff case
    for srs in [Vector{SpatialReaction}(), binding_srs_1, binding_srs_2, binding_srs_3, binding_srs_4]    
        lrs = LatticeReactionSystem(binding_system, srs, grid)
        u0_1 = [:X => 1.0, :Y => 2.0, :XY => 0.0]
        u0_2 = [:X => rand_v_vals(grid), :Y => 2.0, :XY => 0.0]
        u0_3 = [:X => 1.0, :Y => rand_v_vals(grid), :XY => rand_v_vals(grid)]
        u0_4 = [:X => rand_v_vals(grid), :Y => rand_v_vals(grid), :XY => rand_v_vals(grid,3)]
        u0_5 = make_u0_matrix(u0_3, vertices(grid), map(s -> Symbol(s.f), species(lrs.rs)))
        for u0 in [u0_1, u0_2, u0_3, u0_4, u0_5]
            p1 = [:kB => 2.0, :kD => 0.5]
            p2 = [:kB => 2.0, :kD => rand_v_vals(grid)]
            p3 = [:kB => rand_v_vals(grid), :kD => rand_v_vals(grid)]
            p4 = make_u0_matrix(p1, vertices(grid), Symbol.(parameters(lrs.rs)))
            for pV in [p1, p2, p3, p4]
                pE_1 = map(sp -> sp => 0.2, lrs.spatial_params)
                pE_2 = map(sp -> sp => rand(), lrs.spatial_params)
                pE_3 = map(sp -> sp => rand_e_vals(grid, 0.2), lrs.spatial_params)
                pE_4 = make_u0_matrix(pE_3, edges(grid), lrs.spatial_params)
                for pE in [pE_1, pE_2, pE_3, pE_4]
                    oprob = ODEProblem(lrs, u0, (0.0,10.0), (pV, pE))
                    @test SciMLBase.successful_retcode(solve(oprob, Tsit5()))

                    oprob = ODEProblem(lrs, u0, (0.0,10.0), (pV, pE); jac=false)
                    @test SciMLBase.successful_retcode(solve(oprob, Tsit5()))
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
    
    runtime_target = 0.00089
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
    
    runtime_target = 0.451
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
    
    runtime_target = 0.05
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
    
    runtime_target = 140.0
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
    
    runtime_target = 0.00293
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
    
    runtime_target = 1.257
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
    
    runtime_target = 0.023
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
    
    runtime_target = 111.
    runtime = minimum((@benchmark solve($oprob, QNDF())).times)/1000000000
    println("Large grid, mid-sized, non-stiff, system. Runtime: $(runtime), previous standard: $(runtime_target)")
    @test runtime < 1.2*runtime_target
end

states(sigmaB_system)

sigmaB_dif_σB = DiffusionReaction(:DσB, :σB)
sigmaB_dif_w = DiffusionReaction(:Dw, :w)
sigmaB_dif_v = DiffusionReaction(:v, :v)
sigmaB_srs_1 = [sigmaB_dif_σB]
sigmaB_srs_2 = [sigmaB_dif_σB, sigmaB_dif_w, sigmaB_dif_v]


CuH_Amination_diff_1 = DiffusionReaction(:D1, :CuoAc)
CuH_Amination_diff_2 = DiffusionReaction(:D2, :Silane)
CuH_Amination_diff_3 = DiffusionReaction(:D3, :Cu_ELigand)
CuH_Amination_diff_4 = DiffusionReaction(:D4, :Amine)
CuH_Amination_diff_5 = DiffusionReaction(:D5, :CuHLigand)


lrs = LatticeReactionSystem(binding_system, brusselator_srs_2, small_2d_grid)
u0 = [:X => rand_v_vals(lrs.lattice), :Y => rand_v_vals(lrs.lattice), :XY => rand_v_vals(lrs.lattice)]
pV = binding_p
pE = [:dX => 0.1, :dY => 0.2]
oprob = ODEProblem(lrs, u0, (0.0,10.0), (pV,pE); jac=false)

runtime = minimum((@benchmark solve($oprob, Tsit5())).times)

(@benchmark solve(oprob, Tsit5())).times

bm  = @benchmark solve(oprob, Tsit5())

median(bm.times)





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


length(edges(small_directed_cycle))

typeof(zeros(0,3))

pE = matrix_form(pE_in, nE, pE_idxes)

zeros(1,0)

pE_1 = map(sp -> sp => 0.2, lrs.spatial_params)

lrs.spatial_reactions

using Test

grid = [small_2d_grid, short_path, small_directed_cycle][1]
srs = [Vector{SpatialReaction}(), binding_srs_1, binding_srs_2, binding_srs_3, binding_srs_4][2]
lrs = LatticeReactionSystem(binding_system, srs, grid)

SciMLBase.successful_retcode(sol)
sol.retcode

grid = [small_2d_grid, short_path, small_directed_cycle][1]
srs = [Vector{SpatialReaction}(), binding_srs_1, binding_srs_2, binding_srs_3, binding_srs_4][2]
lrs = LatticeReactionSystem(binding_system, srs, grid)

u0_1 = [:X => 1.0, :Y => 2.0, :XY => 0.0]
u0_2 = [:X => rand_v_vals(grid), :Y => 2.0, :XY => 0.0]
u0_3 = [:X => 1.0, :Y => rand_v_vals(grid), :XY => rand_v_vals(grid)]
u0_4 = [:X => rand_v_vals(grid), :Y => rand_v_vals(grid), :XY => rand_v_vals(grid,3)]
u0_5 = make_u0_matrix(u0_3, grid, map(s -> Symbol(s.f), species(lrs.rs)))
u0 = [u0_1, u0_2, u0_3, u0_4, u0_5][3]

p1 = [:kB => 2.0, :kD => 0.5]
p2 = [:kB => 2.0, :kD => rand_v_vals(grid)]
p3 = [:kB => rand_v_vals(grid), :kD => rand_v_vals(grid)]
p4 = make_u0_matrix(p1, grid, Symbol.(parameters(lrs.rs)))
p = p3

pE_1 = map(sp -> sp => 0.2, lrs.spatial_params)
pE_2 = map(sp -> sp => rand(), lrs.spatial_params)
pE_3 = map(sp -> sp => rand_e_vals(grid, 0.2), lrs.spatial_params)
pE_4 = make_u0_matrix(pE_3, grid, lrs.spatial_params)
pE = pE_3

oprob = ODEProblem(lrs, u0, (0.0,10.0), (p, pE))
sol = solve(oprob, Tsit5())







rs = @reaction_network begin
    A, ∅ → X
    1, 2X + Y → 3X
    B, X → Y
    1, X → ∅
end A B
spatial_reactions = [SpatialReaction(:D, ([:X], []), ([], [:X]), ([1], Vector{Int64}()), (Vector{Int64}(), [1]))]
lattice = Graphs.grid([20, 20])
lrs = LatticeReactionSystem(rs, spatial_reactions, lattice)

u0 = [:X => 10 * rand(nv(lattice)), :Y => 10 * rand(nv(lattice))]
tspan = (0.0, 100.0)
p = [:A => 1.0, :B => 4.0, :D => 0.2]

oprob = ODEProblem(lrs, u0, (0.0,100.0), p)
@time sol = solve(oprob, Tsit5())
@time sol = solve(oprob, QNDF())

plot(sol,idxs=[1,11])

@profview  sol = solve(oprob, Tsit5())
@profview  sol = solve(oprob, Rosenbrock23())

@btime solve(oprob, Tsit5())
@btime solve(oprob, QNDF())


sizes = 10:10:200
explit_btimes = []
implicit_btimes = []
@time for s in sizes
    println("Running for size $s.")
    rs = @reaction_network begin
        A, ∅ → X
        1, 2X + Y → 3X
        B, X → Y
        1, X → ∅
    end A B
    spatial_reactions = [SpatialReaction(:D, ([:X], []), ([], [:X]), ([1], Vector{Int64}()), (Vector{Int64}(), [1]))]
    lattice = Graphs.grid([10, s])
    lrs = LatticeReactionSystem(rs, spatial_reactions, lattice)

    u0 = [:X => 10 * rand(nv(lattice)), :Y => 10 * rand(nv(lattice))]
    tspan = (0.0, 100.0)
    p = [:A => 1.0, :B => 4.0, :D => 0.2]

    oprob = ODEProblem(lrs, u0, (0.0,100.0), p)
    b1 = (@benchmark  solve(oprob, Tsit5()))
    b2 = (@benchmark  solve(oprob, QNDF()))
    push!(explit_btimes, median(b1.times))
    push!(implicit_btimes, median(b2.times))
    println("Explicit runtim: $(median(b1.times)/1000000000)")
    println("Implicit runtim: $(median(b2.times)/1000000000)\n")
end

plot(sizes,explit_btimes ./1000000000; label="Explicit runtimes")
plot(sizes,implicit_btimes ./1000000000; label="Implicit runtimes")
plot!(yaxis=:log10)


plot(sizes,implicit_btimes ./ explit_btimes; label="Explicit runtimes")


plot((explit_btimes ./1000000000) ./ sizes; label="Explicit runtimes")
plot((implicit_btimes ./1000000000)./ sizes; label="Implicit runtimes")


u_tmp = [1.0, 2.0]
p = [1.0, 4.0]
ofunc = ODEFunction(convert(ODESystem, lrs.rs); jac=true)
J = zeros(2,2)
ofunc(J,u_tmp,p,0.0)
J

using SparseArrays
J = [1.0 0; 1 0]
J2 = sparse(J)
ofunc(J2,u_tmp,p,0.0)
J2










spatial_params = unique(getfield.(spatial_reactions, :rate))
pV_in, pE_in = Catalyst.split_parameters(p, spatial_params)
u_idxs = Dict(reverse.(enumerate(Symbolics.getname.(states(rs)))))
pV_idxes = Dict(reverse.(enumerate(Symbol.(parameters(rs)))))
pE_idxes = Dict(reverse.(enumerate(spatial_params)))

nS,nV,nE = length.([states(lrs.rs), vertices(lrs.lattice), edges(lrs.lattice)])
u0 = Vector(reshape(Catalyst.matrix_form(u0, nV, u_idxs),1:nS*nV))

pV = Catalyst.matrix_form(pV_in, nV, pV_idxes)
pE = Catalyst.matrix_form(pE_in, nE, pE_idxes)

ofun = Catalyst.build_odefunction(lrs, true, spatial_params)

supertype(typeof(ofun))
ofun isa SciMLBase.AbstractODEFunction{true}



matrix_form(input::Matrix, args...) = input
function matrix_form(input::Vector{Pair{Symbol, Vector{Float64}}}, n, index_dict)
    mapreduce(permutedims, vcat, last.(sort(input, by = i -> index_dict[i[1]])))
end
function matrix_form(input::Vector, n, index_dict)
    matrix_form(map(i -> (i[2] isa Vector) ? i[1] => i[2] : i[1] => fill(i[2], n), input),
                n, index_dict)
end

1































spatial_params = unique(getfield.(spatial_reactions, :rate))
pV_in, pE_in = split_parameters(p, spatial_params)
u_idxs = Dict(reverse.(enumerate(Symbolics.getname.(states(rs)))))
pV_idxes = Dict(reverse.(enumerate(Symbol.(parameters(rs)))))
pE_idxes = Dict(reverse.(enumerate(spatial_params)))

nS, nV, nE = get_sizes(lrs)
u0 = Vector(reshape(matrix_form(u0, nV, u_idxs),1:nS*nV))

pV = matrix_form(pV_in, nV, pV_idxes)
pE = matrix_form(pE_in, nE, pE_idxes)

# As a spatial reaction, but replaces the species (and parameter) symbols with their index.
struct SpatialReactionIndexed
    rate::Int64
    substrates::Tuple{Vector{Int64}, Vector{Int64}}
    products::Tuple{Vector{Int64}, Vector{Int64}}
    substoich::Tuple{Vector{Int64}, Vector{Int64}}
    prodstoich::Tuple{Vector{Int64}, Vector{Int64}}
    netstoich::Tuple{Vector{Pair{Int64,Int64}}, Vector{Pair{Int64,Int64}}}
    only_use_rate::Bool

    function SpatialReactionIndexed(sr::SpatialReaction, species_list::Vector{Symbol}, param_list::Vector{Symbol})
        get_s_idx(species::Symbol) = findfirst(species .== (species_list))
        rate = findfirst(sr.rate .== (param_list))
        substrates = Tuple([get_s_idx.(sr.substrates[i]) for i in 1:2])
        products = Tuple([get_s_idx.(sr.products[i]) for i in 1:2])
        netstoich = Tuple([Pair.(get_s_idx.(first.(sr.netstoich[i])), last.(sr.netstoich[i])) for i in 1:2])
        new(rate, substrates, products, sr.substoich, sr.prodstoich, netstoich, sr.only_use_rate)
    end
end



ofunc = ODEFunction(convert(ODESystem, lrs.rs))
nS,nV,nE = length.([states(lrs.rs), vertices(lrs.lattice), edges(lrs.lattice)])
srs_idxed = [SpatialReactionIndexed(sr, Symbol.(getfield.(states(rs), :f)), spatial_params) for sr in spatial_reactions]

f = build_f(ofunc,nS,nV,srs_idxed)
jac = build_jac(ofunc,nS,nV,srs_idxed)
jac_prototype = build_jac_prototype(nS,nV,srs_idxed)

ofun = ODEFunction(f; jac=jac, jac_prototype=(true ? sparse(jac_prototype) : jac_prototype))

ODEProblem(ofun, u0, tspan, (pV, pE))

function build_f(ofunc,nS,nV,srs_idxed)

    return function(du, u, p, t)
        # Updates for non-spatial reactions.
        for comp_i in 1:nV
            ofunc((@view du[get_indexes(comp_i,nS)]), (@view u[get_indexes(comp_i,nS)]), (@view p[1][:,comp_i]), t)
        end
    
        # Updates for spatial reactions.
        for comp_i in 1:nV
            for comp_j in (lrs.lattice.fadjlist)[comp_i], sr in srs_idxed
                rate = get_rate(sr, p[2], (@view u[get_indexes(comp_i,nS)]), (@view u[get_indexes(comp_j,nS)]))
                for (comp,idx) in [(comp_i,1),(comp_j,2)]
                    for stoich in sr.netstoich[idx]
                        du[get_index(comp,stoich[idx],nS)] += rate * stoich[2]
                    end
                end
            end
        end
    end
end

function build_jac(ofunc,nS,nV,srs_idxed)
    base_zero = zeros(nS*nV,nS*nV)

    return function(J, u, p, t)
        J .= base_zero

        # Updates for non-spatial reactions.
        for comp_i in 1:nV
            ofunc.jac((@view J[get_indexes(comp_i,nS),get_indexes(comp_i,nS)]), (@view u[get_indexes(comp_i,nS)]), (@view p[1][:,comp_i]), t)
        end
    
        # Updates for spatial reactions.
        for comp_i in 1:nV
            for comp_j in (lrs.lattice.fadjlist)[comp_i], sr in srs_idxed
                for (idx1, comp1) in [(1,comp_i),(2,comp_j)], sub in sr.substrates[idx1]
                    rate = get_rate_differential(sr, p[2], sub, (@view u[get_indexes(comp_i,nS)]), (@view u[get_indexes(comp_j,nS)]))
                    for (idx2, comp2) in [(1,comp_i),(2,comp_j)], stoich in sr.netstoich[idx2]
                        J[get_index(comp2,stoich[1],nS,u_idxs),get_index(comp1,sub,nS,u_idxs)] += rate * stoich[2]
                    end
                end
            end
        end
    end
end

function get_rate_differential(sr, pE, diff_species, u_src, u_dst)
    product = pE[sr.rate]
    !isempty(sr.substrates[1]) && for (sub,stoich) in zip(sr.substrates[1], sr.substoich[1])
        if diff_species==sub
            product *= stoich*u_src[sub]^(stoich-1) / factorial(stoich)
        else
            product *= u_src[sub]^stoich / factorial(stoich)
        end
    end
    !isempty(sr.substrates[2]) && for (sub,stoich) in zip(sr.substrates[2], sr.substoich[2])
        product *= u_dst[sub]^stoich / factorial(stoich)
    end
    return product
end

findfirst(isequal(reactions(rs)[1].netstoich[1][1]), states(rs))

function build_jac_prototype(nS,nV,srs_idxed)
    jac_prototype = zeros(nS*nV,nS*nV)

    # Sets non-spatial reactions.
    for comp_i in 1:nV, reaction in reactions(lrs.rs)        
        for substrate in reaction.substrates, ns in reaction.netstoich
            sub_idx = findfirst(isequal(substrate), states(lrs.rs))
            spec_idx = findfirst(isequal(ns[1]), states(lrs.rs))
            jac_prototype[spec_idx, sub_idx] = 1
        end
    end

    for comp_i in 1:nV
        for comp_j in (lrs.lattice.fadjlist)[comp_i], sr in srs_idxed
            for (idx1, comp1) in [(1,comp_i),(2,comp_j)], sub in sr.substrates[idx1]
                for (idx2, comp2) in [(1,comp_i),(2,comp_j)], stoich in sr.netstoich[idx2]
                    jac_prototype[get_index(comp2,stoich[1],nS),get_index(comp1,sub,nS)] = 1
                end
            end
        end
    end
    
    return jac_prototype
end


# Get the rate of a specific reaction.
function get_rate(sr, pE, u_src, u_dst)
    product = pE[sr.rate]
    for (u,idx) in [(u_src,1),(u_dst,2)]
        !isempty(sr.substrates[idx]) && for (sub,stoich) in zip(sr.substrates[idx], sr.substoich[idx])
            product *= u[sub]^stoich / factorial(stoich)
        end
    end
    return product
end


function build_odefunction(lrs::LatticeReactionSystem; sparse=true)
    ofunc = ODEFunction(convert(ODESystem, lrs.rs))
    nS,nV,nE = length.([states(lrs.rs), vertices(lrs.lattice), edges(lrs.lattice)])
    sr_info = [(rate=get_rate(sr, lrs.spatial_reactions), substrate=get_substrate_idxs(sr,lrs.rs), netstoich=get_netstoich_idxs(sr,lrs.rs)) for sr in lrs.spatial_reactions]

    f = build_f(ofunc,nS,nV,sr_info)
    jac = build_jac(ofunc,nS,nV,sr_info)
    jac_prototype = build_jac_prototype(ofunc,nS,nV,sr_info)

    return ODEFunction(f; jac=jac, jac_prototype=(sparse ? sparse(jac_prototype) : jac_prototype))
end

get_index(container::Int64,species::Int64,nS) = (container-1)*nS + species
get_indexes(container::Int64,nS) = (container-1)*nS+1:container*nS


# Splits parameters into those for the compartments and those for the connections.
split_parameters(parameters::Tuple, spatial_params) = parameters
function split_parameters(parameters::Vector, spatial_params)
    filter(p -> !in(p[1], spatial_params), parameters),
    filter(p -> in(p[1], spatial_params), parameters)
end

# Converts species and parameters to matrices form.
matrix_form(input::Matrix, args...) = input
function matrix_form(input::Vector{Pair{Symbol, Vector{Float64}}}, n, index_dict)
    mapreduce(permutedims, vcat, last.(sort(input, by = i -> index_dict[i[1]])))
end
function matrix_form(input::Vector, n, index_dict)
    matrix_form(map(i -> (i[2] isa Vector) ? i[1] => i[2] : i[1] => fill(i[2], n), input),
                n, index_dict)
end

get_sizes(lrs) = length.([states(lrs.rs), vertices(lrs.lattice), edges(lrs.lattice)])

get_index(container::Int64,species::Int64,nS) = (container-1)*nS + species
get_indexes(container::Int64,nS) = (container-1)*nS+1:container*nS








































fieldnames(typeof(reactions(rs)[1]))
reactions(rs)[1].substrates
first.(reactions(rs)[1].netstoich)[1]

fieldnames(typeof(rs))

fieldnames(typeof(first.(reactions(rs)[1].netstoich)[1]))

reaction.substrates

to_sym()
sub_syms = 
for comp_i in 1:nV, reaction in reactions(lrs.rs)
    for substrate in Symbol.(getfield.(reaction.substrates, :f)), reactant in Symbol.(getfield.(first.(reactions(rs)[1].netstoich), :f))
        jac_p[u_idxs[reactant],u_idxs[substrate]] = 1
    end
end

for comp_i in 1:nV
    for comp_j in (lrs.lattice.fadjlist)[comp_i], sr in lrs.spatial_reactions
        for substrate in sr.substrates[1], reactant in first.(sr.netstoich[1])

        end
    end
end

for sub in sr.substrates[1]
    
    for  ns in sr.netstoich[1]
        jpp[get_index(comp_i,ns[1],nS,u_idxs),get_index(comp_i,sub,nS,u_idxs)] = 1
    end
end
for sub in sr.substrates[1], ns in sr.netstoich[2]
    jpp[get_index(comp_j,ns[1],nS,u_idxs),get_index(comp_i,sub,nS,u_idxs)] = 1
end


Symbol(first.(reactions(rs)[1].netstoich)[1].f)

reactions(lrs.rs)[3].substrates

spatial_params = unique(getfield.(spatial_reactions, :rate))
pV_in, pE_in = split_parameters(p, spatial_params)
nV, nE = length.([vertices(lattice), edges(lattice)])
u_idxs = Dict(reverse.(enumerate(Symbolics.getname.(states(rs)))))
pV_idxes = Dict(reverse.(enumerate(Symbol.(parameters(rs)))))
pE_idxes = Dict(reverse.(enumerate(spatial_params)))

u0

u0 = matrix_form(u0, nV, u_idxs)
pV = matrix_form(pV_in, nV, pV_idxes)
pE = matrix_form(pE_in, nE, pE_idxes)

u0_vec = Vector(reshape(u0,1:prod(size(u0))))

ofun = ODEFunction(build_f(lrs, u_idxs, pE_idxes); jac=build_jac(lrs, u_idxs, pE_idxes), jac_prototype=build_jac_prototype(lrs, u_idxs, pE_idxes))

split_parameters(parameters::Tuple, spatial_params) = parameters
function split_parameters(parameters::Vector, spatial_params)
    filter(p -> !in(p[1], spatial_params), parameters),
    filter(p -> in(p[1], spatial_params), parameters)
end;

# Converts species and parameters to matrices form.
matrix_form(input::Matrix, args...) = input
function matrix_form(input::Vector{Pair{Symbol, Vector{Float64}}}, n, index_dict)
    mapreduce(permutedims, vcat, last.(sort(input, by = i -> index_dict[i[1]])))
end
function matrix_form(input::Vector, n, index_dict)
    matrix_form(map(i -> (i[2] isa Vector) ? i[1] => i[2] : i[1] => fill(i[2], n), input),
                n, index_dict)
end;



oprob = ODEProblem(lrs, u0, (0.0,100.0), p)
@btime sol = solve(oprob, Tsit5())
@btime sol = solve(oprob, Rosenbrock23())


@profview solve(oprob, Tsit5())
@profview solve(oprob, Rosenbrock23())

sizes = 10:10:100
explit_btimes = []
implicit_btimes = []
for s in sizes
    rs = @reaction_network begin
        A, ∅ → X
        1, 2X + Y → 3X
        B, X → Y
        1, X → ∅
    end A B
    spatial_reactions = [SpatialReaction(:D, ([:X], []), ([], [:X]), ([1], Vector{Int64}()), (Vector{Int64}(), [1]))]
    lattice = Graphs.grid([10, s])
    lrs = LatticeReactionSystem(rs, spatial_reactions, lattice)
    
    u0 = [:X => 10 * rand(nv(lattice)), :Y => 10 * rand(nv(lattice))]
    tspan = (0.0, 100.0)
    p = [:A => 1.0, :B => 4.0, :D => 0.2]
    
    oprob = ODEProblem(lrs, u0, tspan, p)
    b1 = (@benchmark  solve(oprob, Tsit5()))
    b2 = (@benchmark  solve(oprob, Rosenbrock23()))
    push!(explit_btimes, median(b1.times))
    push!(implicit_btimes, median(b2.times))
end


plot(sizes,explit_btimes; label="Explicit runtimes")
plot!(sizes,implicit_btimes; label="Implicit runtimes")

s = 10
rs = @reaction_network begin
    A, ∅ → X
    1, 2X + Y → 3X
    B, X → Y
    1, X → ∅
end A B
spatial_reactions = [SpatialReaction(:D, ([:X], []), ([], [:X]), ([1], Vector{Int64}()), (Vector{Int64}(), [1]))]
lattice = Graphs.grid([s, s])
lrs = LatticeReactionSystem(rs, spatial_reactions, lattice)

u0 = [:X => 10 * rand(nv(lattice)), :Y => 10 * rand(nv(lattice))]
tspan = (0.0, 100.0)
p = [:A => 1.0, :B => 4.0, :D => 0.2]

oprob = ODEProblem(lrs, u0, tspan, p)
b1 = (@benchmark  solve(oprob, Tsit5()))

b2 = (@benchmark  solve(oprob, Rosenbrock23()))

explit_btimes = map(n -> b_timings(n, Tsit5()), N)
implicit_btimes = map(n -> b_timings(n, Rosenbrock23()), N)

explit_btimes
implicit_btimes

bts = b_timings(100, Tsit5())


plot(N,explit_btimes; label="Explicit runtimes")
plot!(N,implicit_btimes; label="Implicit runtimes")

benchmark_timings

expliti_btimes

@btime sleep(1)

b = (@benchmark sleep(0.1))
median(b.times)
fieldnames(typeof(b))

@time sol = solve(oprob, Rosenbrock23())
plot(sol; idxs=[1,11])
@time sol = solve(oprob, Tsit5())
plot(sol; idxs=[1,11])

@profview solve(oprob, Rosenbrock23())
@profview solve(oprob, Tsit5())

@time sol = solve(oprob, Tsit5())
@profview solve(oprob, Tsit5())





spatial_params = unique(getfield.(spatial_reactions, :rate))
pV_in, pE_in = Catalyst.split_parameters(p, spatial_params)
nV, nE = length.([vertices(lattice), edges(lattice)])
u_idxs = Dict(reverse.(enumerate(Symbolics.getname.(states(rs)))))
pV_idxes = Dict(reverse.(enumerate(Symbol.(parameters(rs)))))
pE_idxes = Dict(reverse.(enumerate(spatial_params)))

u0_mat = Catalyst.matrix_form(u0, nV, u_idxs)
u0_vec = Vector(reshape(u0_mat,1:prod(size(u0_mat))))
pV = Catalyst.matrix_form(pV_in, nV, pV_idxes)
pE = Catalyst.matrix_form(pE_in, nE, pE_idxes)



ofunc = ODEFunction(convert(ODESystem, lrs.rs); jac=true)
u = deepcopy(u0_vec)
du = zeros(length(u0_vec))

nS = length(states(rs))

return function internal___spatial___f(du, u, p, t)
    # Updates for non-spatial reactions.
    for comp_i in 1:nV
        ofunc((@view du[get_indexes(comp_i)]), (@view u[get_indexes(comp_i)]), p[1][:,comp_i], t)
    end

    # Updates for spatial reactions.
    for comp_i in 1:nV
        for comp_j::Int64 in (lrs.lattice.fadjlist::Vector{Vector{Int64}})[comp_i],
            sr::SpatialReaction in lrs.spatial_reactions::Vector{SpatialReaction}
            rate = get_rate(sr, p[2], (@view u[get_indexes(comp_i)]), (@view u[get_indexes(comp_j)]), u_idxs, pE_idxes)
        
            for stoich in sr.netstoich[1]
                du[get_index(comp_i,stoich[1])] += rate * stoich[2]
            end
            for stoich in sr.netstoich[2]
                du[get_index(comp_j,stoich[1])] += rate * stoich[2]
            end
        end
    end
end
function get_rate(sr, pE, u_src, u_dst, u_idxs, pE_idxes)
    product = pE[pE_idxes[sr.rate]]
    !isempty(sr.substrates[1]) && for (sub,stoich) in zip(sr.substrates[1], sr.substoich[1])
        product *= u_src[u_idxs[sub]]^stoich / factorial(stoich)
    end
    !isempty(sr.substrates[2]) && for (sub,stoich) in zip(sr.substrates[2], sr.substoich[2])
        product *= u_dst[u_idxs[sub]]^stoich / factorial(stoich)
    end
    return product
end

fieldnames(typeof(reactions(rs)[1]))

oprob = ODEProblem(internal___spatial___f, u0_vec, (0.0,100.0), (pV, pE))
@time sol = solve(oprob, Rosenbrock23())

plot(sol; idxs=[1,11])

plot(getindex.((sol.u),1))
plot!(getindex.((sol.u),3))

reactions(rs)[2].substrates isa Vector{<:Term}

typeof(reactions(rs)[2].substrates[1])
findfirst(isequal.(reactions(rs)[2].substrates[1],states(rs)))

syms_to_idxs(reactions(rs)[2].substrates, states(rs))
reactions(rs)[2].substrates

zip([1,3], spatial_reactions.netstoich)

spatial_reactions[1].netstoich

odelingToolkit.var_
spatial_reactions[1].netstoich

spatial_reactions[1].substrates

(1 for i in 1:2)

sr = spatial_reactions[1]



syms_to_idxs(sr.substrates[2], states(lrs.rs))

# For a vector of Symbolics or Symbols, find their indexes in an array.
syms_to_idxs(syms, syms_reference::Vector{<:Term}) = [sym_to_idxs(sym, syms_reference) for sym in syms]::Vector{Int64}
# For a Symbolic or Symbol, find its index in an array.
sym_to_idxs(sym::Term, syms_reference::Vector{<:Term}) = findfirst(isequal.(sym, syms_reference))::Int64
sym_to_idxs(sym::Symbol, syms_reference::Vector{<:Term}) = findfirst(isequal.(sym, Symbol.(getfield.(syms_reference,:f))))::Int64




function build_odefunction(lrs::LatticeReactionSystem; sparse=true)
    ofunc = ODEFunction(convert(ODESystem, lrs.rs))
    nS,nV,nE = length.([states(lrs.rs), vertices(lrs.lattice), edges(lrs.lattice)])
    sr_info = [(rate=get_rate(sr, lrs.spatial_reactions), substrate=get_substrate_idxs(sr,lrs.rs), netstoich=get_netstoich_idxs(sr,lrs.rs)) for sr in lrs.spatial_reactions]

    f = build_f(ofunc,nS,nV,sr_info)
    jac = build_jac(ofunc,nS,nV,sr_info)
    jac_prototype = build_jac_prototype(ofunc,nS,nV,sr_info)

    return ODEFunction(f; jac=jac, jac_prototype=(sparse ? sparse(jac_prototype) : jac_prototype))
end




# Creates a function for simulating the spatial ODE with spatial reactions.
function build_f(ofunc, nS, nV, sr_info)

    return function(du, u, p, t)
        # Updates for non-spatial reactions.
        for comp_i in 1:nV
            ofunc((@view du[get_indexes(comp_i,nS)]), (@view u[get_indexes(comp_i,nS)]), (@view p[1][:,comp_i]), t)
        end
    
        # Updates for spatial reactions.
        for comp_i in 1:nV
            for comp_j::Int64 in (lrs.lattice.fadjlist::Vector{Vector{Int64}})[comp_i],
                sr in sr_info




                rate = get_rate(sr, p[2], (@view u[get_indexes(comp_i,nS)]), (@view u[get_indexes(comp_j,nS)]))
            
                for stoich in sr.netstoich[1]
                    du[get_index(comp_i,stoich[1],nS,u_idxs)] += rate * stoich[2]
                end
                for stoich in sr.netstoich[2]
                    du[get_index(comp_j,stoich[1],nS,u_idxs)] += rate * stoich[2]
                end
            end
        end
    end
end

function get_rate(sr, pE, u_src, u_dst, pE_idxes)
    product = pE[sr.rate]
    !isempty(sr.substrates[1]) && for (sub,stoich) in zip(sr.substrates[1], sr.substoich[1])
        product *= u_src[sub]^stoich / factorial(stoich)
    end
    !isempty(sr.substrates[2]) && for (sub,stoich) in zip(sr.substrates[2], sr.substoich[2])
        product *= u_dst[u_idxs[sub]]^stoich / factorial(stoich)
    end
    return product
end


return function internal___spatial___jac(J, u, p, t)
    J .= zeros(length(u),length(u))

    # Updates for non-spatial reactions.
    for comp_i in 1:nV
        ofunc.jac((@view J[get_indexes(comp_i),get_indexes(comp_i)]), (@view u[get_indexes(comp_i)]), p[1][:,comp_i], t)
    end

    # Updates for spatial reactions.
    for comp_i in 1:nV
        for comp_j::Int64 in (lrs.lattice.fadjlist::Vector{Vector{Int64}})[comp_i],
            sr::SpatialReaction in lrs.spatial_reactions::Vector{SpatialReaction}

            for sub in sr.substrates[1]
                rate = get_rate_differential(sr, p[2], sub, (@view u[get_indexes(comp_i)]), (@view u[get_indexes(comp_j)]), u_idxs, pE_idxes)
                for stoich in sr.netstoich[1]
                    J[get_index(comp_i,stoich[1]),get_index(comp_i,sub)] += rate * stoich[2]
                end
                for stoich in sr.netstoich[2]
                    J[get_index(comp_j,stoich[1]),get_index(comp_i,sub)] += rate * stoich[2]
                end
            end

            for sub in sr.substrates[2]
                rate = get_rate_differential(sr, p[2], sub, (@view u[get_indexes(comp_j)]), (@view u[get_indexes(comp_i)]), u_idxs, pE_idxes)
                for stoich in sr.netstoich[1]
                    J[get_index(comp_i,stoich[1]),get_index(comp_j,sub)] += rate * stoich[2]
                end
                for stoich in sr.netstoich[2]
                    J[get_index(comp_j,stoich[1]),get_index(comp_j,sub)] += rate * stoich[2]
                end
            end
        end
    end
end
get_index(container::Int64,species::Symbol) = (container-1)*nS + u_idxs[species]
get_indexes(container::Int64) = (container-1)*nS+1:container*nS

# Get the rate differential of a specific reaction.
function get_rate_differential(sr, pE, diff_species, u_src, u_dst, u_idxs, pE_idxes)
    product = pE[pE_idxes[sr.rate]]
    !isempty(sr.substrates[1]) && for (sub,stoich) in zip(sr.substrates[1], sr.substoich[1])
        (diff_species==sub) && (product *= stoich*u_src[u_idxs[sub]]^(stoich-1) / factorial(stoich))
        (diff_species!=sub) && (product *= u_src[u_idxs[sub]]^stoich / factorial(stoich))
    end
    !isempty(sr.substrates[2]) && for (sub,stoich) in zip(sr.substrates[2], sr.substoich[2])
        product *= u_dst[u_idxs[sub]]^stoich / factorial(stoich)
    end
    return product
end


jpp = zeros(nS*nV,nS*nV)
foreach(i -> jpp[(i-1)*nS+1:i*nS,(i-1)*nS+1:i*nS] = ones(nS,nS), 1:nV)
for comp_i in 1:nV
    for comp_j::Int64 in (lrs.lattice.fadjlist::Vector{Vector{Int64}})[comp_i],
        sr::SpatialReaction in lrs.spatial_reactions::Vector{SpatialReaction}

        for sub in sr.substrates[1], ns in sr.netstoich[1]
            jpp[get_index(comp_i,ns[1]),get_index(comp_i,sub)] = 1
        end
        for sub in sr.substrates[1], ns in sr.netstoich[2]
            jpp[get_index(comp_j,ns[1]),get_index(comp_i,sub)] = 1
        end
        for sub in sr.substrates[2], ns in sr.netstoich[1]
            jpp[get_index(comp_i,ns[1]),get_index(comp_j,sub)] = 1
        end
        for sub in sr.substrates[2], ns in sr.netstoich[2]
            jpp[get_index(comp_j,ns[1]),get_index(comp_j,sub)] = 1
        end
    end
end

jac_prototype = sparse(jpp)
ofun_jac = ODEFunction(internal___spatial___f; jac=internal___spatial___jac, jac_prototype=jac_prototype)

oprob = ODEProblem(ofun_jac, u0_vec, (0.0,100.0), (pV, pE))
@time sol = solve(oprob, Rosenbrock23())

plot(sol; idxs=[1,11])

return function jac_2node(J, u, p, t)
    A,B = p[1][:,1]
    D, = p[2][:,1]
    X1,Y1,X2,Y2 = u
    J .= zeros(4,4)

    J[1,1] = X1*Y1 - (1+B+D)
    J[1,2] = 0.5*X1^2
    J[1,3] = D

    J[2,1] = B - X1*Y1
    J[2,2] = -0.5*X1^2

    J[3,1] = D
    J[3,3] = X2*Y2 - (1+B+D)
    J[3,4] = 0.5*X2^2

    J[4,3] = B - X2*Y2
    J[4,4] = -0.5*X2^2
end

ofun_jac = ODEFunction(internal___spatial___f; jac=jac_2node)
oprob = ODEProblem(ofun_jac, u0_vec, (0.0,100.0), (pV, pE))
@time sol = solve(oprob, Rosenbrock23())

plot(getindex.((sol.u),1))
plot!(getindex.((sol.u),3))

print_m(m) = foreach(i -> println(collect(m)[i,:]), 1:size(m)[1])

u_tmp = [2.035941193990462, 2.7893001714815364, 1.650680925465682, 3.1349007948289445]

J1 = zeros(4,4)
jac_2node(J1, u_tmp, (pV,pE), 0.0)
print_m(J1)

J2 = zeros(4,4)
internal___spatial___jac(J2, u_tmp, (pV,pE), 0.0)
print_m(J2)

print_m(J1.-J2)


ofun_jac = ODEFunction(internal___spatial___f; jac=internal___spatial___jac)
oprob = ODEProblem(ofun_jac, u0_vec, (0.0,100.0), (pV, pE))
@time sol = solve(oprob, Rosenbrock23(); maxiters=100)

plot(getindex.((sol.u),1))
plot!(getindex.((sol.u),3))


return function internal___spatial___jac(J, u, p, t)
    # Updates for non-spatial reactions.
    for comp_i in 1:nV
        ofunc.jac((@view J[get_indexes(comp_i),get_indexes(comp_i)]), (@view u[get_indexes(comp_i)]), p[1][:,comp_i], t)
    end

    # Updates for spatial reactions.
    for comp_i in 1:nV
        for comp_j::Int64 in (lrs.lattice.fadjlist::Vector{Vector{Int64}})[comp_i],
            sr::SpatialReaction in lrs.spatial_reactions::Vector{SpatialReaction}

            for sub in sr.substrates[1]
                rate = get_rate_differential(sr, p[2], sub, (@view u[get_indexes(comp_i)]), (@view u[get_indexes(comp_j)]), u_idxs, pE_idxes)
                for stoich in sr.netstoich[1]
                    J[get_index(comp_i,stoich[1]),get_index(comp_i,sub)] += rate * stoich[2]
                end
                for stoich in sr.netstoich[2]
                    J[get_index(comp_j,stoich[1]),get_index(comp_i,sub)] += rate * stoich[2]
                end
            end

            for sub in sr.substrates[2]
                rate = get_rate_differential(sr, p[2], sub, (@view u[get_indexes(comp_j)]), (@view u[get_indexes(comp_i)]), u_idxs, pE_idxes)
                for stoich in sr.netstoich[1]
                    J[get_index(comp_i,stoich[1]),get_index(comp_j,sub)] += rate * stoich[2]
                end
                for stoich in sr.netstoich[2]
                    J[get_index(comp_j,stoich[1]),get_index(comp_j,sub)] += rate * stoich[2]
                end
            end
        end
    end
end

using FiniteDiff
function my_f(x)
    du = zeros(8)
    internal___spatial___f(du, x, (pV,pE), 0.0)
    return du
end
function my_j(x)
    J = zeros(8,8)
    internal___spatial___jac(J, x, (pV,pE), 0.0)
    return J
end


u0_tmp = 100*rand(8)
maximum(FiniteDiff.finite_difference_jacobian(my_f, u0_cpy) .- my_j(u0_cpy))