# Not actually run in CI, but useful for reference of ODE simulation performance across updates.

### Preparations ###

# Fetch packages.
using OrdinaryDiffEq
using Random, Statistics, SparseArrays, Test

# Sets rnd number.
using StableRNGs
rng = StableRNG(12345)

# Fetch test networks.
include("spatial_test_networks.jl")

### Runtime Checks ###
# Current timings are taken from the SciML CI server.
# Current not used, simply here for reference.
# Useful when attempting to optimise workflow.

# using BenchmarkTools
# runtime_reduction_margin = 10.0

# Small grid, small, non-stiff, system.
let
    lrs = LatticeReactionSystem(SIR_system, SIR_srs_2, small_2d_grid)
    u0 = [:S => 990.0, :I => 20.0 * rand_v_vals(lrs.lattice), :R => 0.0]
    pV = SIR_p
    pE = [:dS => 0.01, :dI => 0.01, :dR => 0.01]
    oprob = ODEProblem(lrs, u0, (0.0, 500.0), (pV, pE); jac = false)
    @test SciMLBase.successful_retcode(solve(oprob, Tsit5()))

    runtime_target = 0.00060
    runtime = minimum((@benchmark solve($oprob, Tsit5())).times) / 1000000000
    println("Small grid, small, non-stiff, system. Runtime: $(runtime), previous standard: $(runtime_target)")
    @test runtime < runtime_reduction_margin * runtime_target
end

# Large grid, small, non-stiff, system.
let
    lrs = LatticeReactionSystem(SIR_system, SIR_srs_2, large_2d_grid)
    u0 = [:S => 990.0, :I => 20.0 * rand_v_vals(lrs.lattice), :R => 0.0]
    pV = SIR_p
    pE = [:dS => 0.01, :dI => 0.01, :dR => 0.01]
    oprob = ODEProblem(lrs, u0, (0.0, 500.0), (pV, pE); jac = false)
    @test SciMLBase.successful_retcode(solve(oprob, Tsit5()))

    runtime_target = 0.26
    runtime = minimum((@benchmark solve($oprob, Tsit5())).times) / 1000000000
    println("Large grid, small, non-stiff, system. Runtime: $(runtime), previous standard: $(runtime_target)")
    @test runtime < runtime_reduction_margin * runtime_target
end

# Small grid, small, stiff, system.

let
    lrs = LatticeReactionSystem(brusselator_system, brusselator_srs_1, small_2d_grid)
    u0 = [:X => rand_v_vals(lrs.lattice, 10), :Y => rand_v_vals(lrs.lattice, 10)]
    pV = brusselator_p
    pE = [:dX => 0.2]
    oprob = ODEProblem(lrs, u0, (0.0, 100.0), (pV, pE))
    @test SciMLBase.successful_retcode(solve(oprob, QNDF()))

    runtime_target = 0.17
    runtime = minimum((@benchmark solve($oprob, QNDF())).times) / 1000000000
    println("Small grid, small, stiff, system. Runtime: $(runtime), previous standard: $(runtime_target)")
    @test runtime < runtime_reduction_margin * runtime_target
end

# Medium grid, small, stiff, system.
let
    lrs = LatticeReactionSystem(brusselator_system, brusselator_srs_1, medium_2d_grid)
    u0 = [:X => rand_v_vals(lrs.lattice, 10), :Y => rand_v_vals(lrs.lattice, 10)]
    pV = brusselator_p
    pE = [:dX => 0.2]
    oprob = ODEProblem(lrs, u0, (0.0, 100.0), (pV, pE))
    @test SciMLBase.successful_retcode(solve(oprob, QNDF()))

    runtime_target = 2.3
    runtime = minimum((@benchmark solve($oprob, QNDF())).times) / 1000000000
    println("Medium grid, small, stiff, system. Runtime: $(runtime), previous standard: $(runtime_target)")
    @test runtime < runtime_reduction_margin * runtime_target
end

# Large grid, small, stiff, system.
let
    lrs = LatticeReactionSystem(brusselator_system, brusselator_srs_1, large_2d_grid)
    u0 = [:X => rand_v_vals(lrs.lattice, 10), :Y => rand_v_vals(lrs.lattice, 10)]
    pV = brusselator_p
    pE = [:dX => 0.2]
    oprob = ODEProblem(lrs, u0, (0.0, 100.0), (pV, pE))
    @test SciMLBase.successful_retcode(solve(oprob, QNDF()))

    runtime_target = 170.0
    runtime = minimum((@benchmark solve($oprob, QNDF())).times) / 1000000000
    println("Large grid, small, stiff, system. Runtime: $(runtime), previous standard: $(runtime_target)")
    @test runtime < runtime_reduction_margin * runtime_target
end

# Small grid, mid-sized, non-stiff, system.
let
    lrs = LatticeReactionSystem(CuH_Amination_system, CuH_Amination_srs_2,
                                small_2d_grid)
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

    runtime_target = 0.0016
    runtime = minimum((@benchmark solve($oprob, Tsit5())).times) / 1000000000
    println("Small grid, mid-sized, non-stiff, system. Runtime: $(runtime), previous standard: $(runtime_target)")
    @test runtime < runtime_reduction_margin * runtime_target
end

# Large grid, mid-sized, non-stiff, system.
let
    lrs = LatticeReactionSystem(CuH_Amination_system, CuH_Amination_srs_2,
                                large_2d_grid)
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

    runtime_target = 0.67
    runtime = minimum((@benchmark solve($oprob, Tsit5())).times) / 1000000000
    println("Large grid, mid-sized, non-stiff, system. Runtime: $(runtime), previous standard: $(runtime_target)")
    @test runtime < runtime_reduction_margin * runtime_target
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

    runtime_target = 0.019
    runtime = minimum((@benchmark solve($oprob, QNDF())).times) / 1000000000
    println("Small grid, mid-sized, stiff, system. Runtime: $(runtime), previous standard: $(runtime_target)")
    @test runtime < runtime_reduction_margin * runtime_target
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

    runtime_target = 35.0
    runtime = minimum((@benchmark solve($oprob, QNDF())).times) / 1000000000
    println("Large grid, mid-sized, stiff, system. Runtime: $(runtime), previous standard: $(runtime_target)")
    @test runtime < runtime_reduction_margin * runtime_target
end
