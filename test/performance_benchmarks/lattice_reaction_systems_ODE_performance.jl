# Not actually run in CI, but useful for reference of ODE simulation performance across updates.

### Preparations ###

# Fetch packages.
using OrdinaryDiffEq
using Random, Statistics, SparseArrays, Test

# Fetch test networks.
include("../spatial_test_networks.jl")

### Runtime Checks ###
# Current timings are taken from the SciML CI server.
# Current not used, simply here for reference.
# Useful when attempting to optimise workflow.

# using BenchmarkTools, Sundials
# runtime_reduction_margin = 10.0

# Small grid, small, non-stiff, system.
let
    lrs = LatticeReactionSystem(SIR_system, SIR_srs_2, small_2d_graph_grid)
    u0 = [:S => 990.0, :I => 20.0 * rand_v_vals(lrs), :R => 0.0]
    pV = SIR_p
    pE = [:dS => 0.01, :dI => 0.01, :dR => 0.01]
    oprob = ODEProblem(lrs, u0, (0.0, 500.0), [pV; pE]; jac = false)
    @test SciMLBase.successful_retcode(solve(oprob, Tsit5()))

    runtime_target = 0.00027
    runtime = minimum((@benchmark solve($oprob, Tsit5())).times) / 1000000000
    println("Small grid, small, non-stiff, system. Runtime: $(runtime), previous standard: $(runtime_target)")
    @test runtime < runtime_reduction_margin * runtime_target
end

# Large grid, small, non-stiff, system.
let
    lrs = LatticeReactionSystem(SIR_system, SIR_srs_2, large_2d_grid)
    u0 = [:S => 990.0, :I => 20.0 * rand_v_vals(lrs), :R => 0.0]
    pV = SIR_p
    pE = [:dS => 0.01, :dI => 0.01, :dR => 0.01]
    oprob = ODEProblem(lrs, u0, (0.0, 500.0), [pV; pE]; jac = false)
    @test SciMLBase.successful_retcode(solve(oprob, Tsit5()))

    runtime_target = 0.12
    runtime = minimum((@benchmark solve($oprob, Tsit5())).times) / 1000000000
    println("Large grid, small, non-stiff, system. Runtime: $(runtime), previous standard: $(runtime_target)")
    @test runtime < runtime_reduction_margin * runtime_target
end

# Small grid, small, stiff, system.
let
    lrs = LatticeReactionSystem(brusselator_system, brusselator_srs_1, small_2d_graph_grid)
    u0 = [:X => rand_v_vals(lrs, 10), :Y => rand_v_vals(lrs, 10)]
    pV = brusselator_p
    pE = [:dX => 0.2]
    oprob = ODEProblem(lrs, u0, (0.0, 100.0), [pV; pE])
    @test SciMLBase.successful_retcode(solve(oprob, CVODE_BDF(linear_solver=:GMRES)))

    runtime_target = 0.013
    runtime = minimum((@benchmark solve($oprob, CVODE_BDF(linear_solver=:GMRES))).times) / 1000000000
    println("Small grid, small, stiff, system. Runtime: $(runtime), previous standard: $(runtime_target)")
    @test runtime < runtime_reduction_margin * runtime_target
end

# Large grid, small, stiff, system.
let
    lrs = LatticeReactionSystem(brusselator_system, brusselator_srs_1, large_2d_grid)
    u0 = [:X => rand_v_vals(lrs, 10), :Y => rand_v_vals(lrs, 10)]
    pV = brusselator_p
    pE = [:dX => 0.2]
    oprob = ODEProblem(lrs, u0, (0.0, 100.0), [pV; pE])
    @test SciMLBase.successful_retcode(solve(oprob, CVODE_BDF(linear_solver=:GMRES)))

    runtime_target = 11.
    runtime = minimum((@benchmark solve($oprob, CVODE_BDF(linear_solver=:GMRES))).times) / 1000000000
    println("Large grid, small, stiff, system. Runtime: $(runtime), previous standard: $(runtime_target)")
    @test runtime < runtime_reduction_margin * runtime_target
end

# Small grid, mid-sized, non-stiff, system.
let
    lrs = LatticeReactionSystem(CuH_Amination_system, CuH_Amination_srs_2,
                                small_2d_graph_grid)
    u0 = [
        :CuoAc => 0.005 .+ rand_v_vals(lrs, 0.005),
        :Ligand => 0.005 .+ rand_v_vals(lrs, 0.005),
        :CuoAcLigand => 0.0,
        :Silane => 0.5 .+ rand_v_vals(lrs, 0.5),
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
    oprob = ODEProblem(lrs, u0, (0.0, 10.0), [pV; pE]; jac = false)
    @test SciMLBase.successful_retcode(solve(oprob, Tsit5()))

    runtime_target = 0.0012
    runtime = minimum((@benchmark solve($oprob, Tsit5())).times) / 1000000000
    println("Small grid, mid-sized, non-stiff, system. Runtime: $(runtime), previous standard: $(runtime_target)")
    @test runtime < runtime_reduction_margin * runtime_target
end

# Large grid, mid-sized, non-stiff, system.
let
    lrs = LatticeReactionSystem(CuH_Amination_system, CuH_Amination_srs_2,
                                large_2d_grid)
    u0 = [
        :CuoAc => 0.005 .+ rand_v_vals(lrs, 0.005),
        :Ligand => 0.005 .+ rand_v_vals(lrs, 0.005),
        :CuoAcLigand => 0.0,
        :Silane => 0.5 .+ rand_v_vals(lrs, 0.5),
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
    oprob = ODEProblem(lrs, u0, (0.0, 10.0), [pV; pE]; jac = false)
    @test SciMLBase.successful_retcode(solve(oprob, Tsit5()))

    runtime_target = 0.56
    runtime = minimum((@benchmark solve($oprob, Tsit5())).times) / 1000000000
    println("Large grid, mid-sized, non-stiff, system. Runtime: $(runtime), previous standard: $(runtime_target)")
    @test runtime < runtime_reduction_margin * runtime_target
end

# Small grid, mid-sized, stiff, system.
let
    lrs = LatticeReactionSystem(sigmaB_system, sigmaB_srs_2, small_2d_graph_grid)
    u0 = [
        :w => 0.5 .+ rand_v_vals(lrs, 0.5),
        :w2 => 0.5 .+ rand_v_vals(lrs, 0.5),
        :w2v => 0.5 .+ rand_v_vals(lrs, 0.5),
        :v => 0.5 .+ rand_v_vals(lrs, 0.5),
        :w2v2 => 0.5 .+ rand_v_vals(lrs, 0.5),
        :vP => 0.5 .+ rand_v_vals(lrs, 0.5),
        :σB => 0.5 .+ rand_v_vals(lrs, 0.5),
        :w2σB => 0.5 .+ rand_v_vals(lrs, 0.5),
        :vPp => 0.0,
        :phos => 0.4,
    ]
    pV = sigmaB_p
    pE = [:DσB => 0.1, :Dw => 0.1, :Dv => 0.1]
    oprob = ODEProblem(lrs, u0, (0.0, 50.0), [pV; pE])
    @test SciMLBase.successful_retcode(solve(oprob, CVODE_BDF(linear_solver=:GMRES)))

    runtime_target = 0.61
    runtime = minimum((@benchmark solve($oprob, CVODE_BDF(linear_solver=:GMRES))).times) / 1000000000
    println("Small grid, mid-sized, stiff, system. Runtime: $(runtime), previous standard: $(runtime_target)")
    @test runtime < runtime_reduction_margin * runtime_target
end

# Large grid, mid-sized, stiff, system.
let
    lrs = LatticeReactionSystem(sigmaB_system, sigmaB_srs_2, large_2d_grid)
    u0 = [
        :w => 0.5 .+ rand_v_vals(lrs, 0.5),
        :w2 => 0.5 .+ rand_v_vals(lrs, 0.5),
        :w2v => 0.5 .+ rand_v_vals(lrs, 0.5),
        :v => 0.5 .+ rand_v_vals(lrs, 0.5),
        :w2v2 => 0.5 .+ rand_v_vals(lrs, 0.5),
        :vP => 0.5 .+ rand_v_vals(lrs, 0.5),
        :σB => 0.5 .+ rand_v_vals(lrs, 0.5),
        :w2σB => 0.5 .+ rand_v_vals(lrs, 0.5),
        :vPp => 0.0,
        :phos => 0.4,
    ]
    pV = sigmaB_p
    pE = [:DσB => 0.1, :Dw => 0.1, :Dv => 0.1]
    oprob = ODEProblem(lrs, u0, (0.0, 10.0), [pV; pE]) # Time reduced from 50.0 (which casues Julai to crash).
    @test SciMLBase.successful_retcode(solve(oprob, CVODE_BDF(linear_solver=:GMRES)))

    runtime_target = 59.
    runtime = minimum((@benchmark solve($oprob, CVODE_BDF(linear_solver=:GMRES))).times) / 1000000000
    println("Large grid, mid-sized, stiff, system. Runtime: $(runtime), previous standard: $(runtime_target)")
    @test runtime < runtime_reduction_margin * runtime_target
end
