# Fetch packages.
using Catalyst, JumpProcesses, OrdinaryDiffEq
using Random, Statistics, SparseArrays, Test
using Graphs

### Creates a Model ###
SIR_system = @reaction_network begin
    α, S + I --> 2I
    β, I --> R
end

SIR_srs_1 = [DiffusionReaction(:dS, :S)]
SIR_srs_2 = diffusion_reactions([(:dS, :S), (:dI, :I), (:dR, :R)])

grid = Graphs.grid([2, 2])

lrs1 = LatticeReactionSystem(SIR_system, SIR_srs_1, grid)
lrs2 = LatticeReactionSystem(SIR_system, SIR_srs_2, grid)

@test length(compartment_parameters(lrs1)) == length(compartment_parameters(lrs2)) == 2
@test length(diffusion_parameters(lrs1)) == length(diffusion_species(lrs1)) == 1
@test length(diffusion_parameters(lrs2)) == length(diffusion_species(lrs2)) == 3

### Check Parameters and u0 Handling ###
l = Graphs.nv(grid)
SIR_pC = [:α => 0.1 / 1000, :β => 0.01]
SIR_pD = [:dS => 0.01]
SIR_u0 = [:S => 999.0, :I => 1.0, :R => 0.0]
pC_syms = [:α, :β]
pD_syms = [:dS]
u0_syms = [:S, :I, :R]

# Split parameters.
@test Catalyst.split_parameters((SIR_pC, SIR_pD)) == (SIR_pC, SIR_pD)
@test Catalyst.split_parameters([SIR_pC; SIR_pD], [:α, :β], [:dS]) == (SIR_pC, SIR_pD)

# Convert input to vector of vector (of values).
@test Catalyst.lattice_process_input([:S => 1, :I => 2, :R => 3], u0_syms, l) == [[1], [2], [3]]
@test Catalyst.lattice_process_input([:S => 1, :I => [2, 2, 2, 2], :R => 3], u0_syms, l) == [[1], [2, 2, 2, 2], [3]]
@test Catalyst.lattice_process_input([1 1 1 1; 2 2 2 2; 3 3 3 3], u0_syms, l) == [[1, 1, 1, 1], [2, 2, 2, 2], [3, 3, 3, 3]]
@test Catalyst.lattice_process_input([1, 2, 3], u0_syms, l) == [[1], [2], [3]]
@test Catalyst.lattice_process_input([1, [2, 2, 2, 2], 3], u0_syms, l) == [[1], [2, 2, 2, 2], [3]]
@test Catalyst.lattice_process_input([1, [1,2,3,4],1], u0_syms, l)

# Fetching teh diffusion rates.
@test isequal(Catalyst.compute_all_diffusion_rates([[1.0], [2.0, 2.0, 2.0, 2.0], [3.0]], [[0.1]], lrs1), Pair.(diffusion_species(lrs1),[[0.1]]))
@test isequal(Catalyst.compute_all_diffusion_rates([[1.0], [2.0, 2.0, 2.0, 2.0], [3.0]], [[0.1], [0.2], [0.3]], lrs2), Pair.(diffusion_species(lrs2),[[0.1], [0.2], [0.3]]))
