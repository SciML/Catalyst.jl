using Catalyst, BenchmarkTools
using OrdinaryDiffEq, JumpProcesses, StableRNGs

const SUITE = BenchmarkGroup()

# SIR reaction network
rn = @reaction_network begin
    k1, S + I --> 2I
    k2, I --> R
end

u0 = [:S => 999, :I => 1, :R => 0]
tspan = (0.0, 250.0)
p = [:k1 => 1.0e-4, :k2 => 0.01]

# =============================================================================
# Model construction
# =============================================================================

SUITE["model"] = BenchmarkGroup()

SUITE["model"]["reaction_network"] = @benchmarkable @reaction_network begin
    k1, S + I --> 2I
    k2, I --> R
end
SUITE["model"]["complete"] = @benchmarkable Catalyst.complete($rn)

rn_complete = Catalyst.complete(rn)

# =============================================================================
# Problem construction
# =============================================================================

SUITE["problem"] = BenchmarkGroup()

SUITE["problem"]["odeproblem"] = @benchmarkable ODEProblem(
    $rn_complete, $u0, $tspan, $p
)
SUITE["problem"]["jumpproblem"] = @benchmarkable JumpProblem(
    $rn_complete, $u0, $tspan, $p; aggregator = Direct()
)

# =============================================================================
# Solves
# =============================================================================

SUITE["solve"] = BenchmarkGroup()

odeprob = ODEProblem(rn_complete, u0, tspan, p)
jprob = JumpProblem(
    rn_complete, u0, tspan, p; aggregator = Direct(), rng = StableRNG(12345)
)

SUITE["solve"]["ode_tsit5"] = @benchmarkable solve($odeprob, Tsit5())
SUITE["solve"]["jump_direct"] = @benchmarkable solve($jprob, SSAStepper())
