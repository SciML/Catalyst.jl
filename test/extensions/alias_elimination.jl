### Prepares Tests ###

# Fetch packages.
using Catalyst, BifurcationKit, Logging, StructuralIdentifiability, Test
import HomotopyContinuation

# A Schlögl-type network (three positive steady states) written with aliases, and by hand.
# The eliminated names (`X`, `k5`) are used in all inputs to the aliased network.
rn = @reaction_network begin
    @parameters k4
    @aliases begin
        X ~ Y
        k5 ~ k4
    end
    k1, 2X --> 3X
    k2, 3Y --> 2Y
    k3, 0 --> X
    k5, Y --> 0
end
rn_merged = @reaction_network begin
    k1, 2Y --> 3Y
    k2, 3Y --> 2Y
    k3, 0 --> Y
    k4, Y --> 0
end
ps = [:k1 => 3.0, :k2 => 0.5, :k3 => 4.0, :k5 => 6.0]
ps_merged = [:k1 => 3.0, :k2 => 0.5, :k3 => 4.0, :k4 => 6.0]


### Tests ###

# HomotopyContinuation: the same steady states, in the resolved system's unknown order (checked
# by passing each to `steady_state_stability`).
let
    sols = sort(hc_steady_states(rn, ps))
    sols_merged = sort(hc_steady_states(rn_merged, ps_merged))
    @test length(sols) == 3
    @test sols ≈ sols_merged
    for (sol, sol_merged) in zip(sols, sols_merged)
        @test steady_state_stability(sol, rn, ps) == steady_state_stability(sol_merged, rn_merged, ps_merged)
    end
end

# BifurcationKit: the bifurcation parameter, plot variable and inputs may use eliminated names.
let
    bprob = BifurcationProblem(rn, [:X => 1.0], ps, :k5; plot_var = :X)
    bprob_merged = BifurcationProblem(rn_merged, [:Y => 1.0], ps_merged, :k4; plot_var = :Y)
    @test BifurcationKit.getparam(bprob) == BifurcationKit.getparam(bprob_merged) == 6.0
    for u in ([0.5], [3.0], [10.0])
        @test bprob.VF.F(u, bprob.params) ≈ bprob_merged.VF.F(u, bprob_merged.params)
    end
end

# StructuralIdentifiability: measured quantities, known parameters and functions to check may
# use eliminated names; results are keyed by the functions as given.
let
    loglevel = Logging.Error
    out = assess_identifiability(rn; measured_quantities = [:X], known_p = [:k5],
        funcs_to_check = [rn.k5, rn.k1, rn.k3], loglevel)
    out_merged = assess_identifiability(rn_merged; measured_quantities = [:Y], known_p = [:k4],
        funcs_to_check = [rn_merged.k4, rn_merged.k1, rn_merged.k3], loglevel)
    @test isequal(collect(keys(out)), [rn.k5, rn.k1, rn.k3])
    @test collect(values(out)) == collect(values(out_merged))

    out_local = assess_local_identifiability(rn; measured_quantities = [:X], known_p = [:k5],
        funcs_to_check = [rn.k5, rn.k1, rn.k3], loglevel)
    out_local_merged = assess_local_identifiability(rn_merged; measured_quantities = [:Y],
        known_p = [:k4], funcs_to_check = [rn_merged.k4, rn_merged.k1, rn_merged.k3], loglevel)
    @test isequal(collect(keys(out_local)), [rn.k5, rn.k1, rn.k3])
    @test collect(values(out_local)) == collect(values(out_local_merged))

    # Both names of an aliased parameter get their result.
    out_both = assess_local_identifiability(rn; measured_quantities = [:X], known_p = [:k5],
        funcs_to_check = [rn.k5, rn.k4, rn.k1], loglevel)
    @test isequal(collect(keys(out_both)), [rn.k5, rn.k4, rn.k1])
    @test collect(values(out_both)) == collect(values(out_local))[[1, 1, 2]]

    # The identifiable functions are given in the surviving symbols, which are named as in the
    # hand-merged network.
    funcs = find_identifiable_functions(rn; measured_quantities = [:X], known_p = [:k5], loglevel)
    funcs_merged = find_identifiable_functions(rn_merged; measured_quantities = [:Y],
        known_p = [:k4], loglevel)
    @test issetequal(string.(funcs), string.(funcs_merged))
end
