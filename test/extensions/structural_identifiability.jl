### Prepares Tests ###

# Fetch packages.
using Catalyst, Logging, StructuralIdentifiability, Test

# Sets the `loglevel`.
loglevel = Logging.Error

# Helper function for checking that results are correct identifiability calls from different packages.
# Converts the output dicts from StructuralIdentifiability functions from "weird symbol => stuff" to "symbol => stuff" (the output have some strange meta data which prevents equality checks, this enables this).
# Structural identifiability also provides variables like x (rather than x(t)). This is a bug, but we have to convert to make it work (now just remove any (t) to make them all equal).
function sym_dict(dict_in)
    dict_out = Dict{Symbol, Any}()
    for key in keys(dict_in)
        sym_key = Symbol(key)
        sym_key = Symbol(replace(String(sym_key), "(t)" => ""))
        dict_out[sym_key] = dict_in[key]
    end
    return dict_out
end


### Basic Tests ###

# Tests for Goodwin model (model with both global, local, and non identifiable components).
# Tests for system using Catalyst function (in this case, Michaelis-Menten function)
let
    # Identifiability analysis for Catalyst model.
    goodwind_oscillator_catalyst = @reaction_network begin
        (mmr(P, pₘ, 1), dₘ), 0 <--> M
        (pₑ * M, dₑ), 0 <--> E
        (pₚ * E, dₚ), 0 <--> P
    end
    gi_1 = assess_identifiability(goodwind_oscillator_catalyst; measured_quantities = [:M], loglevel)
    li_1 = assess_local_identifiability(goodwind_oscillator_catalyst; measured_quantities = [:M], loglevel)
    ifs_1 = find_identifiable_functions(goodwind_oscillator_catalyst; measured_quantities = [:M], loglevel)

    # Identifiability analysis for Catalyst converted to StructuralIdentifiability.jl model.
    si_catalyst_ode = make_si_ode(goodwind_oscillator_catalyst; measured_quantities = [:M])
    gi_2 = assess_identifiability(si_catalyst_ode; loglevel)
    li_2 = assess_local_identifiability(si_catalyst_ode; loglevel)
    ifs_2 = find_identifiable_functions(si_catalyst_ode; loglevel)

    # Identifiability analysis for StructuralIdentifiability.jl model (declare this overwrites e.g. X2 variable etc.).
    goodwind_oscillator_si = @ODEmodel(
        M'(t) = pₘ / (1 + P(t)) - dₘ * M(t),
        E'(t) = -dₑ * E(t) + pₑ * M(t),
        P'(t) = -dₚ * P(t) + pₚ * E(t),
        y1(t) = M(t)
    )
    gi_3 = assess_identifiability(goodwind_oscillator_si; loglevel)
    li_3 = assess_local_identifiability(goodwind_oscillator_si; loglevel)
    ifs_3 = find_identifiable_functions(goodwind_oscillator_si; loglevel)

    # Check outputs.
    @test sym_dict(gi_1) == sym_dict(gi_2) == sym_dict(gi_3)
    @test sym_dict(li_1) == sym_dict(li_2) == sym_dict(li_3)
    @test length(ifs_1) == length(ifs_2) == length(ifs_3)

    # Checks output to manually checked correct answers.
    @test isequal(collect(keys(gi_1)), [unknowns(goodwind_oscillator_catalyst); parameters(goodwind_oscillator_catalyst)])
    @test isequal(collect(values(gi_1)), [:globally, :nonidentifiable, :globally, :globally, :globally, :nonidentifiable, :locally, :nonidentifiable, :locally])
    @test isequal(collect(keys(li_1)), [unknowns(goodwind_oscillator_catalyst); parameters(goodwind_oscillator_catalyst)])
    @test isequal(collect(values(li_1)), [1, 0, 1, 1, 1, 0, 1, 0, 1])
end

# Tests on a made-up reaction network with mix of identifiable and non-identifiable components.
# Tests for symbolics input.
# Tests using known_p argument.
let
    # Identifiability analysis for Catalyst model.
    rs_catalyst = @reaction_network begin
        (p1, d), 0 <--> X1
        k1, X1 --> X2
        (k2f, k2b), X2 <--> X3
        k3, X3 --> X4
        d, X4 --> 0
    end
    @unpack X2, X3 = rs_catalyst
    gi_1 = assess_identifiability(rs_catalyst; measured_quantities = [X2, X3], known_p = [:k2f], loglevel)
    li_1 = assess_local_identifiability(rs_catalyst; measured_quantities = [X2, X3], known_p = [:k2f], loglevel)
    ifs_1 = find_identifiable_functions(rs_catalyst; measured_quantities = [X2, X3], known_p = [:k2f], loglevel)

    # Identifiability analysis for Catalyst converted to StructuralIdentifiability.jl model.
    rs_ode = make_si_ode(rs_catalyst; measured_quantities = [X2, X3], known_p = [:k2f])
    gi_2 = assess_identifiability(rs_ode; loglevel)
    li_2 = assess_local_identifiability(rs_ode; loglevel)
    ifs_2 = find_identifiable_functions(rs_ode; loglevel)

    # Identifiability analysis for StructuralIdentifiability.jl model (declare this overwrites e.g. X2 variable etc.).
    rs_si = @ODEmodel(
        X1'(t) = p1 - d * X1(t) - k1 * X1(t),
        X2'(t) = k1 * X1(t) + k2b * X3(t) - k2f * X2(t),
        X3'(t) = -k2b * X3(t) + k2f * X2(t) - k3 * X3(t),
        X4'(t) = d * X4(t) + k3 * X3(t),
        y1(t) = X2,
        y2(t) = X3,
        y3(t) = k2f
    )
    gi_3 = assess_identifiability(rs_si; loglevel)
    li_3 = assess_local_identifiability(rs_si; loglevel)
    ifs_3 = find_identifiable_functions(rs_si; loglevel)

    # Check outputs.
    @test sym_dict(gi_1) == sym_dict(gi_2) == sym_dict(gi_3)
    @test sym_dict(li_1) == sym_dict(li_2) == sym_dict(li_3)
    @test length(ifs_1) == length(ifs_2) == length(ifs_3)

    # Checks output to manually checked correct answers.
    @test isequal(collect(keys(gi_1)), [unknowns(rs_catalyst); parameters(rs_catalyst)])
    @test isequal(collect(values(gi_1)), [:nonidentifiable, :globally, :globally, :nonidentifiable, :nonidentifiable, :nonidentifiable, :nonidentifiable, :globally, :globally, :globally])
    @test isequal(collect(keys(li_1)), [unknowns(rs_catalyst); parameters(rs_catalyst)])
    @test isequal(collect(values(li_1)), [0, 1, 1, 0, 0, 0, 0, 1, 1, 1])
end

# Tests on a made-up reaction network with mix of identifiable and non-identifiable components.
# Tests for system with conserved quantity.
# Tests for symbolics known_p
# Tests using an equation for measured quantity.
let
    # Identifiability analysis for Catalyst model.
    rs_catalyst = @reaction_network begin
        p, 0 --> X1
        k1, X1 --> X2
        k2, X2 --> X3
        k3, X3 --> X4
        k3, X3 --> X5
        d, (X4, X5) --> 0
        (kA * X3, kD), Yi <--> Ya
    end
    @unpack X1, X2, X3, X4, k1, k2, Yi, Ya, k1, kD = rs_catalyst
    gi_1 = assess_identifiability(rs_catalyst; measured_quantities = [X1 + Yi, Ya], known_p = [k1, kD], loglevel)
    li_1 = assess_local_identifiability(rs_catalyst; measured_quantities = [X1 + Yi, Ya], known_p = [k1, kD], loglevel)
    ifs_1 = find_identifiable_functions(rs_catalyst; measured_quantities = [X1 + Yi, Ya], known_p = [k1, kD], loglevel)

    # Identifiability analysis for Catalyst converted to StructuralIdentifiability.jl model.
    rs_ode = make_si_ode(rs_catalyst; measured_quantities = [X1 + Yi, Ya], known_p = [k1, kD], remove_conserved = false)
    gi_2 = assess_identifiability(rs_ode; loglevel)
    li_2 = assess_local_identifiability(rs_ode; loglevel)
    ifs_2 = find_identifiable_functions(rs_ode; loglevel)

    # Identifiability analysis for StructuralIdentifiability.jl model (declare this overwrites e.g. X2 variable etc.).
    rs_si = @ODEmodel(
        X1'(t) = p - k1 * X1(t),
        X2'(t) = k1 * X1(t) - k2 * X2(t),
        X3'(t) = k2 * X2(t) - 2k3 * X3(t),
        X4'(t) = -d * X4(t) + k3 * X3(t),
        X5'(t) = -d * X5(t) + k3 * X3(t),
        Yi'(t) = kD * Ya(t) - kA * Yi(t) * X3(t),
        Ya'(t) = -kD * Ya(t) + kA * Yi(t) * X3(t),
        y1(t) = X1 + Yi,
        y2(t) = Ya,
        y3(t) = k1,
        y4(t) = kD
    )
    gi_3 = assess_identifiability(rs_si; loglevel)
    li_3 = assess_local_identifiability(rs_si; loglevel)
    ifs_3 = find_identifiable_functions(rs_si; loglevel)

    # Check outputs.
    @test sym_dict(gi_1) == sym_dict(gi_2) == sym_dict(gi_3)
    @test sym_dict(li_1) == sym_dict(li_2) == sym_dict(li_3)
    @test length(ifs_1[2:end]) == length(ifs_2) == length(ifs_3) # In the first case, the conservation law parameter is also identifiable.

    # Checks output to manually checked correct answers.
    correct_gi = Pair.([unknowns(rs_catalyst); parameters(rs_catalyst)], [:globally, :locally, :locally, :nonidentifiable, :nonidentifiable, :globally, :globally, :globally, :globally, :locally, :locally, :nonidentifiable, :locally, :globally])
    correct_li = Pair.([unknowns(rs_catalyst); parameters(rs_catalyst)], [1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1])
    @test issetequal(gi_1, correct_gi)
    @test issetequal(li_1, correct_li)
end

# Tests that various inputs types work.
let
    goodwind_oscillator_catalyst = @reaction_network begin
        (mmr(P, pₘ, 1), dₘ), 0 <--> M
        (pₑ * M, dₑ), 0 <--> E
        (pₚ * E, dₚ), 0 <--> P
    end
    @unpack M, E, P, pₑ, pₚ, pₘ = goodwind_oscillator_catalyst
    si_catalyst_ode = make_si_ode(goodwind_oscillator_catalyst; measured_quantities = [:M])
    si_catalyst_ode = make_si_ode(goodwind_oscillator_catalyst; known_p = [:pₑ], ignore_no_measured_warn = true)
    si_catalyst_ode = make_si_ode(goodwind_oscillator_catalyst; measured_quantities = [:M], known_p = [:pₑ])
    si_catalyst_ode = make_si_ode(goodwind_oscillator_catalyst; measured_quantities = [:M, :E], known_p = [:pₑ])
    si_catalyst_ode = make_si_ode(goodwind_oscillator_catalyst; measured_quantities = [:M], known_p = [:pₑ, :pₚ])
    si_catalyst_ode = make_si_ode(goodwind_oscillator_catalyst; measured_quantities = [:M, :E], known_p = [:pₑ, :pₚ])
    si_catalyst_ode = make_si_ode(goodwind_oscillator_catalyst; measured_quantities = [M])
    si_catalyst_ode = make_si_ode(goodwind_oscillator_catalyst; known_p = [pₑ], ignore_no_measured_warn = true)
    si_catalyst_ode = make_si_ode(goodwind_oscillator_catalyst; measured_quantities = [M], known_p = [pₑ])
    si_catalyst_ode = make_si_ode(goodwind_oscillator_catalyst; measured_quantities = [M, E], known_p = [pₑ])
    si_catalyst_ode = make_si_ode(goodwind_oscillator_catalyst; measured_quantities = [M], known_p = [pₑ, pₚ])
    si_catalyst_ode = make_si_ode(goodwind_oscillator_catalyst; measured_quantities = [M, E], known_p = [pₑ, pₚ])
    si_catalyst_ode = make_si_ode(goodwind_oscillator_catalyst; measured_quantities = [M + pₑ])
    si_catalyst_ode = make_si_ode(goodwind_oscillator_catalyst; measured_quantities = [M + E, pₑ * M], known_p = [:pₑ])
    si_catalyst_ode = make_si_ode(goodwind_oscillator_catalyst; measured_quantities = [pₑ, pₚ], known_p = [pₑ])

    # Tests using model.component style (have to make system complete first).
    gw_osc_complt = complete(goodwind_oscillator_catalyst)
    @test make_si_ode(gw_osc_complt; measured_quantities = [gw_osc_complt.M]) isa ODE
    @test make_si_ode(gw_osc_complt; known_p = [gw_osc_complt.pₑ], ignore_no_measured_warn = true) isa ODE
    @test make_si_ode(gw_osc_complt; measured_quantities = [gw_osc_complt.M], known_p = [gw_osc_complt.pₑ]) isa ODE
    @test make_si_ode(gw_osc_complt; measured_quantities = [gw_osc_complt.M, gw_osc_complt.E], known_p = [gw_osc_complt.pₑ]) isa ODE
    @test make_si_ode(gw_osc_complt; measured_quantities = [gw_osc_complt.M], known_p = [gw_osc_complt.pₑ, gw_osc_complt.pₚ]) isa ODE
    @test make_si_ode(gw_osc_complt; measured_quantities = [gw_osc_complt.M], known_p = [:pₚ]) isa ODE
    @test make_si_ode(gw_osc_complt; measured_quantities = [gw_osc_complt.M * gw_osc_complt.E]) isa ODE
end

# Tests for hierarchical model with conservation laws at both top and internal levels.
let
    # Identifiability analysis for Catalyst model.
    rs1 = @network_component rs1 begin
        (k1, k2), X1 <--> X2
    end
    rs2 = @network_component rs2 begin
        (k3, k4), X3 <--> X4
    end
    @named rs_catalyst = compose(rs1, [rs2])
    rs_catalyst = complete(rs_catalyst)
    @unpack X1, X2, k1, k2 = rs1
    gi_1 = assess_identifiability(rs_catalyst; measured_quantities = [X1, X2, rs2.X3], known_p = [k1], loglevel)
    li_1 = assess_local_identifiability(rs_catalyst; measured_quantities = [X1, X2, rs2.X3], known_p = [k1], loglevel)
    ifs_1 = find_identifiable_functions(rs_catalyst; measured_quantities = [X1, X2, rs2.X3], known_p = [k1], loglevel)

    # Identifiability analysis for Catalyst converted to StructuralIdentifiability.jl model.
    rs_ode = make_si_ode(rs_catalyst; measured_quantities = [X1, X2, rs2.X3], known_p = [k1])
    gi_2 = assess_identifiability(rs_ode; loglevel)
    li_2 = assess_local_identifiability(rs_ode; loglevel)
    ifs_2 = find_identifiable_functions(rs_ode; loglevel)

    # Identifiability analysis for StructuralIdentifiability.jl model (declare this overwrites e.g. X2 variable etc.).
    rs_si = @ODEmodel(
        X1'(t) = -k1 * X1(t) + k2 * X2(t),
        X2'(t) = k1 * X1(t) - k2 * X2(t),
        rs2₊X3'(t) = -rs2₊k3 * rs2₊X3(t) + rs2₊k4 * rs2₊X4(t),
        rs2₊X4'(t) = rs2₊k3 * rs2₊X3(t) - rs2₊k4 * rs2₊X4(t),
        y1(t) = X1,
        y2(t) = X2,
        y3(t) = rs2₊X3,
        y4(t) = k1
    )
    gi_3 = assess_identifiability(rs_si; loglevel)
    li_3 = assess_local_identifiability(rs_si; loglevel)
    ifs_3 = find_identifiable_functions(rs_si; loglevel)

    # Check outputs.
    @test sym_dict(gi_1) == sym_dict(gi_3)
    @test sym_dict(li_1) == sym_dict(li_3)
    @test (length(ifs_1) - 2) == (length(ifs_2) - 2) == length(ifs_3) # In the first case, the conservation law parameter is also identifiable.

    # Checks output for the SI converted version of the catalyst model.
    # For nested systems with conservation laws, conserved quantities like Γ[1], cannot be replaced back.
    # Hence, here you display identifiability for `Γ[1]` instead of X2.
    gi_1_no_cq = filter(x -> !occursin("X2", String(x[1])) && !occursin("X4", String(x[1])), sym_dict(gi_1))
    gi_2_no_cq = filter(x -> !occursin("Γ", String(x[1])), sym_dict(gi_2))
    li_1_no_cq = filter(x -> !occursin("X2", String(x[1])) && !occursin("X4", String(x[1])), sym_dict(li_1))
    li_2_no_cq = filter(x -> !occursin("Γ", String(x[1])), sym_dict(li_2))
    @test gi_1_no_cq == gi_2_no_cq
    @test li_1_no_cq == li_2_no_cq
end

# Tests directly on reaction systems with known identifiability structures.
# Test provided by Alexander Demin.
let
    rs = @reaction_network begin
        k1, x1 --> x2
    end
    # Measure the source
    id_report = assess_identifiability(rs; measured_quantities = [:x1], loglevel)
    @test sym_dict(id_report) == Dict(
        :x1 => :globally,
        :x2 => :nonidentifiable,
        :k1 => :globally
    )
    # Measure the target instead
    id_report = assess_identifiability(rs; measured_quantities = [:x2], loglevel)
    @test sym_dict(id_report) == Dict(
        :x1 => :globally,
        :x2 => :globally,
        :k1 => :globally
    )

    # Example from
    #   Identifiability of chemical reaction networks
    #   DOI: 10.1007/s10910-007-9307-x
    # The rate constants a, b, c are not identifiable even if all of the species
    # are observed.
    rs = @reaction_network begin
        a, A0 --> 2A1
        b, A0 --> 2A2
        c, A0 --> A1 + A2
    end
    id_report = assess_identifiability(rs; measured_quantities = [:A0, :A1, :A2], loglevel)
    @test sym_dict(id_report) == Dict(
        :A0 => :globally,
        :A1 => :globally,
        :A2 => :globally,
        :a => :nonidentifiable,
        :b => :nonidentifiable,
        :c => :nonidentifiable
    )

    # Test with no parameters
    rs = @reaction_network begin
        1, x1 --> x2
        1, x2 --> x3
    end
    id_report = assess_identifiability(rs; measured_quantities = [:x3], loglevel)
    @test sym_dict(id_report) == Dict(
        :x1 => :globally,
        :x2 => :globally,
        :x3 => :globally,
    )
    @test length(find_identifiable_functions(rs; measured_quantities = [:x3], loglevel)) == 1
end


### Other Tests ###

# Checks that identifiability can be assessed for coupled CRN/DAE systems.
# `remove_conserved = false` is used to remove info print statement from log.
let
    rs = @reaction_network begin
        @parameters k c1 c2
        @variables C(t)
        @equations begin
            D(V) ~ k * X - V
            C ~ (c1 + c2) * X / V
        end
        (p / V, d / V), 0 <--> X
    end
    @unpack p, d, k, c1, c2 = rs

    # Tests identifiability assessment when all unknowns are measured.
    remove_conserved = false
    gi_1 = assess_identifiability(rs; measured_quantities = [:X, :V, :C], loglevel, remove_conserved)
    li_1 = assess_local_identifiability(rs; measured_quantities = [:X, :V, :C], loglevel, remove_conserved)
    ifs_1 = find_identifiable_functions(rs; measured_quantities = [:X, :V, :C], loglevel, remove_conserved)
    @test sym_dict(gi_1) == Dict(
        [
            :X => :globally, :C => :globally, :V => :globally, :k => :globally,
            :c1 => :nonidentifiable, :c2 => :nonidentifiable, :p => :globally, :d => :globally,
        ]
    )
    @test sym_dict(li_1) == Dict([:X => 1, :C => 1, :V => 1, :k => 1, :c1 => 0, :c2 => 0, :p => 1, :d => 1])
    @test issetequal(ifs_1, [d, p, k, c1 + c2])

    # Tests identifiability assessment when only variables are measured.
    # Checks that a parameter in an equation can be set as known.
    gi_2 = assess_identifiability(rs; measured_quantities = [:V, :C], known_p = [:c1], loglevel, remove_conserved)
    li_2 = assess_local_identifiability(rs; measured_quantities = [:V, :C], known_p = [:c1], loglevel, remove_conserved)
    ifs_2 = find_identifiable_functions(rs; measured_quantities = [:V, :C], known_p = [:c1], loglevel, remove_conserved)
    @test sym_dict(gi_2) == Dict(
        [
            :X => :nonidentifiable, :C => :globally, :V => :globally, :k => :nonidentifiable,
            :c1 => :globally, :c2 => :nonidentifiable, :p => :nonidentifiable, :d => :globally,
        ]
    )
    @test sym_dict(li_2) == Dict([:X => 0, :C => 1, :V => 1, :k => 0, :c1 => 1, :c2 => 0, :p => 0, :d => 1])
    @test issetequal(ifs_2, [d, c1, k * p, c1 * p + c2 * p])
end

# Checks that identifiability functions cannot be applied to non-complete `ReactionSystems`s.
let
    # Create model.
    incomplete_network = @network_component begin
        (p, d), 0 <--> X
    end
    measured_quantities = [:X]

    # Computes bifurcation diagram.
    @test_throws Exception assess_identifiability(incomplete_network; measured_quantities, loglevel)
    @test_throws Exception assess_local_identifiability(incomplete_network; measured_quantities, loglevel)
    @test_throws Exception find_identifiable_functions(incomplete_network; measured_quantities, loglevel)
end
