### Fetch Packages ###

using Catalyst, Test
using StructuralIdentifiability


### Helper Function ###

# Converts the output dicts from StructuralIdentifiability functions from "weird symbol => stuff" to "symbol => stuff" (the output have some strange meta data which prevents equality checks, this enables this).
# Structural identifiability also provides variables like x (rather than x(t)). This is a bug, but we have to convert to make it work (now just remove any (t) to make them all equal).
function sym_dict(dict_in)
    dict_out = Dict{Symbol,Any}()
    for key in keys(dict_in)
        sym_key = Symbol(key)
        sym_key = Symbol(replace(String(sym_key), "(t)" => ""))
        dict_out[sym_key] = dict_in[key]
    end    
    return dict_out
end


### Run Tests ###

# Tests for Goodwin model (model with both global, local, and non identifiable components).
# Tests for system using Catalyst function (in this case, Michaelis-Menten function)
let 
    # Identifiability analysis for Catalyst model.
    goodwind_oscillator_catalyst = @reaction_network begin
        (mmr(P,pₘ,1), dₘ), 0 <--> M
        (pₑ*M,dₑ), 0 <--> E
        (pₚ*E,dₚ), 0 <--> P
    end
    gi_1 = assess_identifiability(goodwind_oscillator_catalyst; measured_quantities=[:M])
    li_1 = assess_local_identifiability(goodwind_oscillator_catalyst; measured_quantities=[:M])
    ifs_1 = find_identifiable_functions(goodwind_oscillator_catalyst; measured_quantities=[:M])

    # Identifiability analysis for Catalyst converted to StructuralIdentifiability.jl model.
    si_catalyst_ode = make_si_ode(goodwind_oscillator_catalyst; measured_quantities=[:M])
    gi_2 = assess_identifiability(si_catalyst_ode)
    li_2 = assess_local_identifiability(si_catalyst_ode)
    ifs_2 = find_identifiable_functions(si_catalyst_ode)

    # Identifiability analysis for StructuralIdentifiability.jl model (declare this overwrites e.g. X2 variable etc.).
    goodwind_oscillator_si = @ODEmodel(
        M'(t) = pₘ / (1 + P(t)) - dₘ*M(t), 
        E'(t) = -dₑ*E(t) + pₑ*M(t),
        P'(t) = -dₚ*P(t) + pₚ*E(t),
        y1(t) = M(t)
    )
    gi_3 = assess_identifiability(goodwind_oscillator_si)
    li_3 = assess_local_identifiability(goodwind_oscillator_si)
    ifs_3 = find_identifiable_functions(goodwind_oscillator_si)

    # Check outputs.
    @test sym_dict(gi_1) == sym_dict(gi_2) == sym_dict(gi_3)
    @test sym_dict(li_1) == sym_dict(li_2) == sym_dict(li_3)
    @test length(ifs_1) == length(ifs_2) == length(ifs_3)     

    # Checks output to manually checked correct answers.
    @test isequal(collect(keys(gi_1)), [states(goodwind_oscillator_catalyst); parameters(goodwind_oscillator_catalyst)])
    @test isequal(collect(values(gi_1)), [:globally, :nonidentifiable, :globally, :globally, :globally, :nonidentifiable, :locally, :nonidentifiable, :locally])
    @test isequal(collect(keys(li_1)), [states(goodwind_oscillator_catalyst); parameters(goodwind_oscillator_catalyst)])
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
        (k2f,k2b), X2 <--> X3
        k3, X3 --> X4
        d, X4 --> 0
    end
    @unpack X2, X3 = rs_catalyst    
    gi_1 = assess_identifiability(rs_catalyst; measured_quantities=[X2, X3], known_p=[:k2f])
    li_1 = assess_local_identifiability(rs_catalyst; measured_quantities=[X2, X3], known_p=[:k2f])
    ifs_1 = find_identifiable_functions(rs_catalyst; measured_quantities=[X2, X3], known_p=[:k2f])
    
    # Identifiability analysis for Catalyst converted to StructuralIdentifiability.jl model.
    rs_ode = make_si_ode(rs_catalyst; measured_quantities=[X2, X3], known_p=[:k2f])
    gi_2 = assess_identifiability(rs_ode)
    li_2 = assess_local_identifiability(rs_ode)
    ifs_2 = find_identifiable_functions(rs_ode)

    # Identifiability analysis for StructuralIdentifiability.jl model (declare this overwrites e.g. X2 variable etc.).
    rs_si = @ODEmodel(
        X1'(t) = p1 - d*X1(t) - k1*X1(t), 
        X2'(t) = k1*X1(t) + k2b*X3(t) - k2f*X2(t),
        X3'(t) = -k2b*X3(t) + k2f*X2(t) - k3*X3(t), 
        X4'(t) = d*X4(t) + k3*X3(t),
        y1(t) = X2,
        y2(t) = X3,
        y3(t) = k2f
    )
    gi_3 = assess_identifiability(rs_si)
    li_3 = assess_local_identifiability(rs_si)
    ifs_3 = find_identifiable_functions(rs_si)

    # Check outputs.
    @test sym_dict(gi_1) == sym_dict(gi_2) == sym_dict(gi_3)
    @test sym_dict(li_1) == sym_dict(li_2) == sym_dict(li_3)
    @test length(ifs_1) == length(ifs_2) == length(ifs_3)   

    # Checks output to manually checked correct answers. 
    @test isequal(collect(keys(gi_1)),[states(rs_catalyst); parameters(rs_catalyst)])
    @test isequal(collect(values(gi_1)),[:nonidentifiable, :globally, :globally, :nonidentifiable, :nonidentifiable, :nonidentifiable, :nonidentifiable, :globally, :globally, :globally])
    @test isequal(collect(keys(li_1)),[states(rs_catalyst); parameters(rs_catalyst)])
    @test isequal(collect(values(li_1)),[0, 1, 1, 0, 0, 0, 0, 1, 1, 1])  
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
        d, (X4,X5) --> 0
        (kA*X3, kD), Yi <--> Ya
    end
    @unpack X1, X2, X3, X4, k1, k2, Yi, Ya, k1, kD = rs_catalyst
    gi_1 = assess_identifiability(rs_catalyst; measured_quantities=[X1 + Yi, Ya], known_p=[k1, kD])
    li_1 = assess_local_identifiability(rs_catalyst; measured_quantities=[X1 + Yi, Ya], known_p=[k1, kD])
    ifs_1 = find_identifiable_functions(rs_catalyst; measured_quantities=[X1 + Yi, Ya], known_p=[k1, kD])

    # Identifiability analysis for Catalyst converted to StructuralIdentifiability.jl model.
    rs_ode = make_si_ode(rs_catalyst; measured_quantities=[X1 + Yi, Ya], known_p=[k1, kD], remove_conserved=false)
    gi_2 = assess_identifiability(rs_ode)
    li_2 = assess_local_identifiability(rs_ode)
    ifs_2 = find_identifiable_functions(rs_ode)

    # Identifiability analysis for StructuralIdentifiability.jl model (declare this overwrites e.g. X2 variable etc.).
    rs_si = @ODEmodel(
        X1'(t) = p - k1*X1(t),
        X2'(t) = k1*X1(t) - k2*X2(t),
        X3'(t) = k2*X2(t) - 2k3*X3(t),
        X4'(t) = -d*X4(t) + k3*X3(t),
        X5'(t) = -d*X5(t) + k3*X3(t),
        Yi'(t) = kD*Ya(t) - kA*Yi(t)*X3(t),
        Ya'(t) = -kD*Ya(t) + kA*Yi(t)*X3(t),
        y1(t) = X1 + Yi,
        y2(t) = Ya,
        y3(t) = k1,
        y4(t) = kD
    )
    gi_3 = assess_identifiability(rs_si)
    li_3 = assess_local_identifiability(rs_si)
    ifs_3 = find_identifiable_functions(rs_si)
    
    # Check outputs.
    @test sym_dict(gi_1) == sym_dict(gi_2) == sym_dict(gi_3)
    @test sym_dict(li_1) == sym_dict(li_2) == sym_dict(li_3)
    @test length(ifs_1[2:end]) == length(ifs_2) == length(ifs_3) # In the first case, the conservation law parameter is also identifiable.

    # Checks output to manually checked correct answers. 
    correct_gi = Pair.([states(rs_catalyst); parameters(rs_catalyst)], [:globally, :locally, :locally, :nonidentifiable, :nonidentifiable, :globally, :globally, :globally, :globally, :locally, :locally, :nonidentifiable, :locally, :globally])
    correct_li = Pair.([states(rs_catalyst); parameters(rs_catalyst)], [1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1])
    @test issetequal(gi_1, correct_gi)
    @test issetequal(li_1, correct_li)
end

# Tests that various inputs types work.
let 
    goodwind_oscillator_catalyst = @reaction_network begin
        (mmr(P,pₘ,1), dₘ), 0 <--> M
        (pₑ*M,dₑ), 0 <--> E
        (pₚ*E,dₚ), 0 <--> P
    end
    @unpack M, E, P, pₑ, pₚ, pₘ = goodwind_oscillator_catalyst
    si_catalyst_ode = make_si_ode(goodwind_oscillator_catalyst; measured_quantities=[:M])
    si_catalyst_ode = make_si_ode(goodwind_oscillator_catalyst; known_p=[:pₑ])
    si_catalyst_ode = make_si_ode(goodwind_oscillator_catalyst; measured_quantities=[:M], known_p=[:pₑ])
    si_catalyst_ode = make_si_ode(goodwind_oscillator_catalyst; measured_quantities=[:M, :E], known_p=[:pₑ])
    si_catalyst_ode = make_si_ode(goodwind_oscillator_catalyst; measured_quantities=[:M], known_p=[:pₑ, :pₚ])
    si_catalyst_ode = make_si_ode(goodwind_oscillator_catalyst; measured_quantities=[:M, :E], known_p=[:pₑ, :pₚ])
    si_catalyst_ode = make_si_ode(goodwind_oscillator_catalyst; measured_quantities=[M])
    si_catalyst_ode = make_si_ode(goodwind_oscillator_catalyst; known_p=[pₑ])
    si_catalyst_ode = make_si_ode(goodwind_oscillator_catalyst; measured_quantities=[M], known_p=[pₑ])
    si_catalyst_ode = make_si_ode(goodwind_oscillator_catalyst; measured_quantities=[M, E], known_p=[pₑ])
    si_catalyst_ode = make_si_ode(goodwind_oscillator_catalyst; measured_quantities=[M], known_p=[pₑ, pₚ])
    si_catalyst_ode = make_si_ode(goodwind_oscillator_catalyst; measured_quantities=[M, E], known_p=[pₑ, pₚ])
    si_catalyst_ode = make_si_ode(goodwind_oscillator_catalyst; measured_quantities=[M + pₑ])
    si_catalyst_ode = make_si_ode(goodwind_oscillator_catalyst; measured_quantities=[M + E, pₑ*M], known_p=[:pₑ])
    si_catalyst_ode = make_si_ode(goodwind_oscillator_catalyst; measured_quantities=[pₑ, pₚ], known_p=[pₑ])    

    # Tests using model.component style (have to make system complete first).
    gw_osc_complt = complete(goodwind_oscillator_catalyst)
    @test make_si_ode(gw_osc_complt; measured_quantities=[gw_osc_complt.M]) isa ODE
    @test make_si_ode(gw_osc_complt; known_p=[gw_osc_complt.pₑ]) isa ODE
    @test make_si_ode(gw_osc_complt; measured_quantities=[gw_osc_complt.M], known_p=[gw_osc_complt.pₑ]) isa ODE
    @test make_si_ode(gw_osc_complt; measured_quantities=[gw_osc_complt.M, gw_osc_complt.E], known_p=[gw_osc_complt.pₑ]) isa ODE
    @test make_si_ode(gw_osc_complt; measured_quantities=[gw_osc_complt.M], known_p=[gw_osc_complt.pₑ, gw_osc_complt.pₚ]) isa ODE
    @test make_si_ode(gw_osc_complt; measured_quantities=[gw_osc_complt.M], known_p = [:pₚ]) isa ODE
    @test make_si_ode(gw_osc_complt; measured_quantities=[gw_osc_complt.M*gw_osc_complt.E]) isa ODE
end

# Check that `prob_threshold` alternative kwarg works.
let 
    rs = @reaction_network begin
        p, X --> 0
    end

    assess_identifiability(rs_catalyst; measured_quantities=[rs.X], prob_thres=0.9)
    assess_identifiability(rs_catalyst; measured_quantities=[rs.X], prob_thres=0.999)
end

# Tests for hierarchical model with conservation laws at both top and internal levels.
let
    # Identifiability analysis for Catalyst model.
    rs1 = @reaction_network rs1 begin
        (k1, k2), X1 <--> X2
    end
    rs2 = @reaction_network rs2 begin
        (k3, k4), X3 <--> X4
    end
    @named rs_catalyst = compose(rs1, [rs2])
    @unpack X1, X2, k1, k2 = rs1
    gi_1 = assess_identifiability(rs_catalyst; measured_quantities=[X1, X2, rs2.X3], known_p=[k1])
    li_1 = assess_local_identifiability(rs_catalyst; measured_quantities=[X1, X2, rs2.X3], known_p=[k1])
    ifs_1 = find_identifiable_functions(rs_catalyst; measured_quantities=[X1, X2, rs2.X3], known_p=[k1])

    # Identifiability analysis for Catalyst converted to StructuralIdentifiability.jl model.
    rs_ode = make_si_ode(rs_catalyst; measured_quantities=[X1, X2, rs2.X3], known_p=[k1])
    gi_2 = assess_identifiability(rs_ode)
    li_2 = assess_local_identifiability(rs_ode)
    ifs_2 = find_identifiable_functions(rs_ode)

    # Identifiability analysis for StructuralIdentifiability.jl model (declare this overwrites e.g. X2 variable etc.).
    rs_si = @ODEmodel(
        X1'(t) = -k1*X1(t) + k2*X2(t), 
        X2'(t) = k1*X1(t) - k2*X2(t),
        rs2₊X3'(t) = -rs2₊k3*rs2₊X3(t) + rs2₊k4*rs2₊X4(t), 
        rs2₊X4'(t) = rs2₊k3*rs2₊X3(t) - rs2₊k4*rs2₊X4(t),
        y1(t) = X1,
        y2(t) = X2,
        y3(t) = rs2₊X3,
        y4(t) = k1
    )
    gi_3 = assess_identifiability(rs_si)
    li_3 = assess_local_identifiability(rs_si)
    ifs_3 = find_identifiable_functions(rs_si)
        
    # Check outputs.
    @test sym_dict(gi_1) == sym_dict(gi_3)
    @test sym_dict(li_1) == sym_dict(li_3)
    @test length(ifs_1)-2 == length(ifs_2)-2 == length(ifs_3) # In the first case, the conservation law parameter is also identifiable.

    # Checks output for the SI converted version of the catalyst model.
    # For nested systems with conservation laws, conserved quantities like Γ[1], cannot be replaced back.
    # Hence, here you display identifiability for `Γ[1]` instead of X2.
    gi_1_no_cq = filter(x -> !occursin("X2",String(x[1])) && !occursin("X4",String(x[1])), sym_dict(gi_1))
    gi_2_no_cq = filter(x -> !occursin("Γ",String(x[1])), sym_dict(gi_2))
    li_1_no_cq = filter(x -> !occursin("X2",String(x[1])) && !occursin("X4",String(x[1])), sym_dict(li_1))
    li_2_no_cq = filter(x -> !occursin("Γ",String(x[1])), sym_dict(li_2))
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
    id_report = assess_identifiability(rs, measured_quantities = [:x1])
    @test sym_dict(id_report) == Dict(
        :x1 => :globally,
        :x2 => :nonidentifiable,
        :k1 => :globally
    )
    # Measure the target instead
    id_report = assess_identifiability(rs, measured_quantities = [:x2])
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
    id_report = assess_identifiability(rs, measured_quantities = [:A0, :A1, :A2])
    @test sym_dict(id_report) == Dict(
        :A0 => :globally,
        :A1 => :globally,
        :A2 => :globally,
        :a  => :nonidentifiable,
        :b  => :nonidentifiable,
        :c  => :nonidentifiable
    )

    # Test with no parameters
    rs = @reaction_network begin
        1, x1 --> x2
        1, x2 --> x3
    end
    id_report = assess_identifiability(rs, measured_quantities = [:x3])
    @test sym_dict(id_report) == Dict(
        :x1 => :globally,
        :x2 => :globally,
        :x3 => :globally,
    )
    @test length(find_identifiable_functions(rs, measured_quantities = [:x3])) == 1
end