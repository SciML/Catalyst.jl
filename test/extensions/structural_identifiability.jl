### Fetch Packages ###

using Catalyst, Test
using StructuralIdentifiability


### Helper Function ###

# Converts the output dicts from StructuralIdentifiability functions from "weird symbol => stuff" to "symbol => stuff" (the output have some strange meta data which prevents equality checks, this enables this).
function sym_dict(dict_in)
    dict_out = Dict{Symbol,Any}()
    for key in keys(dict_in)
        dict_out[Symbol(key)] = dict_in[key]
    end
    return dict_out
end


### Run Tests ###

# Tests for Goodwin model (model with both global, local, and non identifiable components).
# Tests for system using Catalyst function (in this case, michaelis-menten function)
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
    sym_dict(gi_1) == sym_dict(gi_2) == sym_dict(gi_3)
    sym_dict(li_1) == sym_dict(li_2) == sym_dict(li_3)
    length(ifs_1) == length(ifs_2) == length(ifs_3)     
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
    sym_dict(gi_1) == sym_dict(gi_2) == sym_dict(gi_3)
    sym_dict(li_1) == sym_dict(li_2) == sym_dict(li_3)
    length(ifs_1) == length(ifs_2) == length(ifs_3)     
end

# Tests on a made-up reaction network with mix of identifiable and non-identifiable components.
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
    rs_ode = make_si_ode(rs_catalyst; measured_quantities=[X1 + Yi, Ya], known_p=[k1, kD])
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
        y1(t) = X1 + Ya,
        y2(t) = X2,
        y3(t) = k1,
        y4(t) = k2
    )
    gi_3 = assess_identifiability(rs_si)
    li_3 = assess_local_identifiability(rs_si)
    ifs_3 = find_identifiable_functions(rs_si)
    
    # Check outputs.
    sym_dict(gi_1) == sym_dict(gi_2) == sym_dict(gi_3)
    sym_dict(li_1) == sym_dict(li_2) == sym_dict(li_3)
    length(ifs_1) == length(ifs_2) == length(ifs_3)    
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
end