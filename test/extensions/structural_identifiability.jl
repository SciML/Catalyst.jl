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
    goodwind_oscillator_catalyst = @reaction_network begin
        (mmr(P,pₘ,1), dₘ), 0 <--> M
        (pₑ*M,dₑ), 0 <--> E
        (pₚ*E,dₚ), 0 <--> P
    end
    goodwind_oscillator_si = @ODEmodel(
        M'(t) = pₘ / (1 + P(t)) - dₘ*M(t), 
        E'(t) = -dₑ*E(t) + pₑ*M(t),
        P'(t) = -dₚ*P(t) + pₚ*E(t),
        y1(t) = M(t)
    )
    si_catalyst_ode = make_si_ode(goodwind_oscillator_catalyst; measured_quantities=[:M])
    
    gi_1 = assess_identifiability(goodwind_oscillator_catalyst; measured_quantities=[:M])
    gi_2 = assess_identifiability(goodwind_oscillator_si)
    gi_3 = assess_identifiability(si_catalyst_ode)
    sym_dict(gi_1) == sym_dict(gi_2) == sym_dict(gi_3)
    
    li_1 = assess_local_identifiability(goodwind_oscillator_catalyst; measured_quantities=[:M])
    li_2 = assess_local_identifiability(goodwind_oscillator_si)
    li_3 = assess_local_identifiability(si_catalyst_ode)
    sym_dict(li_1) == sym_dict(li_2) == sym_dict(li_3)
    
    ifs1 = find_identifiable_functions(goodwind_oscillator_catalyst; measured_quantities=[:M])
    ifs2 = find_identifiable_functions(goodwind_oscillator_si)
    ifs3 = find_identifiable_functions(si_catalyst_ode)
    length(ifs1) == length(ifs2) == length(ifs3)    
end

# Tests on a made-up reaction network with mix of identifiable and non-identifiable components.
# Tests for symbolics input.
# Tests using known_p argument.
let 
    rs_catalyst = @reaction_network begin
        (p1, d), 0 <--> X1
        k1, X1 --> X2
        (k2f,k2b), X2 <--> X3
        k3, X3 --> X4
        d, X4 --> 0
    end
    @unpack X2, X3 = rs_catalyst
    rs_si = @ODEmodel(
        X1'(t) = p1 - d*X1(t) - k1*X1(t), 
        X2'(t) = k1*X1(t) + k2b*X3(t) - k2f*X2(t),
        X3'(t) = -k2b*X3(t) + k2f*X2(t) - k3*X3(t), 
        X4'(t) = d*X4(t) + k3*X3(t),
        y1(t) = X2,
        y2(t) = X3,
        y3(t) = k2f
    )
    rs_ode = make_si_ode(rs_catalyst; measured_quantities=[X2, X3], known_p=[:k2f])
    
    known_quantities = make_measured_quantities(rs, measured_quantities, known_p)
    
    gi_1 = assess_identifiability(rs_catalyst; measured_quantities=[X2, X3], known_p=[:k2f])
    gi_2 = assess_identifiability(rs_si)
    gi_3 = assess_identifiability(rs_ode)
    sym_dict(gi_1) == sym_dict(gi_2) == sym_dict(gi_3)
    
    li_1 = assess_local_identifiability(rs_catalyst; measured_quantities=[X2, X3], known_p=[:k2f])
    li_2 = assess_local_identifiability(rs_si)
    li_3 = assess_local_identifiability(rs_ode)
    sym_dict(li_1) == sym_dict(li_2) == sym_dict(li_3)
    
    ifs1 = find_identifiable_functions(rs_catalyst; measured_quantities=[X2, X3], known_p=[:k2f])
    ifs2 = find_identifiable_functions(rs_si)
    ifs3 = find_identifiable_functions(rs_ode)
    length(ifs1) == length(ifs2) == length(ifs3)    
end

# Tests on a made-up reaction network with mix of identifiable and non-identifiable components.
# Tests for symbolics known_p
# Tests using an equation for measured quantity.
let 
    rs_catalyst = @reaction_network begin
        p, 0 --> X1
        k1, X1 --> X2
        k2, X2 --> X3
        k3, X3 --> X4
        k3, X3 --> X5
        d, (X4,X5) --> 0
        (kA*X3, kD), Yi <--> Ya
    end
    @unpack X1, X2, X3, X4, k1, k2, Yi, Ya = rs_catalyst
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
    rs_ode = make_si_ode(rs_catalyst; measured_quantities=[X + Yi, YaX2], known_p=[k1, kD])
    
    known_quantities = make_measured_quantities(rs, measured_quantities, known_p)
    
    gi_1 = assess_identifiability(rs_catalyst; measured_quantities=[X + Yi, YaX2], known_p=[k1, kD])
    gi_2 = assess_identifiability(rs_si)
    gi_3 = assess_identifiability(rs_ode)
    sym_dict(gi_1) == sym_dict(gi_2) == sym_dict(gi_3)
    
    li_1 = assess_local_identifiability(rs_catalyst; measured_quantities=[X + Yi, YaX2], known_p=[k1, kD])
    li_2 = assess_local_identifiability(rs_si)
    li_3 = assess_local_identifiability(rs_ode)
    sym_dict(li_1) == sym_dict(li_2) == sym_dict(li_3)
    
    ifs1 = find_identifiable_functions(rs_catalyst; measured_quantities=[X + Yi, YaX2], known_p=[k1, kD])
    ifs2 = find_identifiable_functions(rs_si)
    ifs3 = find_identifiable_functions(rs_ode)
    length(ifs1) == length(ifs2) == length(ifs3)    
end