# generate an inplace rhs for ReactionNetwork ODE models
function gen_odefun_inplace(rn::DiffEqBase.AbstractReactionNetwork)
    @unpack reactions, syms_to_ints, params_to_ints = rn
    f_expr = get_f(reactions, syms_to_ints)
    eval(make_func(f_expr, syms_to_ints, params_to_ints))
end

# generate noise function and datatype matrix for ReactionNetwork SDE models
function gen_noisefun(rn::DiffEqBase.AbstractReactionNetwork)
    @unpack reactions, syms_to_ints, params_to_ints, scale_noise = rn
    g_expr = get_g(reactions, syms_to_ints, scale_noise)
    g = eval(make_func(g_expr, syms_to_ints, params_to_ints))
    p_matrix = zeros(length(syms_to_ints), length(reactions))
    
    g,p_matrix
end