# parse rhs for ODEs from the ReactionNetwork
function gen_odefun_inplace!(rn::DiffEqBase.AbstractReactionNetwork)
    @unpack reactions, syms_to_ints, params_to_ints = rn
    
    f_expr                 = get_f(reactions, syms_to_ints)
    rn.properties[:f]      = eval(make_func(f_expr, syms_to_ints, params_to_ints))
    rn.properties[:f_func] = [eval(element.args[2]) for element in f_expr]    
    rn.properties[:symjac] = eval(Expr(:quote, calculate_jac(deepcopy(f_rhs), syms))
    
    nothing
end

# generate an ODEFunction for the ReactionNetwork
function get_odefun!(rn::DiffEqBase.AbstractReactionNetwork)

    # if have already defined an ODEFun, just return it
    if !haskey(rn.properties, :ODEFun)

        # check if already parsed rhs
        if !haskey(rn.properties, :f) 
            gen_odefun_inplace!(rn)        
        end
        rn.properties[:ODEFun] = ODEFunction(r.properties[:f]; syms=rn.syms)
    end

    rn.properties[:ODEFun]
end

# parse noise function and datatype matrix for ReactionNetwork SDE models
function gen_noisefun!(rn::DiffEqBase.AbstractReactionNetwork)
    @unpack reactions, syms_to_ints, params_to_ints, scale_noise = rn

    g_expr                   = get_g(reactions, syms_to_ints, scale_noise)
    rn.properties[:g]        = eval(make_func(g_expr, syms_to_ints, params_to_ints))
    rn.properties[:p_matrix] = zeros(length(syms_to_ints), length(reactions))

    nothing
end


# generate a SDEFunction for the ReactionNetwork
function gen_sdefun!(rn::DiffEqBase.AbstractReactionNetwork)

    # check that the ODE rhs has been parsed
    if !haskey(rn.properties, :f) 
        gen_odefun_inplace!(rn)
    end

    # now check for noise function
    if !haskey(rn.properties, :SDEFun) 
        if !haskey(rn.properties, :g)
            gen_noisefun!(rn)
        end
        rn.properties[:SDEFun] = SDEFunction(rn.properties[:f], rn.properties[:g]; syms=rn.syms)
    end

    rn.properties[:SDEFun]
end
