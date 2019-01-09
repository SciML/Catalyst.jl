# parse rhs for ODEs from the ReactionNetwork
function gen_odefun_inplace!(rn::DiffEqBase.AbstractReactionNetwork)
    @unpack reactions, syms_to_ints, params_to_ints = rn
    
    f_expr     = get_f(reactions, syms_to_ints)
    f          = make_func(f_expr, syms_to_ints, params_to_ints)
    f_rhs      = [element.args[2] for element in f_expr]
    symjac     = Expr(:quote, calculate_jac(deepcopy(f_rhs), rn.syms))
    f_symfuncs = hcat([SymEngine.Basic(f) for f in f_rhs])
    
    rn.properties[:f]          = eval(f)
    rn.properties[:f_func]     = f_rhs
    rn.properties[:symjac]     = eval(symjac)
    rn.properties[:f_symfuncs] = f_symfuncs
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
        rn.properties[:ODEFun] = ODEFunction(rn.properties[:f]; syms=rn.syms)
    end

    rn.properties[:ODEFun]
end

# parse noise function and datatype matrix for ReactionNetwork SDE models
function gen_noisefun!(rn::DiffEqBase.AbstractReactionNetwork)
    @unpack reactions, syms_to_ints, params_to_ints, scale_noise = rn

    g_expr   = get_g(reactions, syms_to_ints, scale_noise)
    g        = make_func(g_expr, syms_to_ints, params_to_ints)
    g_funcs  = [element.args[2] for element in g_expr]
    p_matrix = zeros(length(syms_to_ints), length(reactions))

    rn.properties[:g]        = eval(g)
    rn.properties[:p_matrix] = p_matrix
    rn.properties[:g_func]   = g_funcs
    nothing
end


# generate a SDEFunction for the ReactionNetwork
function get_sdefun!(rn::DiffEqBase.AbstractReactionNetwork)

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


# parse jump functions from the ReactionNetwork
function gen_jumpfun!(rn::DiffEqBase.AbstractReactionNetwork)
    @unpack reactions, syms_to_ints, params_to_ints = rn
    
    (jump_rate_expr, jump_affect_expr, jumps, regular_jumps) = get_jumps(reactions, syms_to_ints, params_to_ints)

    rn.properties[:jumps] = eval(jumps)
    rn.properties[:regular_jumps] = eval(regular_jumps)
    rn.properties[:jump_rate_expr] = eval(Meta.quot(jump_rate_expr))
    rn.properties[:jump_affect_expr] = eval(Meta.quot(jump_affect_expr)) 

    nothing
end