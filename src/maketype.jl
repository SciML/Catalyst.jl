function maketype(name, syms, scale_noise; params = Symbol[], 
                reactions = Vector{ReactionStruct}(undef, 0),
                syms_to_ints = OrderedDict{Symbol,Int}(),
                params_to_ints = OrderedDict{Symbol,Int}())
    typeex = :(mutable struct $name <: DiffEqBase.AbstractReactionNetwork
        syms::Vector{Symbol}
        scale_noise::Symbol
        params::Vector{Symbol}
        reactions::Vector{ReactionStruct}
        syms_to_ints::OrderedDict{Symbol,Int}
        params_to_ints::OrderedDict{Symbol,Int}
    end)

    # Make the default constructor
    constructorex = :($(name)(;
                    $(Expr(:kw, :syms, syms)),
                    $(Expr(:kw, :scale_noise, Meta.quot(scale_noise))),
                    $(Expr(:kw, :params, params)),
                    $(Expr(:kw, :reactions, reactions)),
                    $(Expr(:kw, :syms_to_ints, syms_to_ints)),
                    $(Expr(:kw, :params_to_ints, params_to_ints))) =
                    $(name)(syms, scale_noise, params, reactions, syms_to_ints, params_to_ints)) |> esc

    # Make the type instance using the default constructor
    typeex, constructorex
end

mutable struct ODEReactionNetwork{Q,R,S,T,V,W <: DiffEqBase.AbstractReactionNetwork} 
    f_expr::Q
    f::R
    f_func::S
    symjac::T
    f_symfuncs::V
    rn::W
end

function ODEReactionNetwork(rn::DiffEqBase.AbstractReactionNetwork)
    @unpack reactions, syms_to_ints, params_to_ints = rn
        
    f_expr     = get_f(reactions, syms_to_ints)
    f          = eval(make_func(f_expr, syms_to_ints, params_to_ints))
    f_rhs      = [element.args[2] for element in f_expr]
    symjac     = eval( Expr(:quote, calculate_jac(deepcopy(f_rhs), rn.syms)) )
    f_symfuncs = hcat([SymEngine.Basic(f) for f in f_rhs])

    ODEReactionNetwork(f_expr, f, f_rhs, symjac, f_symfuncs, rn)
end

mutable struct SDEReactionNetwork{Q,R,S,T,U <: ODEReactionNetwork}
    g_expr::Q
    g::R
    p_matrix::S
    g_func::T
    odern::U
end

function SDEReactionNetwork(rn::DiffEqBase.AbstractReactionNetwork; odern=nothing)
    @unpack reactions, syms_to_ints, params_to_ints, scale_noise = rn

    # first construct an ODE reaction network
    oderxn = (odern == nothing) ? ODEReactionNetwork(rn) : odern

    g_expr   = get_g(reactions, syms_to_ints, scale_noise)
    g        = eval(make_func(g_expr, syms_to_ints, params_to_ints))
    g_func   = [element.args[2] for element in g_expr]
    p_matrix = zeros(length(syms_to_ints), length(reactions))

    SDEReactionNetwork(g_expr, g, p_matrix, g_func, oderxn)
end


mutable struct JumpReactionNetwork{Q,R,S,T,U <: DiffEqBase.AbstractReactionNetwork}
    jump_rate_expr::Q
    jump_affect_expr::R
    jumps::S
    regular_jumps::T
    rn::U
end

function JumpReactionNetwork(rn, args...)
    @unpack reactions, syms_to_ints, params_to_ints = rn

    # parse the jumps
    (jump_rate_expr, jump_affect_expr, jumps, regular_jumps) = get_jumps(reactions, syms_to_ints, params_to_ints)

    JumpReactionNetwork(jump_rate_expr, jump_affect_expr, jumps, regular_jumps, rn)
end