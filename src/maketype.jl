function maketype(abstracttype, 
                  name,
                  f,
                  f_func,
                  f_symfuncs,
                  g,
                  g_func,
                  jumps,
                  regular_jumps,
                  jump_rate_expr,
                  jump_affect_expr,
                  p_matrix,
                  syms,
                  scale_noise;
                  params = Symbol[],
                  pfuncs=Vector{Expr}(undef,0),
                  symjac=Matrix{Expr}(undef,0,0),
                  reactions=Vector{ReactionStruct}(undef,0),
                  syms_to_ints = OrderedDict{Symbol,Int}(),
                  params_to_ints = OrderedDict{Symbol,Int}()              
                  )

    typeex = :(mutable struct $name <: $(abstracttype)
        f::Union{Function,Nothing}
        f_func::Union{Vector{Expr},Nothing}
        f_symfuncs::Union{Matrix{SymEngine.Basic},Nothing}
        g::Union{Function,Nothing}
        g_func::Union{Vector{Any},Nothing}
        jumps::Union{Tuple{DiffEqJump.AbstractJump,Vararg{DiffEqJump.AbstractJump}},Nothing}
        regular_jumps::Union{RegularJump,Nothing}
        jump_rate_expr::Union{Tuple{Any,Vararg{Any}},Nothing}
        jump_affect_expr::Union{Tuple{Vector{Expr},Vararg{Vector{Expr}}},Nothing}
        p_matrix::Union{Array{Float64,2},Nothing}
        syms::Vector{Symbol}
        params::Vector{Symbol}
        symjac::Union{Matrix{Expr},Nothing}
        reactions::Vector{ReactionStruct}
        syms_to_ints::OrderedDict{Symbol,Int}
        params_to_ints::OrderedDict{Symbol,Int}
        scale_noise::Symbol
    end)
    # Make the default constructor
    constructorex = :($(name)(;
                $(Expr(:kw,:f,f)),
                $(Expr(:kw,:f_func,f_func)),
                $(Expr(:kw,:g,g)),
                $(Expr(:kw,:g_func,g_func)),
                $(Expr(:kw,:jumps,jumps)),
                $(Expr(:kw,:regular_jumps,regular_jumps)),
                $(Expr(:kw,:jump_rate_expr,jump_rate_expr)),
                $(Expr(:kw,:jump_affect_expr,jump_affect_expr)),
                $(Expr(:kw,:p_matrix,p_matrix)),
                $(Expr(:kw,:f_symfuncs,f_symfuncs)),
                $(Expr(:kw,:syms,syms)),
                $(Expr(:kw,:params,params)),
                $(Expr(:kw,:symjac,symjac)),
                $(Expr(:kw,:reactions,reactions)),
                $(Expr(:kw,:syms_to_ints, syms_to_ints)),
                $(Expr(:kw,:params_to_ints, params_to_ints)), 
                $(Expr(:kw,:scale_noise, Meta.quot(scale_noise)))) =
                $(name)(
                        f,
                        f_func,
                        f_symfuncs,
                        g,
                        g_func,
                        jumps,
                        regular_jumps,
                        jump_rate_expr,
                        jump_affect_expr,
                        p_matrix,
                        syms,
                        params,
                        symjac,
                        reactions,
                        syms_to_ints,
                        params_to_ints,
                        scale_noise
                        )) |> esc

    # Make the type instance using the default constructor
    typeex,constructorex
end

function add_ode_funs!(rn::MinReactionNetwork)
    @unpack reactions, syms_to_ints, params_to_ints = rn

    f_expr        = get_f(reactions, syms_to_ints)
    rn.f          = eval(make_func(f_expr, syms_to_ints, params_to_ints))
    rn.f_func     = [element.args[2] for element in f_expr]
    rn.symjac     = eval( Expr(:quote, calculate_jac(deepcopy(rn.f_func), rn.syms)) )
    rn.f_symfuncs = hcat([SymEngine.Basic(f) for f in rn.f_func])

    nothing
end

function add_sde_funs!(rn::MinReactionNetwork)
    @unpack reactions, syms_to_ints, params_to_ints, scale_noise = rn

    # first construct an ODE reaction network
    if rn.f == nothing 
        gen_ode!(rn)
    end

    g_expr      = get_g(reactions, syms_to_ints, scale_noise)
    rn.g        = eval(make_func(g_expr, syms_to_ints, params_to_ints))
    rn.g_func   = [element.args[2] for element in g_expr]
    rn.p_matrix = zeros(length(syms_to_ints), length(reactions))

    nothing
end

function add_jump_funs!(rn::MinReactionNetwork)
    @unpack reactions, syms_to_ints, params_to_ints = rn

    # parse the jumps
    (jump_rate_expr, jump_affect_expr, jumps, regular_jumps) = get_jumps(reactions, syms_to_ints, params_to_ints)

    rn.jump_rate_expr   = jump_rate_expr
    rn.jump_affect_expr = jump_affect_expr
    rn.jumps            = eval(jumps)
    rn.regular_jumps    = eval(regular_jumps)

    nothing
end
