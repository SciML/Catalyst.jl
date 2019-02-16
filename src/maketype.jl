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
                  jac = nothing,
                  paramjac = nothing,
                  jac_prototype = nothing,
                  symjac=Matrix{Expr}(undef,0,0),
                  reactions=Vector{ReactionStruct}(undef,0),
                  syms_to_ints = OrderedDict{Symbol,Int}(),
                  params_to_ints = OrderedDict{Symbol,Int}(),
                  odefun = nothing,
                  sdefun = nothing,
                  make_polynomial = nothing,
                  fixed_concentrations = Dict{Symbol,Polynomial}(),
                  homotopy_continuation_template = nothing,
                  equilibratium_polynomial = nothing,
                  test_f = nothing,
                  is_polynomial_system = make_poly_system()
                  )

    typeex = :(mutable struct $name <: $(abstracttype)
        f::Union{Function,Nothing}
        f_func::Union{Vector{Expr},Nothing}
        f_symfuncs::Union{Matrix{SymEngine.Basic},Nothing}
        g::Union{Function,Nothing}
        g_func::Union{Vector{Any},Nothing}
        jumps::Union{Tuple{Vararg{DiffEqJump.AbstractJump}},Nothing}
        regular_jumps::Union{RegularJump,Nothing}
        jump_rate_expr::Union{Tuple{Any,Vararg{Any}},Nothing}
        jump_affect_expr::Union{Tuple{Vector{Expr},Vararg{Vector{Expr}}},Nothing}
        p_matrix::Union{Array{Float64,2},Nothing}
        syms::Vector{Symbol}
        params::Vector{Symbol}
        jac::Union{Function,Nothing}
        paramjac::Union{Function,Nothing}
        jac_prototype::Nothing
        symjac::Union{Matrix{Expr},Nothing}
        reactions::Vector{ReactionStruct}
        syms_to_ints::OrderedDict{Symbol,Int}
        params_to_ints::OrderedDict{Symbol,Int}
        scale_noise::Symbol
        odefun::Union{ODEFunction,Nothing}
        sdefun::Union{SDEFunction,Nothing}
        make_polynomial::Union{Function,Nothing}
        fixed_concentrations::Dict{Symbol,Polynomial}
        homotopy_continuation_template::Union{Tuple{Array{Complex{Float64},1},Array{Array{Complex{Float64},1},1}},Nothing}
        equilibratium_polynomial::Union{Vector,Nothing}
        test_f::Any
        is_polynomial_system::Bool
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
                $(Expr(:kw,:jac,jac)),
                $(Expr(:kw,:paramjac,paramjac)),
                $(Expr(:kw,:jac_prototype,jac_prototype)),
                $(Expr(:kw,:reactions,reactions)),
                $(Expr(:kw,:syms_to_ints, syms_to_ints)),
                $(Expr(:kw,:params_to_ints, params_to_ints)),
                $(Expr(:kw,:scale_noise, Meta.quot(scale_noise))),
                $(Expr(:kw,:odefun, odefun)),
                $(Expr(:kw,:sdefun, sdefun)),
                $(Expr(:kw,:make_polynomial, make_polynomial)),
                $(Expr(:kw,:fixed_concentrations, fixed_concentrations)),
                $(Expr(:kw,:homotopy_continuation_template, homotopy_continuation_template)),
                $(Expr(:kw,:equilibratium_polynomial, equilibratium_polynomial)),
                $(Expr(:kw,:test_f, test_f)),
                $(Expr(:kw,:is_polynomial_system, is_polynomial_system))) =
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
                        jac,
                        paramjac,
                        jac_prototype,
                        symjac,
                        reactions,
                        syms_to_ints,
                        params_to_ints,
                        scale_noise,
                        odefun,
                        sdefun,
                        make_polynomial,
                        fixed_concentrations,
                        homotopy_continuation_template,
                        equilibratium_polynomial,
                        test_f,
                        is_polynomial_system
                        )) |> esc

    # Make the type instance using the default constructor
    typeex,constructorex
end

# type function expressions
function gentypefun_exprs(name; esc_exprs=true, gen_inplace=true, gen_outofplace=true, gen_constructor=true)
    exprs = Vector{Expr}(undef,0)

    ## Overload the type so that it can act as a function.
    if gen_inplace
        overloadex = :(((f::$name))(du, u, p, t::Number) = (f.f(du, u, p, t); nothing))
        push!(exprs,overloadex)
    end

    ## Add a method which allocates the `du` and returns it instead of being inplace
    if gen_outofplace
        overloadex = :(((f::$name))(u,p,t::Number) = (du=similar(u); f(du,u,p,t); du))
        push!(exprs,overloadex)
    end

    # export type constructor
    if gen_constructor
        def_const_ex = :(($name)())
        push!(exprs,def_const_ex)
    end

    # escape expressions for macros
    if esc_exprs
        for i in eachindex(exprs)
            exprs[i] = exprs[i] |> esc
        end
    end

    exprs
end

function addodes!(rn::DiffEqBase.AbstractReactionNetwork; kwargs...)
    @unpack reactions, syms_to_ints, params_to_ints, syms = rn

    (f_expr, f, f_rhs, symjac, jac, paramjac, f_symfuncs) = genode_exprs(reactions, syms_to_ints, params_to_ints, syms; kwargs...)
    rn.f          = eval(f)
    rn.f_func     = f_rhs
    rn.jac        = eval(jac)
    rn.paramjac   = eval(paramjac)
    rn.symjac     = eval(symjac)
    rn.f_symfuncs = f_symfuncs
    rn.odefun     = ODEFunction(rn.f; jac=rn.jac, jac_prototype=nothing, paramjac=rn.paramjac, syms=rn.syms)

    # functor for evaluating f
    functor_exprs = gentypefun_exprs(typeof(rn), esc_exprs=false, gen_constructor=false)
    eval( expr_arr_to_block(functor_exprs) )

    nothing
end

function addsdes!(rn::DiffEqBase.AbstractReactionNetwork)
    @unpack reactions, syms_to_ints, params_to_ints, scale_noise = rn

    # first construct an ODE reaction network
    if rn.f == nothing
        addodes!(rn)
    end

    (g_expr, g, g_funcs, p_matrix) = gensde_exprs(reactions, syms_to_ints, params_to_ints, scale_noise)
    rn.g        = eval(g)
    rn.g_func   = g_funcs
    rn.p_matrix = p_matrix
    rn.sdefun   = SDEFunction(rn.f, rn.g; jac=rn.jac, jac_prototype=nothing, paramjac=rn.paramjac, syms=rn.syms)

    nothing
end

function addjumps!(rn::DiffEqBase.AbstractReactionNetwork;
                                    build_jumps=true,
                                    build_regular_jumps=true,
                                    minimal_jumps=false)

    @unpack reactions, syms_to_ints, params_to_ints = rn

    # parse the jumps
    (jump_rate_expr, jump_affect_expr, jumps, regular_jumps) = get_jumps(reactions,
                                                                    syms_to_ints,
                                                                    params_to_ints;
                                                                    minimal_jumps=minimal_jumps)

    rn.jump_rate_expr   = jump_rate_expr
    rn.jump_affect_expr = jump_affect_expr
    rn.jumps            = build_jumps ? eval(jumps) : nothing
    rn.regular_jumps    = build_regular_jumps ? eval(regular_jumps) : nothing

    nothing
end
