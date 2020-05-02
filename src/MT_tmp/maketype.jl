function MT_maketype(abstracttype,
                  name,
                  reaction_system;
                  ode_system = nothing
                  # f,
                  # f_func,
                  # f_symfuncs,
                  # g,
                  # g_func,
                  # jumps,
                  # regular_jumps,
                  # jump_rate_expr,
                  # jump_affect_expr,
                  # p_matrix,
                  # syms,
                  # scale_noise;
                  # params = Symbol[],
                  # jac = nothing,
                  # paramjac = nothing,
                  # symjac=Matrix{ExprValues}(undef,0,0),
                  # reactions=Vector{ReactionStruct}(undef,0),
                  # syms_to_ints = OrderedDict{Symbol,Int}(),
                  # params_to_ints = OrderedDict{Symbol,Int}(),
                  # odefun = nothing,
                  # sdefun = nothing,
                  # equilibrate_content = nothing
                  )

    typeex = :(mutable struct $name <: $(abstracttype)
        reaction_system::ReactionSystem
        ode_system::Union{ODESystem,Nothing}
        # f::Union{Function,Nothing}
        # f_func::Union{Vector{ExprValues},Nothing}
        # f_symfuncs::Union{Matrix{SymEngine.Basic},Nothing}
        # g::Union{Function,Nothing}
        # g_func::Union{Vector{ExprValues},Nothing}
        # jumps::Union{Tuple{Vararg{DiffEqJump.AbstractJump}},Nothing}
        # regular_jumps::Union{RegularJump,Nothing}
        # jump_rate_expr::Union{Tuple{ExprValues,Vararg{ExprValues}},Nothing}
        # jump_affect_expr::Union{Tuple{Vector{Expr},Vararg{Vector{Expr}}},Nothing}
        # p_matrix::Union{Array{Float64,2},Nothing}
        # syms::Vector{Symbol}
        # params::Vector{Symbol}
        # jac::Union{Function,Nothing}
        # paramjac::Union{Function,Nothing}
        # symjac::Union{Matrix{ExprValues},Nothing}
        # reactions::Vector{ReactionStruct}
        # syms_to_ints::OrderedDict{Symbol,Int}
        # params_to_ints::OrderedDict{Symbol,Int}
        # scale_noise::Symbol
        # odefun::Union{ODEFunction,Nothing}
        # sdefun::Union{SDEFunction,Nothing}
        # equilibrate_content::Union{EquilibrateContent,Nothing}
    end)
    # Make the default constructor
    constructorex = :($(name)(;
                $(Expr(:kw,:reaction_system, reaction_system)),
                $(Expr(:kw,:ode_system,ode_system))) =
                # $(Expr(:kw,:f,f)),
                # $(Expr(:kw,:f_func,f_func)),
                # $(Expr(:kw,:g,g)),
                # $(Expr(:kw,:g_func,g_func)),
                # $(Expr(:kw,:jumps,jumps)),
                # $(Expr(:kw,:regular_jumps,regular_jumps)),
                # $(Expr(:kw,:jump_rate_expr,jump_rate_expr)),
                # $(Expr(:kw,:jump_affect_expr,jump_affect_expr)),
                # $(Expr(:kw,:p_matrix,p_matrix)),
                # $(Expr(:kw,:f_symfuncs,f_symfuncs)),
                # $(Expr(:kw,:syms,syms)),
                # $(Expr(:kw,:params,params)),
                # $(Expr(:kw,:symjac,symjac)),
                # $(Expr(:kw,:jac,jac)),
                # $(Expr(:kw,:paramjac,paramjac)),
                # $(Expr(:kw,:reactions,reactions)),
                # $(Expr(:kw,:syms_to_ints, syms_to_ints)),
                # $(Expr(:kw,:params_to_ints, params_to_ints)),
                # $(Expr(:kw,:scale_noise, Meta.quot(scale_noise))),
                # $(Expr(:kw,:odefun, odefun)),
                # $(Expr(:kw,:sdefun, sdefun)),
                # $(Expr(:kw,:equilibrate_content, equilibrate_content))) =
                $(name)(
                        reaction_system,
                        ode_system
                        # f,
                        # f_func,
                        # f_symfuncs,
                        # g,
                        # g_func,
                        # jumps,
                        # regular_jumps,
                        # jump_rate_expr,
                        # jump_affect_expr,
                        # p_matrix,
                        # syms,
                        # params,
                        # jac,
                        # paramjac,
                        # symjac,
                        # reactions,
                        # syms_to_ints,
                        # params_to_ints,
                        # scale_noise,
                        # odefun,
                        # sdefun,
                        # equilibrate_content
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
