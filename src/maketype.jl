function maketype(name,
                  f,
                  f_expr,
                  f_symfuncs,
                  g,
                  g_expr,
                  jumps,
                  jump_rate_expr,
                  jump_affect_expr,
                  p_matrix,
                  syms;
                  params = Symbol[],
                  pfuncs=Vector{Expr}(0),
                  symjac=Matrix{SymEngine.Basic}(0,0),
                  )

    typeex = :(mutable struct $name <: AbstractReaction
        f::Function
        f_expr::Vector{Expr}
        f_symfuncs::Matrix{SymEngine.Basic}
        g::Function
        g_expr::Expr
        jumps::Tuple{ConstantRateJump,Vararg{ConstantRateJump}}
        jump_rate_expr::Tuple{Expr,Vararg{Expr}}
        jump_affect_expr::Tuple{Vector{Expr},Vararg{Vector{Expr}}}
        p_matrix::Array{Float64,2}
        syms::Vector{Symbol}
        params::Vector{Symbol}
        symjac::Matrix{SymEngine.Basic}
    end)
    # Make the default constructor
    constructorex = :($(name)(;
                  $(Expr(:kw,:f,f)),
                  $(Expr(:kw,:f_expr,f_expr)),
                  $(Expr(:kw,:g,g)),
                  $(Expr(:kw,:g_expr,g_expr)),
                  $(Expr(:kw,:jumps,jumps)),
                  $(Expr(:kw,:jump_rate_expr,jump_affect_expr)),
                  $(Expr(:kw,:jump_affect_expr,jump_affect_expr)),
                  $(Expr(:kw,:p_matrix,p_matrix)),
                  $(Expr(:kw,:f_symfuncs,f_symfuncs)),
                  $(Expr(:kw,:syms,syms)),
                  $(Expr(:kw,:params,params)),
                  $(Expr(:kw,:symjac,symjac))) =
                  $(name)(
                      f,
                      f_expr,
                      f_symfuncs,
                      g,
                      g_expr,
                      jumps,
                      jump_affect_expr,
                      jump_affect_expr,
                      p_matrix,
                      syms,
                      params,
                      symjac,
                      )) |> esc

                      #f_funcs,symfuncs,pfuncs,syms,symjac) |> esc

    # Make the type instance using the default constructor
    typeex,constructorex
end
