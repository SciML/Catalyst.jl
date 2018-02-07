function maketype(name,
                  f,
                  f_funcs,
                  f_symfuncs,
                  g,
                  g_funcs,
                  jumps,
                  jumps_expr,
                  p_matrix,
                  syms;
                  params = Symbol[],
                  pfuncs=Vector{Expr}(0),
                  symjac=Matrix{SymEngine.Basic}(0,0),
                  )

    typeex = :(mutable struct $name <: AbstractReaction
        f::Function
        f_funcs::Vector{Expr}
        f_symfuncs::Matrix{SymEngine.Basic}
        g::Function
        g_funcs::Expr
        jumps::Tuple{ConstantRateJump,Vararg{ConstantRateJump}}
        jumps_expr::Expr
        p_matrix::Array{Float64,2}
        syms::Vector{Symbol}
        params::Vector{Symbol}
        symjac::Matrix{SymEngine.Basic}
    end)
    # Make the default constructor
    constructorex = :($(name)(;
                  $(Expr(:kw,:f,f)),
                  $(Expr(:kw,:f_funcs,f_funcs)),
                  $(Expr(:kw,:g,g)),
                  $(Expr(:kw,:g_funcs,g_funcs)),
                  $(Expr(:kw,:jumps,jumps)),
                  $(Expr(:kw,:jumps_expr,jumps_expr)),
                  $(Expr(:kw,:p_matrix,p_matrix)),
                  $(Expr(:kw,:f_symfuncs,f_symfuncs)),
                  $(Expr(:kw,:syms,syms)),
                  $(Expr(:kw,:params,params)),
                  $(Expr(:kw,:symjac,symjac))) =
                  $(name)(
                      f,
                      f_funcs,
                      f_symfuncs,
                      g,
                      g_funcs,
                      jumps,
                      jumps_expr,
                      p_matrix,
                      syms,
                      params,
                      symjac,
                      )) |> esc

                      #f_funcs,symfuncs,pfuncs,syms,symjac) |> esc

    # Make the type instance using the default constructor
    typeex,constructorex
end
