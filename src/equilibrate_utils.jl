
#Generates the expression which is then converted into a function which generates polynomials for given parameters. Parameters not given can be given at a later stage, but parameters in the exponential must be given here.
function get_equilibration(params::Vector{Symbol},f_expr::Vector{Expr})
    func_expr = :((;TO___BE___REMOVED=to___be___removed) -> [])
    foreach(i -> push!(func_expr.args[1].args[1].args,Expr(:kw,$(params[i]),:(internal___polyvar___p[$i]))), 1:length(params))
    deleteat!(func_expr.args[1].args[1].args,1)
    foreach(poly->push!(func_expr.args[2].args[2].args,poly), f_expr)
    return func_expr
end
