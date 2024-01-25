# Returns the length of a expression tuple, or 1 if it is not an expression tuple (probably a  Symbol/Numerical).
function tup_leng(ex::ExprValues)
    (typeof(ex) == Expr && ex.head == :tuple) && (return length(ex.args))
    return 1
end

#Gets the ith element in a expression tuple, or returns the input itself if it is not an expression tuple (probably a  Symbol/Numerical).
function get_tup_arg(ex::ExprValues, i::Int)
    (tup_leng(ex) == 1) && (return ex)
    return ex.args[i]
end

# In variable/species/parameters on the forms like:
# X
# X = 1.0
# X, [metadata=true]
# X = 1.0, [metadata=true]
# X(t)
# X(t) = 1.0
# X(t), [metadata=true]
# X(t) = 1.0, [metadata=true]
# Finds the: Variable name (X), Independent variable name(s) ([t]), default value (2.0), and metadata (:([metadata=true])).
# If a field does not exist (e.g. independent variable in `X, [metadata=true]`), gives nothing.
# The independent variables are given as a vector (empty if none given).
# Does not support e.g. "X [metadata=true]" (when metadata does not have a comma before).
function find_varinfo_in_declaration(expr)
    # Case: X
    (expr isa Symbol) && (return expr, [], nothing, nothing)                          
    # Case: X(t)         
    (expr.head == :call) && (return expr.args[1], expr.args[2:end], nothing, nothing)                   
    if expr.head == :(=)
        # Case: X = 1.0
        (expr.args[1] isa Symbol) && (return expr.args[1], [], expr.args[2], nothing)  
        # Case: X(t) = 1.0        
        (expr.args[1].head == :call) && (return expr.args[1].args[1], expr.args[1].args[2:end], expr.args[2].args[1], nothing) 
    end
    if expr.head == :tuple
        # Case: X, [metadata=true]
        (expr.args[1] isa Symbol) && (return expr.args[1], [], nothing, expr.args[2])          
        # Case: X(t), [metadata=true]
        (expr.args[1].head == :call) && (return expr.args[1].args[1], expr.args[1].args[2:end], nothing, expr.args[2]) 
        if (expr.args[1].head == :(=)) 
            # Case: X = 1.0, [metadata=true]
            (expr.args[1].args[1] isa Symbol) && (return expr.args[1].args[1], [], expr.args[1].args[2], expr.args[2]) 
            # Case: X(t) = 1.0, [metadata=true]
            (expr.args[1].args[1].head == :call) && (return expr.args[1].args[1].args[1], expr.args[1].args[1].args[2:end], expr.args[1].args[2].args[1], expr.args[2]) 
        end
    end
    error("Unable to detect the variable declared in expression: $expr.")
end

# Converts an expression of the forms:
# X
# X = 1.0
# X, [metadata=true]
# X = 1.0, [metadata=true]
# To the form:
# X(t)
# X(t) = 1.0
# X(t), [metadata=true]
# X(t) = 1.0, [metadata=true]
# (In this example the independent variable :t was inserted).
# Here, the iv is a iv_expr, which can be anything, which is inserted
function insert_independent_variable(expr_in, iv_expr)
    # If expr is a symbol, just attach the iv. If not we have to create a new expr and mutate it. 
    # Because Symbols (a possible input) cannot be mutated, this function cannot mutate the input
    # (would have been easier if Expr input was guaranteed).
    (expr_in isa Symbol) && (return Expr(:call, expr_in, iv_expr))
    expr = deepcopy(expr_in)

    # Loops through possible cases.
    if expr.head == :(=) 
        # Case: :(X = 1.0)
        expr.args[1] = Expr(:call, expr.args[1], iv_expr)
    elseif expr.head == :tuple
        if expr.args[1] isa Symbol 
            # Case: :(X, [metadata=true])
            expr.args[1] = Expr(:call, expr.args[1], iv_expr)
        elseif (expr.args[1].head == :(=)) && (expr.args[1].args[1] isa Symbol)
            # Case: :(X = 1.0, [metadata=true])
            expr.args[1].args[1] = Expr(:call, expr.args[1].args[1], iv_expr)
        end
    end
    return expr
end


### Old Stuff ###

#This will be called whenever a function stored in funcdict is called.
# function replace_names(expr, old_names, new_names)
#     mapping = Dict(zip(old_names, new_names))
#     MacroTools.postwalk( x -> x in old_names ? x= mapping[x] : x, expr)
# end

# below are currently not used

# function splitplus!(ex)
#     dosplit = ex.head == :(=) && ex.args[2] isa Expr && ex.args[2].head == :call && ex.args[2].args[1] == :(+)
#     if dosplit
#       summands = ex.args[2].args[2:end]
#       ex.args[2] = foldl((x,y)->(:(($x + $y))), summands)
#     end
#     dosplit
# end

# #Creates an expression which can be evaluated to an actual function. Input is an array of expressions where each entry is a line in the function. Uses the array of expressions generated in either get_f or get_g.
# function make_func(func_expr::Vector{Expr},reactants::OrderedDict{Symbol,Int}, parameters::OrderedDict{Symbol,Int})
#     system = Expr(:block)
#     for func_line in deepcopy(func_expr)
#         ex = recursive_replace!(func_line, (reactants,:internal_var___u), (parameters, :internal_var___p))
#         splitplus!(ex)
#         push!(system.args,ex)
#     end
#     push!(system.args, :(nothing))
#     return :((internal_var___du,internal_var___u,internal_var___p,t) -> @inbounds $system)
# end

# """
# clean_subtractions(ex::Expr)
# Replace additions of negative terms with subtractions.
# This is a fairly stupid function, which is designed for a specific problem
# with reaction networks. It is neither recursive nor very general.
# Return :: cleaned out expression

# From Latexify.jl with permission:
# [see](https://github.com/JuliaDiffEq/Catalyst.jl/issues/89#issuecomment-462147882)
# """
# function clean_subtractions(ex::Expr)
#     ex.args[1] != :+ && return ex

#     term = ex.args[2]

#     ### Sort out the first term
#     if term isa Expr && length(term.args) >= 3 && term.args[1:2] == [:*, -1]
#         result = :(- *($(term.args[3:end]...)))
#     else
#         result = :($term)
#     end

#     ### Sort out the other terms
#     for term in ex.args[3:end]
#         if term isa Expr && length(term.args) >= 3 && term.args[1:2] == [:*, -1]
#             result = :($result - *($(term.args[3:end]...)))
#         else
#             result = :($result + $term)
#         end
#     end
#     return result
# end

# #Recursively traverses an expression and removes things like X^1, 1*X. Will not actually have any effect on the expression when used as a function, but will make it much easier to look at it for debugging, as well as if it is transformed to LaTeX code.
# function recursive_clean!(expr::ExprValues)
#     (expr == :no___noise___scaling) && (return 1)
#     (typeof(expr)!=Expr) && (return expr)
#     for i = 1:length(expr.args)
#         expr.args[i] = recursive_clean!(expr.args[i])
#     end
#     (expr.args[1] == :^) && (expr.args[3] == 1) && (return expr.args[2])
#     if expr.args[1] == :*
#         in(0,expr.args) && (return 0)
#         i = 1
#         while (i = i + 1) <= length(expr.args)
#              if (typeof(expr.args[i]) == Expr) && (expr.args[i].head == :call) && (expr.args[i].args[1] == :*)
#                  for arg in expr.args[i].args
#                      (arg != :*) && push!(expr.args, arg)
#                  end
#              end
#         end
#         for i = length(expr.args):-1:2
#             (typeof(expr.args[i]) == Expr) && (expr.args[i].head == :call) && (expr.args[i].args[1] == :*) && deleteat!(expr.args,i)
#             (expr.args[i] == 1) && deleteat!(expr.args,i)
#         end
#         (length(expr.args) == 2) && (return expr.args[2])                   # We have a multiplication of only one thing, return only that thing.
#         (length(expr.args) == 1) && (return 1)                              # We have only * and no real arguments.
#         (length(expr.args) == 3) && (expr.args[2] == -1) && return :(-$(expr.args[3]))
#         (length(expr.args) == 3) && (expr.args[3] == -1) && return :(-$(expr.args[2]))
#     end
#     if expr.head == :call
#         (expr.args[1] == :/) && (expr.args[3] == 1) && (return expr.args[2])
#         haskey(funcdict, expr.args[1]) && return funcdict[expr.args[1]](expr.args[2:end])
#         in(expr.args[1],hill_name) && return hill(expr)
#         in(expr.args[1],hillR_name) && return hillR(expr)
#         in(expr.args[1],mm_name) && return mm(expr)
#         in(expr.args[1],mmR_name) && return mmR(expr)
#         (expr.args[1] == :binomial) && (expr.args[3] == 1) && return expr.args[2]
#         #@isdefined($(expr.args[1])) || error("Function $(expr.args[1]) not defined.")
#     end
#     return expr
# end

# #Recursively traverses an expression and replace instances of variables and parameters with things that the DifferentialEquations packages simulation algorithms can understand. E.g., X --> u[1], kB1 --> p[1] etc.
# function recursive_replace!(expr::ExprValues, replace_requests::Tuple{OrderedDict{Symbol,Int},Symbol}...)
#     if typeof(expr) == Symbol
#         for rr in replace_requests
#             (haskey(rr[1],expr)) && (return :($(rr[2])[$(rr[1][expr])]))
#         end
#     elseif typeof(expr) == Expr
#         for i = 1:length(expr.args)
#             expr.args[i] = recursive_replace!(expr.args[i], replace_requests...)
#         end
#     end
#     return expr
# end

# #Recursively traverses an expression and replaces a symbol with another.
# function recursive_replace!(expr::ExprValues, replace_requests::Dict{Symbol,Symbol})
#     if typeof(expr) == Symbol
#         haskey(replace_requests,expr) && return replace_requests[expr]
#     elseif typeof(expr) == Expr
#         for i = 1:length(expr.args)
#             expr.args[i] = recursive_replace!(expr.args[i], replace_requests)
#         end
#     end
#     return expr
# end

# #Recursive Contains, checks whether an expression contains a certain symbol.
# function recursive_contains(s,ex)
#     (typeof(ex)!=Expr) && (return s==ex)
#     for arg in ex.args
#         recursive_contains(s,arg) && (return true)
#     end
#     return false
# end

# #Parses an expression, and returns a set with all symbols in the expression, which is also a part of the provided vector with symbols (syms).
# function recursive_content(ex,syms::Vector{Symbol},content::Vector{Symbol})
#     if ex isa Symbol
#         in(ex,syms) && push!(content,ex)
#     elseif ex isa Expr
#         foreach(arg -> recursive_content(arg,syms,content), ex.args)
#     end
#     return content
# end

# function recursive_content(ex,symsmap::OrderedDict{Symbol,Int},content::Vector{Symbol})
#     if ex isa Symbol
#         haskey(symsmap,ex) && push!(content,ex)
#     elseif ex isa Expr
#         foreach(arg -> recursive_content(arg,symsmap,content), ex.args)
#     end
#     return content
# end

# #Turns an array of expressions into an expression block with corresponding expressions.
# function expr_arr_to_block(exprs)
#     block = :(begin end)
#     foreach(expr -> push!(block.args, expr), exprs)
#     return block
#   end
