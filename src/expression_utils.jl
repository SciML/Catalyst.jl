### Expression Escaping ###

# Function that handles variable interpolation.
function esc_dollars!(ex)
    if ex isa Expr
        if ex.head == :$
            return esc(:($(ex.args[1])))
        else
            for i in 1:length(ex.args)
                ex.args[i] = esc_dollars!(ex.args[i])
            end
        end
    end
    ex
end

# Checks if an expression is an escaped expression (e.g. on the form `$(Expr(:escape, :Y))`)
function is_escaped_expr(expr)
    return (expr isa Expr) && (expr.head == :escape) && (length(expr.args) == 1)
end


### Generic Expression Manipulation ###

# Recursively traverses an expression and replaces special function call like "hill(...)" with the actual corresponding expression.
function recursive_expand_functions!(expr::ExprValues)
    (typeof(expr) != Expr) && (return expr)
    foreach(i -> expr.args[i] = recursive_expand_functions!(expr.args[i]),
            1:length(expr.args))
    if expr.head == :call
        !isdefined(Catalyst, expr.args[1]) && (expr.args[1] = esc(expr.args[1]))
    end
    expr
end

# Returns the length of a expression tuple, or 1 if it is not an expression tuple (probably a Symbol/Numerical).
function tup_leng(ex::ExprValues)
    (typeof(ex) == Expr && ex.head == :tuple) && (return length(ex.args))
    return 1
end

# Gets the ith element in a expression tuple, or returns the input itself if it is not an expression tuple
# (probably a  Symbol/Numerical).
function get_tup_arg(ex::ExprValues, i::Int)
    (tup_leng(ex) == 1) && (return ex)
    return ex.args[i]
end


### Parameters/Species/Variables Symbols Correctness Checking ###

# Throws an error when a forbidden symbol is used.
function forbidden_symbol_check(v)
    !isempty(intersect(forbidden_symbols_error, v)) &&
        error("The following symbol(s) are used as species or parameters: " *
              ((map(s -> "'" * string(s) * "', ",
                    intersect(forbidden_symbols_error, v))...)) *
              "this is not permited.")
    nothing
end

# Throws an error when a forbidden variable is used (a forbidden symbol that is not `:t`).
function forbidden_variable_check(v)
    !isempty(intersect(forbidden_variables_error, v)) &&
        error("The following symbol(s) are used as variables: " *
              ((map(s -> "'" * string(s) * "', ",
                    intersect(forbidden_variables_error, v))...)) *
              "this is not permited.")
end

function unique_symbol_check(syms)
    allunique(syms) ||
        error("Reaction network independent variables, parameters, species, and variables must all have distinct names, but a duplicate has been detected. ")
    nothing
end


### Catalyst-specific Expressions Manipulation ###

# Some options takes input on form that is either `@option ...` or `@option begin ... end`. 
# This transforms input of the latter form to the former (with only one line in the `begin ... end` block)
function option_block_form(expr)
    (expr.head == :block) && return expr
    return Expr(:block, expr)
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
    # Handles the $(Expr(:escape, :Y)) case:
    is_escaped_expr(expr) && (return find_varinfo_in_declaration(expr.args[1]))

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
