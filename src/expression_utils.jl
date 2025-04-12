### Expression Escaping ###

# Function that handles variable interpolation.
function esc_dollars!(ex)
    # If we do not have an expression: recursion has finished and we return the input.
    (ex isa Expr) || (return ex)

    # If we have encountered an interpolation, perform the appropriate modification, else recur.
    if ex.head == :$
        return esc(:($(ex.args[1])))
    else
        for i in eachindex(ex.args)
            ex.args[i] = esc_dollars!(ex.args[i])
        end
    end
    ex
end

# Checks if an expression is an escaped expression (e.g. on the form `$(Expr(:escape, :Y))`)
function is_escaped_expr(expr)
    return Meta.isexpr(expr, :escape) && (length(expr.args) == 1)
end

### Parameters/Species/Variables Symbols Correctness Checking ###

# Throws an error when a forbidden symbol is used.
function forbidden_symbol_check(syms)
    used_forbidden_syms = intersect(forbidden_symbols_error, syms)
    isempty(used_forbidden_syms) ||
        error("The following symbol(s) are used as species or parameters: $used_forbidden_syms, this is not permitted.")
end


### Catalyst-specific Expressions Manipulation ###

# Many option inputs can be on a form `@option input` or `@option begin ... end`. In both these
# cases we want to retrieve the third argument in the option expression. Further more, we wish
# to throw an error if there is more inputs (suggesting e.g. multiple inputs on a single line).
# Note that there are only some options for which we wish to make this check.
function get_block_option(expr)
    (length(expr.args) < 3) &&
        error("The $(expr.args[1]) option's input was misformatted (full declaration: `$expr`). It seems that it has no inputs, whereas some input is expected.")
    (length(expr.args) > 3) &&
        error("The $(expr.args[1]) option's input was misformatted (full declaration: `$expr`). Potentially, it has multiple inputs on a single line, in which case these should be split across multiple lines using a `begin ... end` block.")
    return expr.args[3]
end

# Some options takes input on form that is either `@option ...` or `@option begin ... end`.
# This transforms input of the latter form to the former (with only one line in the `begin ... end` block)
function option_block_form(expr)
    (expr.head == :block) && return expr
    return Expr(:block, expr)
end

# In variable/species/parameters on the forms like (examples, 16 alternatives in total):
# X
# X = 1.0, [misc=5]
# X(t)[1:2]
# X(t) = 1.0, [misc=5]
# X(t)[1:2] = 1.0, [misc = 5]
# Finds the: Variable name (X), Independent variable name(s) ([t]), indexes (1:2),
# default value (1.0), and metadata (:([metadata=true])).
# Information that does not exist (e.g. independent variable in `X, [metadata=true]`), is
# returned as `nothing`.
# The independent variables are returned as a vector (empty if none given).
function find_varinfo_in_declaration(expr::ExprValues)
    # Initialises values. Step by step, reads one and scales it away from the expression
    metadata = default = idxs = ivs = nothing

    # Reads and removes metadata (e.g. `[misc = 5]` in `:(X(t)[1:2] = 1.0, [misc = 5])`).
    if Meta.isexpr(expr, :tuple)
        metadata = expr.args[2]
        expr = expr.args[1]
    end
    # Reads and removes metadata (e.g. `1.0` in `:(X(t)[1:2] = 1.0)`).
    if Meta.isexpr(expr, :(=))
        default = expr.args[2]
        expr = expr.args[1]
    end
    # Reads and removes indexes (e.g. `[1:2]` in `:(X(t)[1:2])`).
    if Meta.isexpr(expr, :ref)
        idxs = expr.args[2]
        expr = expr.args[1]
    end
    # Reads and removes independent variables (e.g. `t` in `:(X(t))`).
    if Meta.isexpr(expr, :call)
        ivs = expr.args[2:end]
        expr = expr.args[1]
    end

    (expr isa Symbol) ||
        error("Erroneous expression encountered in `find_varinfo_in_declaration` (got `$expr` after processing, this should be a symbol).")
    return (;ivs, idxs, default, metadata)
end

# Sometimes (when declared on a single line, e.g. `@variables X [misc = 4] Y [misc = 5]`)
# the symbolic variable are metadata never form a single expression, but are simply separate
# (but adjacent) expressions in a longer expression array. This dispatch extracts the same
# information, but from this case. The input is the vector of all expressions, and the index
# of the symbolic variable which information we wish to extract.
function find_varinfo_in_declaration(exprs::Vector, idx::Integer)
    if (idx != length(exprs)) && Meta.isexpr(exprs[idx + 1], :vect)
        expr, (ivs, idxs, default, metadata) = find_varinfo_in_declaration(exprs[idx])
        return expr => (ivs, idxs, default, exprs[idx + 1])
    end
    return find_varinfo_in_declaration(exprs[idx])
end

# Checks an expression that declares several symbolic variables (either using `@variables ...`
# or using `@variables begin ... end`). Returns a dictionary with each information, using the
# form used in find_varinfo_in_declaration.
function find_all_varinfo_in_declaration(exprs)
    # (1)) Removes the macro call (2) Handles the `@variables begin .. end` case
    # (3) Find all indexes with variables (4) Extract their information.
    exprs = exprs.args[3:end]
    _head_equisexpral(exprs[1], :block) && (exprs = exprs[1].args)
    var_idxs = findall(!Meta.isexpr(expr, :vect) for expr in exprs)
    return Dict([find_varinfo_in_declaration(exprs, var_idx) for var_idx in var_idxs])
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
