### Custom CRN Function Implementations ###
"""
    mm(X,v,K) = v*X / (X + K)

A Michaelis-Menten rate function.
"""
mm(X, v, K) = v * X / (X + K)
@register_symbolic mm(X, v, K)
@register_derivative mm(X, v, K) 1 (v * K) / (X + K)^2
@register_derivative mm(X, v, K) 2 X / (X + K)
@register_derivative mm(X, v, K) 3 -v * X / (X + K)^2

# Registers the repressing Michaelis-Menten function.
"""
    mmr(X,v,K) = v*K / (X + K)

A repressive Michaelis-Menten rate function.
"""
mmr(X, v, K) = v * K / (X + K)
@register_symbolic mmr(X, v, K)
@register_derivative mmr(X, v, K) 1 -(v * K) / (X + K)^2
@register_derivative mmr(X, v, K) 2 K / (X + K)
@register_derivative mmr(X, v, K) 3 v * X / (X + K)^2

"""
    hill(X,v,K,n) = v*(X^n) / (X^n + K^n)

A Hill rate function.
"""
hill(X, v, K, n) = v * (X^n) / (X^n + K^n)
@register_symbolic hill(X, v, K, n)
@register_derivative hill(X, v, K, n) 1  v * n * (K^n) * (X^(n - 1)) / (X^n + K^n)^2
@register_derivative hill(X, v, K, n) 2 (X^n) / (X^n + K^n)
@register_derivative hill(X, v, K, n) 3 -v * n * (X^n) * (K^(n - 1)) / (X^n + K^n)^2
@register_derivative hill(X, v, K, n) 4 v * (X^n) * (K^n) * (log(X) - log(K)) / (X^n + K^n)^2

"""
    hillr(X,v,K,n) = v*(K^n) / (X^n + K^n)

A repressive Hill rate function.
"""
hillr(X, v, K, n) = v * (K^n) / (X^n + K^n)
@register_symbolic hillr(X, v, K, n)
@register_derivative hillr(X, v, K, n) 1  -v * n * (K^n) * (X^(n - 1)) / (X^n + K^n)^2
@register_derivative hillr(X, v, K, n) 2 (K^n) / (X^n + K^n)
@register_derivative hillr(X, v, K, n) 3 v * n * (X^n) * (K^(n - 1)) / (X^n + K^n)^2
@register_derivative hillr(X, v, K, n) 4 v * (X^n) * (K^n) * (log(K) - log(X)) / (X^n + K^n)^2

"""
    hillar(X,Y,v,K,n) = v*(X^n) / (X^n + Y^n + K^n)

An activation/repressing Hill rate function.
"""
hillar(X, Y, v, K, n) = v * (X^n) / (X^n + Y^n + K^n)
@register_symbolic hillar(X, Y, v, K, n)
@register_derivative hillar(X, Y, v, K, n) 1 v * n * (X^(n - 1)) * (Y^n + K^n) / (X^n + Y^n + K^n)^2
@register_derivative hillar(X, Y, v, K, n) 2 -v * n * (Y^(n - 1)) * (X^n) / (X^n + Y^n + K^n)^2
@register_derivative hillar(X, Y, v, K, n) 3 (X^n) / (X^n + Y^n + K^n)
@register_derivative hillar(X, Y, v, K, n) 4 -v * n * (K^(n - 1)) * (X^n) / (X^n + Y^n + K^n)^2
@register_derivative hillar(X, Y, v, K, n) 5 v * (X^n) * (log(X) * (Y^n + K^n) - (Y^n) * log(Y) - 
    (K^n) * log(K)) / (X^n + Y^n + K^n)^2


# Tuple storing all registered function (for use in various functionalities).
const registered_funcs = (mm, mmr, hill, hillr, hillar)

### Custom CRN FUnction-related Functions ###

"""
expand_registered_functions(in)

Takes an expression, and expands registered function expressions. E.g. `mm(X,v,K)` is replaced
with v*X/(X+K). Currently supported functions: `mm`, `mmr`, `hill`, `hillr`, and `hill`. Can
be applied to a reaction system, a reaction, an equation, or a symbolic expression. The input
is not modified, while an output with any functions expanded is returned. If applied to a 
reaction system model, any cached network properties are reset.
"""
function expand_registered_functions(expr)
    if hasnode(is_catalyst_function, expr)
        expr = replacenode(expr, expand_catalyst_function)
    end
    return expr
end

# Checks whether an expression corresponds to a catalyst function call (e.g. `mm(X,v,K)`).
function is_catalyst_function(expr)
    iscall(expr) || (return false)
    return operation(expr) in registered_funcs
end

# If the input expression corresponds to a catalyst function call (e.g. `mm(X,v,K)`), returns
# it in its expanded form. If not, returns the input expression.
function expand_catalyst_function(expr)
    is_catalyst_function(expr) || (return expr)
    args = sorted_arguments(expr)
    if operation(expr) == Catalyst.mm
        return args[2] * args[1] / (args[1] + args[3])
    elseif operation(expr) == Catalyst.mmr
        return args[2] * args[3] / (args[1] + args[3])
    elseif operation(expr) == Catalyst.hill
        return args[2] * (args[1]^args[4]) / ((args[1])^args[4] + (args[3])^args[4])
    elseif operation(expr) == Catalyst.hillr
        return args[2] * (args[3]^args[4]) / ((args[1])^args[4] + (args[3])^args[4])
    elseif operation(expr) == Catalyst.hillar
        return args[3] * (args[1]^args[5]) /
               ((args[1])^args[5] + (args[2])^args[5] + (args[4])^args[5])
    end
end

# If applied to a Reaction, return a reaction with its rate modified.
function expand_registered_functions(rx::Reaction)
    Reaction(expand_registered_functions(rx.rate), rx.substrates, rx.products,
        rx.substoich, rx.prodstoich, rx.netstoich, rx.only_use_rate, rx.metadata)
end

# If applied to a Equation, returns it with it applied to lhs and rhs.
function expand_registered_functions(eq::Equation)
    return expand_registered_functions(eq.lhs) ~ expand_registered_functions(eq.rhs)
end

# If applied to a continuous event, returns it applied to eqs and affect.
function expand_registered_functions(ce::MT.SymbolicContinuousCallback)
    eqs = expand_registered_functions(ce.eqs)
    affect = expand_registered_functions(ce.affect)
    return MT.SymbolicContinuousCallback(eqs, affect)
end

# If applied to a discrete event, returns it applied to condition and affects.
function expand_registered_functions(de::MT.SymbolicDiscreteCallback)
    condition = expand_registered_functions(de.condition)
    affects = expand_registered_functions(de.affects)
    return MT.SymbolicDiscreteCallback(condition, affects)
end

# If applied to a vector, applies it to every element in the vector.
function expand_registered_functions(vec::Vector)
    return [Catalyst.expand_registered_functions(element) for element in vec]
end

# If applied to a ReactionSystem, applied function to all Reactions and other Equations, and return updated system.
# Currently, `ModelingToolkitBase.has_X_events` returns `true` even if event vector is empty (hence
# this function cannot be used).
function expand_registered_functions(rs::ReactionSystem)
    @set! rs.eqs = Catalyst.expand_registered_functions(get_eqs(rs))
    @set! rs.rxs = Catalyst.expand_registered_functions(get_rxs(rs))
    if !isempty(MT.get_continuous_events(rs))
        @set! rs.continuous_events = Catalyst.expand_registered_functions(MT.get_continuous_events(rs))
    end
    if !isempty(MT.get_discrete_events(rs))
        @set! rs.discrete_events = Catalyst.expand_registered_functions(MT.get_discrete_events(rs))
    end
    reset_networkproperties!(rs)
    return rs
end
