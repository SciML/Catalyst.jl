### Custom CRN Function Implementations ###
"""
    mm(X,v,K) = v*X / (X + K)

A Michaelis-Menten rate function.
"""
mm(X, v, K) = v * X / (X + K)
@register_symbolic mm(X, v, K);
function Symbolics.derivative(::typeof(mm), args::NTuple{3, Any}, ::Val{1})
    return (args[2] * args[3]) / (args[1] + args[3])^2
end
function Symbolics.derivative(::typeof(mm), args::NTuple{3, Any}, ::Val{2})
    return args[1] / (args[1] + args[3])
end
function Symbolics.derivative(::typeof(mm), args::NTuple{3, Any}, ::Val{3})
    return -args[2] * args[1] / (args[1] + args[3])^2
end

# Registers the repressing Michaelis-Menten function.
"""
    mmr(X,v,K) = v*K / (X + K)

A repressive Michaelis-Menten rate function.
"""
mmr(X, v, K) = v * K / (X + K)
@register_symbolic mmr(X, v, K);
function Symbolics.derivative(::typeof(mmr), args::NTuple{3, Any}, ::Val{1})
    return -(args[2] * args[3]) / (args[1] + args[3])^2
end
function Symbolics.derivative(::typeof(mmr), args::NTuple{3, Any}, ::Val{2})
    return args[3] / (args[1] + args[3])
end
function Symbolics.derivative(::typeof(mmr), args::NTuple{3, Any}, ::Val{3})
    return args[2] * args[1] / (args[1] + args[3])^2
end

"""
    hill(X,v,K,n) = v*(X^n) / (X^n + K^n)

A Hill rate function.
"""
hill(X, v, K, n) = v * (X^n) / (X^n + K^n)
@register_symbolic hill(X, v, K, n);
function Symbolics.derivative(::typeof(hill), args::NTuple{4, Any}, ::Val{1})
    return args[2] * args[4] * (args[3]^args[4]) * (args[1]^(args[4] - 1)) /
        (args[1]^args[4] + args[3]^args[4])^2
end
function Symbolics.derivative(::typeof(hill), args::NTuple{4, Any}, ::Val{2})
    return (args[1]^args[4]) / (args[1]^args[4] + args[3]^args[4])
end
function Symbolics.derivative(::typeof(hill), args::NTuple{4, Any}, ::Val{3})
    return -args[2] * args[4] * (args[1]^args[4]) * (args[3]^(args[4] - 1)) /
        (args[1]^args[4] + args[3]^args[4])^2
end
function Symbolics.derivative(::typeof(hill), args::NTuple{4, Any}, ::Val{4})
    return args[2] * (args[1]^args[4]) * (args[3]^args[4]) * (log(args[1]) - log(args[3])) /
        (args[1]^args[4] + args[3]^args[4])^2
end

"""
    hillr(X,v,K,n) = v*(K^n) / (X^n + K^n)

A repressive Hill rate function.
"""
hillr(X, v, K, n) = v * (K^n) / (X^n + K^n)
@register_symbolic hillr(X, v, K, n);
function Symbolics.derivative(::typeof(hillr), args::NTuple{4, Any}, ::Val{1})
    return -args[2] * args[4] * (args[3]^args[4]) * (args[1]^(args[4] - 1)) /
        (args[1]^args[4] + args[3]^args[4])^2
end
function Symbolics.derivative(::typeof(hillr), args::NTuple{4, Any}, ::Val{2})
    return (args[3]^args[4]) / (args[1]^args[4] + args[3]^args[4])
end
function Symbolics.derivative(::typeof(hillr), args::NTuple{4, Any}, ::Val{3})
    return args[2] * args[4] * (args[1]^args[4]) * (args[3]^(args[4] - 1)) /
        (args[1]^args[4] + args[3]^args[4])^2
end
function Symbolics.derivative(::typeof(hillr), args::NTuple{4, Any}, ::Val{4})
    return args[2] * (args[1]^args[4]) * (args[3]^args[4]) * (log(args[3]) - log(args[1])) /
        (args[1]^args[4] + args[3]^args[4])^2
end

"""
    hillar(X,Y,v,K,n) = v*(X^n) / (X^n + Y^n + K^n)

An activation/repressing Hill rate function.
"""
hillar(X, Y, v, K, n) = v * (X^n) / (X^n + Y^n + K^n)
@register_symbolic hillar(X, Y, v, K, n);
function Symbolics.derivative(::typeof(hillar), args::NTuple{5, Any}, ::Val{1})
    return args[3] * args[5] * (args[1]^(args[5] - 1)) * (args[2]^args[5] + args[4]^args[5]) /
        (args[1]^args[5] + args[2]^args[5] + args[4]^args[5])^2
end
function Symbolics.derivative(::typeof(hillar), args::NTuple{5, Any}, ::Val{2})
    return -args[3] * args[5] * (args[2]^(args[5] - 1)) * (args[1]^args[5]) /
        (args[1]^args[5] + args[2]^args[5] + args[4]^args[5])^2
end
function Symbolics.derivative(::typeof(hillar), args::NTuple{5, Any}, ::Val{3})
    return (args[1]^args[5]) / (args[1]^args[5] + args[2]^args[5] + args[4]^args[5])
end
function Symbolics.derivative(::typeof(hillar), args::NTuple{5, Any}, ::Val{4})
    return -args[3] * args[5] * (args[3]^(args[5] - 1)) * (args[1]^args[5]) /
        (args[1]^args[5] + args[2]^args[5] + args[4]^args[5])^2
end
function Symbolics.derivative(::typeof(hillar), args::NTuple{5, Any}, ::Val{5})
    return args[3] * (args[1]^args[5]) *
        (
        log(args[1]) * (args[2]^args[5] + args[4]^args[5]) - (args[2]^args[5]) * log(args[2]) -
            (args[4]^args[5]) * log(args[4])
    ) /
        (args[1]^args[5] + args[2]^args[5] + args[4]^args[5])^2
end

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
    return Reaction(
        expand_registered_functions(rx.rate), rx.substrates, rx.products,
        rx.substoich, rx.prodstoich, rx.netstoich, rx.only_use_rate, rx.metadata
    )
end

# If applied to a Equation, returns it with it applied to lhs and rhs.
function expand_registered_functions(eq::Equation)
    return expand_registered_functions(eq.lhs) ~ expand_registered_functions(eq.rhs)
end

# If applied to a continuous event, returns it applied to eqs and affect.
function expand_registered_functions(ce::ModelingToolkit.SymbolicContinuousCallback)
    eqs = expand_registered_functions(ce.eqs)
    affect = expand_registered_functions(ce.affect)
    return ModelingToolkit.SymbolicContinuousCallback(eqs, affect)
end

# If applied to a discrete event, returns it applied to condition and affects.
function expand_registered_functions(de::ModelingToolkit.SymbolicDiscreteCallback)
    condition = expand_registered_functions(de.condition)
    affects = expand_registered_functions(de.affects)
    return ModelingToolkit.SymbolicDiscreteCallback(condition, affects)
end

# If applied to a vector, applies it to every element in the vector.
function expand_registered_functions(vec::Vector)
    return [Catalyst.expand_registered_functions(element) for element in vec]
end

# If applied to a ReactionSystem, applied function to all Reactions and other Equations, and return updated system.
# Currently, `ModelingToolkit.has_X_events` returns `true` even if event vector is empty (hence
# this function cannot be used).
function expand_registered_functions(rs::ReactionSystem)
    @set! rs.eqs = Catalyst.expand_registered_functions(get_eqs(rs))
    @set! rs.rxs = Catalyst.expand_registered_functions(get_rxs(rs))
    if !isempty(ModelingToolkit.get_continuous_events(rs))
        @set! rs.continuous_events = Catalyst.expand_registered_functions(ModelingToolkit.get_continuous_events(rs))
    end
    if !isempty(ModelingToolkit.get_discrete_events(rs))
        @set! rs.discrete_events = Catalyst.expand_registered_functions(ModelingToolkit.get_discrete_events(rs))
    end
    reset_networkproperties!(rs)
    return rs
end
