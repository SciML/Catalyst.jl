### Custom CRN Function Implementations ###
"""
    mm(X,v,K) = v*X / (X + K)

A Michaelis-Menten rate function.
"""
mm(X, v, K) = v * X / (X + K)
@register_symbolic mm(X, v, K);
function Symbolics.derivative(::typeof(mm), args::NTuple{3, Any}, ::Val{1})
    (args[2] * args[3]) / (args[1] + args[3])^2
end
function Symbolics.derivative(::typeof(mm), args::NTuple{3, Any}, ::Val{2})
    args[1] / (args[1] + args[3])
end
function Symbolics.derivative(::typeof(mm), args::NTuple{3, Any}, ::Val{3})
    -args[2] * args[1] / (args[1] + args[3])^2
end

# Registers the repressing Michaelis-Menten function.
"""
    mmr(X,v,K) = v*K / (X + K)

A repressive Michaelis-Menten rate function.
"""
mmr(X, v, K) = v * K / (X + K)
@register_symbolic mmr(X, v, K);
function Symbolics.derivative(::typeof(mmr), args::NTuple{3, Any}, ::Val{1})
    -(args[2] * args[3]) / (args[1] + args[3])^2
end
function Symbolics.derivative(::typeof(mmr), args::NTuple{3, Any}, ::Val{2})
    args[3] / (args[1] + args[3])
end
function Symbolics.derivative(::typeof(mmr), args::NTuple{3, Any}, ::Val{3})
    args[2] * args[1] / (args[1] + args[3])^2
end

"""
    hill(X,v,K,n) = v*(X^n) / (X^n + K^n)

A Hill rate function.
"""
hill(X, v, K, n) = v * (X^n) / (X^n + K^n)
@register_symbolic hill(X, v, K, n);
function Symbolics.derivative(::typeof(hill), args::NTuple{4, Any}, ::Val{1})
    args[2] * args[4] * (args[3]^args[4]) * (args[1]^(args[4] - 1)) /
    (args[1]^args[4] + args[3]^args[4])^2
end
function Symbolics.derivative(::typeof(hill), args::NTuple{4, Any}, ::Val{2})
    (args[1]^args[4]) / (args[1]^args[4] + args[3]^args[4])
end
function Symbolics.derivative(::typeof(hill), args::NTuple{4, Any}, ::Val{3})
    -args[2] * args[4] * (args[1]^args[4]) * (args[3]^(args[4] - 1)) /
    (args[1]^args[4] + args[3]^args[4])^2
end
function Symbolics.derivative(::typeof(hill), args::NTuple{4, Any}, ::Val{4})
    args[2] * (args[1]^args[4]) * (args[3]^args[4]) * (log(args[1]) - log(args[3])) /
    (args[1]^args[4] + args[3]^args[4])^2
end

"""
    hillr(X,v,K,n) = v*(K^n) / (X^n + K^n)

A repressive Hill rate function.
"""
hillr(X, v, K, n) = v * (K^n) / (X^n + K^n)
@register_symbolic hillr(X, v, K, n);
function Symbolics.derivative(::typeof(hillr), args::NTuple{4, Any}, ::Val{1})
    -args[2] * args[4] * (args[3]^args[4]) * (args[1]^(args[4] - 1)) /
    (args[1]^args[4] + args[3]^args[4])^2
end
function Symbolics.derivative(::typeof(hillr), args::NTuple{4, Any}, ::Val{2})
    (args[3]^args[4]) / (args[1]^args[4] + args[3]^args[4])
end
function Symbolics.derivative(::typeof(hillr), args::NTuple{4, Any}, ::Val{3})
    args[2] * args[4] * (args[1]^args[4]) * (args[3]^(args[4] - 1)) /
    (args[1]^args[4] + args[3]^args[4])^2
end
function Symbolics.derivative(::typeof(hillr), args::NTuple{4, Any}, ::Val{4})
    args[2] * (args[1]^args[4]) * (args[3]^args[4]) * (log(args[3]) - log(args[1])) /
    (args[1]^args[4] + args[3]^args[4])^2
end

"""
    hillar(X,Y,v,K,n) = v*(X^n) / (X^n + Y^n + K^n)

An activation/repressing Hill rate function.
"""
hillar(X, Y, v, K, n) = v * (X^n) / (X^n + Y^n + K^n)
@register_symbolic hillar(X, Y, v, K, n);
function Symbolics.derivative(::typeof(hillar), args::NTuple{5, Any}, ::Val{1})
    args[3] * args[5] * (args[1]^(args[5] - 1)) * (args[2]^args[5] + args[4]^args[5]) /
    (args[1]^args[5] + args[2]^args[5] + args[4]^args[5])^2
end
function Symbolics.derivative(::typeof(hillar), args::NTuple{5, Any}, ::Val{2})
    -args[3] * args[5] * (args[2]^(args[5] - 1)) * (args[1]^args[5]) /
    (args[1]^args[5] + args[2]^args[5] + args[4]^args[5])^2
end
function Symbolics.derivative(::typeof(hillar), args::NTuple{5, Any}, ::Val{3})
    (args[1]^args[5]) / (args[1]^args[5] + args[2]^args[5] + args[4]^args[5])
end
function Symbolics.derivative(::typeof(hillar), args::NTuple{5, Any}, ::Val{4})
    -args[3] * args[5] * (args[3]^(args[5] - 1)) * (args[1]^args[5]) /
    (args[1]^args[5] + args[2]^args[5] + args[4]^args[5])^2
end
function Symbolics.derivative(::typeof(hillar), args::NTuple{5, Any}, ::Val{5})
    args[3] * (args[1]^args[5]) *
    (log(args[1]) * (args[2]^args[5] + args[4]^args[5]) - (args[2]^args[5]) * log(args[2]) -
     (args[4]^args[5]) * log(args[4])) /
    (args[1]^args[5] + args[2]^args[5] + args[4]^args[5])^2
end

### Custom CRN FUnction-related Functions ###

"""
expand_registered_functions(expr)

Takes an expression, and expands registered function expressions. E.g. `mm(X,v,K)` is replaced with v*X/(X+K). Currently supported functions: `mm`, `mmr`, `hill`, `hillr`, and `hill`.
"""
function expand_registered_functions(expr)
    istree(expr) || return expr
    args = arguments(expr)
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
    for i in 1:length(args)
        args[i] = expand_registered_functions(args[i])
    end
    return expr
end
# If applied to a Reaction, return a reaction with its rate modified.
function expand_registered_functions(rx::Reaction)
    Reaction(expand_registered_functions(rx.rate), rx.substrates, rx.products, 
        rx.substoich, rx.prodstoich, rx.netstoich, rx.only_use_rate, rx.metadata)
end
# If applied to a Equation, returns it with it applied to lhs and rhs
function expand_registered_functions(eq::Equation)
    return expand_registered_functions(eq.lhs) ~ expand_registered_functions(eq.rhs)
end
# If applied to a ReactionSystem, applied function to all Reactions and other Equations, and return updated system.
function expand_registered_functions(rs::ReactionSystem)
    @set! rs.eqs = [Catalyst.expand_registered_functions(eq) for eq in get_eqs(rs)]
    @set! rs.rxs = [Catalyst.expand_registered_functions(rx) for rx in get_rxs(rs)]
    return rs
end
