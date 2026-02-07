### Unit Inference Helpers ###

# SymbolicDimensions-preserving unit inference for Catalyst's validate().
# Unlike MTKBase's `get_unit` (which expands to concrete SI Dimensions via
# `uexpand`, introducing floating-point precision loss), this preserves
# the user-specified unit system (e.g., M, μM) by reading raw metadata
# and performing exact SymbolicDimensions arithmetic.

const DQ = DynamicQuantities

# Sentinel for unitless quantities in SymbolicDimensions space.
# We construct a Quantity{Float64, SymbolicDimensions} with all-zero exponents.
const SYM_UNITLESS = DQ.Quantity(1.0, DQ.SymbolicDimensions())

"""
    catalyst_get_unit(x, [noise_units])

Get the unit of a symbolic expression, preserving SymbolicDimensions.

Unlike MTKBase's `get_unit`, this does not expand to concrete SI Dimensions,
avoiding floating-point precision loss when using non-SI units like M or μM.
Returns a `DynamicQuantities.Quantity{Float64, SymbolicDimensions}`.

The optional `noise_units` argument is a `Dict` mapping unwrapped noise variable
symbols to their effective units. This is built by `validate(rs::ReactionSystem)`
from the system's brownian and poissonian variables:
- Brownians get units `time^(-1/2)` (from the Wiener process derivative `dW/dt`).
- Poissonians get units from their associated rate parameter (via `getpoissonianrate`).
"""
function catalyst_get_unit(x, noise_units=nothing)
    # Unwrap Symbolics wrappers (Num, CallAndWrap, etc.) to the underlying symbolic
    x = Symbolics.unwrap(x)
    if x isa DQ.AbstractQuantity
        return x
    elseif x isa SymbolicT
        return _cgu_symbolic(x, noise_units)
    else
        return SYM_UNITLESS
    end
end

# Main symbolic dispatch
function _cgu_symbolic(x, noise_units)
    # 1. Constant (numeric wrapped in symbolic)
    if SymbolicUtils.isconst(x)
        return catalyst_get_unit(value(x), noise_units)
    end

    # 2. Literal unit metadata on this node
    raw = getmetadata(x, MT.VariableUnit, nothing)
    if raw isa DQ.AbstractQuantity
        return _ensure_sym(raw)
    end

    # 3. Bare symbol without unit metadata
    if SymbolicUtils.issym(x)
        # Check noise variable unit overrides (brownians/poissonians)
        if noise_units !== nothing && haskey(noise_units, x)
            return noise_units[x]
        end
        return SYM_UNITLESS
    end

    # 4. Addition: return first term's unit (validation of term consistency
    #    is done separately in validate())
    if SymbolicUtils.isadd(x)
        return catalyst_get_unit(SymbolicUtils.arguments(x)[1], noise_units)
    end

    # 5. Power: base^exponent
    if SymbolicUtils.ispow(x)
        pargs = SymbolicUtils.arguments(x)
        base_u = catalyst_get_unit(pargs[1], noise_units)
        if base_u == SYM_UNITLESS
            return SYM_UNITLESS
        end
        exp_val = pargs[2]
        if SymbolicUtils.isconst(exp_val)
            return base_u ^ value(exp_val)
        else
            return (1 * base_u) ^ exp_val
        end
    end

    # 6. Function calls
    if iscall(x)
        op = SymbolicUtils.operation(x)
        args = SymbolicUtils.arguments(x)

        # Dependent variable calls like A(t) — look up metadata on the term
        if SymbolicUtils.issym(op) ||
           (iscall(op) && iscall(SymbolicUtils.operation(op)))
            raw2 = getmetadata(x, MT.VariableUnit, nothing)
            return raw2 isa DQ.AbstractQuantity ? _ensure_sym(raw2) : SYM_UNITLESS
        end

        # Indexed access: getindex(array_var, i...) — unit comes from the array.
        # Handles both plain getindex (op === getindex) and nested call forms.
        if op === getindex
            return catalyst_get_unit(args[1], noise_units)
        end
        if iscall(op) && SymbolicUtils.operation(op) === Base.getindex
            gp = SymbolicUtils.arguments(op)[1]
            raw3 = getmetadata(gp, MT.VariableUnit, nothing)
            return raw3 isa DQ.AbstractQuantity ? _ensure_sym(raw3) : SYM_UNITLESS
        end

        # Actual operations (*, /, +, -, Differential, registered functions, etc.)
        return _cgu_op(op, args, noise_units)
    end

    return SYM_UNITLESS
end

# Multiplication: multiply all argument units
function _cgu_op(::typeof(*), args, noise_units)
    result = SYM_UNITLESS
    for a in args
        result *= catalyst_get_unit(a, noise_units)
    end
    return result
end

# Division
function _cgu_op(::typeof(/), args, noise_units)
    return catalyst_get_unit(args[1], noise_units) / catalyst_get_unit(args[2], noise_units)
end

# Addition/subtraction: return first argument's unit
# (term consistency checked separately in validate)
function _cgu_op(::Union{typeof(+), typeof(-)}, args, noise_units)
    return catalyst_get_unit(args[1], noise_units)
end

# Differential: unit(arg) / unit(iv)^order
function _cgu_op(op::Differential, args, noise_units)
    return catalyst_get_unit(args[1], noise_units) / catalyst_get_unit(op.x, noise_units)^op.order
end

# Generic fallback: call the function with unit arguments and use oneunit to
# strip the numeric coefficient. This handles registered functions like mm, hill, etc.
function _cgu_op(op, args, noise_units)
    unit_args = [catalyst_get_unit(a, noise_units) for a in args]
    try
        result = op(unit_args...)
        return oneunit(result)
    catch e
        @warn "Cannot infer unit for $op with units $unit_args: $e"
        return SYM_UNITLESS
    end
end

# Ensure we have SymbolicDimensions form. If user used u"..." instead of us"...",
# we keep it as-is (concrete Dimensions) — this won't break arithmetic but
# won't benefit from exact comparison.
function _ensure_sym(u::DQ.AbstractQuantity)
    return u
end

"""
    _validate_unit_expr(expr, label, [noise_units]) -> Bool

Walk a symbolic expression tree and check that all additive terms have consistent
units. Returns `true` if valid, `false` if any addition has terms with mismatched
units. Issues warnings for each mismatch found.
"""
function _validate_unit_expr(expr, label, noise_units=nothing)
    x = expr isa Symbolics.Num ? Symbolics.unwrap(expr) : expr
    x isa SymbolicT || return true
    return _vue_walk(x, label, noise_units)
end

function _vue_walk(x, label, noise_units)
    valid = true
    if SymbolicUtils.isadd(x)
        args = SymbolicUtils.arguments(x)
        first_unit = catalyst_get_unit(args[1], noise_units)
        for i in 2:length(args)
            other_unit = catalyst_get_unit(args[i], noise_units)
            if first_unit != other_unit
                valid = false
                @warn "$label: additive terms have mismatched units [$first_unit] and [$other_unit] in $x."
            end
        end
        # Also recurse into each term
        for a in args
            a_unwrapped = a isa Symbolics.Num ? Symbolics.unwrap(a) : a
            if a_unwrapped isa SymbolicT
                valid &= _vue_walk(a_unwrapped, label, noise_units)
            end
        end
    elseif iscall(x)
        for a in SymbolicUtils.arguments(x)
            a_unwrapped = a isa Symbolics.Num ? Symbolics.unwrap(a) : a
            if a_unwrapped isa SymbolicT
                valid &= _vue_walk(a_unwrapped, label, noise_units)
            end
        end
    end
    return valid
end

"""
    _validate_equation(eq::Equation; noise_units=nothing) -> Bool

Validate units in an `Equation`. Checks:
1. Internal consistency of additive terms on each side.
2. That LHS and RHS have matching overall units.

The optional `noise_units` keyword is a `Dict` mapping noise variable symbols
to their effective units (built by `_build_noise_units` during system validation).
"""
function _validate_equation(eq::Equation; noise_units=nothing)
    valid = true
    eq_str = string(eq)
    valid &= _validate_unit_expr(eq.lhs, string(eq_str, ": lhs"), noise_units)
    valid &= _validate_unit_expr(eq.rhs, string(eq_str, ": rhs"), noise_units)
    lhs_unit = catalyst_get_unit(eq.lhs, noise_units)
    rhs_unit = catalyst_get_unit(eq.rhs, noise_units)
    if lhs_unit != rhs_unit
        valid = false
        @warn(string("Equation unit mismatch in ", eq,
            ": lhs has units ", lhs_unit, ", rhs has units ", rhs_unit, "."))
    end
    return valid
end

"""
    _build_noise_units(rs) -> Dict

Build a dictionary mapping noise variable symbols to their effective units
for use in equation unit validation.

- Brownians: `dW/dt` has units `time^(-1/2)` (Wiener process derivative).
- Poissonians: `dN/dt` has the same units as the associated rate parameter
  (accessed via `MT.getpoissonianrate`).
"""
function _build_noise_units(rs)
    timeunit = catalyst_get_unit(get_iv(rs))
    noise_units = Dict{Any, Any}()
    for b in MT.get_brownians(rs)
        noise_units[Symbolics.unwrap(b)] = timeunit^(-1//2)
    end
    for p in MT.get_poissonians(rs)
        rate = MT.getpoissonianrate(p)
        noise_units[Symbolics.unwrap(p)] = rate !== nothing ? catalyst_get_unit(rate) : SYM_UNITLESS
    end
    return noise_units
end
