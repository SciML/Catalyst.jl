### Unit Inference Helpers ###

# SymbolicDimensions-preserving unit inference for Catalyst's unit validation.
# Unlike MTKBase's `get_unit` (which expands to concrete SI Dimensions via
# `uexpand`, introducing floating-point precision loss), this preserves
# the user-specified unit system (e.g., M, μM) by reading raw metadata
# and performing exact SymbolicDimensions arithmetic.

const DQ = DynamicQuantities

# Sentinel for unitless quantities in SymbolicDimensions space.
# We construct a Quantity{Float64, SymbolicDimensions} with all-zero exponents.
const SYM_UNITLESS = DQ.Quantity(1.0, DQ.SymbolicDimensions())
const Conditional = Union{typeof(ifelse)}
const Comparison = Union{typeof.([==, !=, <, <=, >, >=])...}
const UnitValidationValue = Union{DQ.AbstractQuantity, SymbolicT}

"""
    UnitValidationIssue

A single unit-validation diagnostic entry.
"""
struct UnitValidationIssue
    """
    Machine-readable issue category.
    Current values include: `:species_unit_mismatch`,
    `:reaction_species_unit_mismatch`, `:reaction_rate_unit_mismatch`,
    `:additive_term_unit_mismatch`, `:equation_unit_mismatch`,
    `:reaction_side_unit_mismatch`, `:comparison_unit_mismatch`,
    `:conditional_condition_unit_mismatch`,
    `:conditional_branch_unit_mismatch`, `:exponent_unit_mismatch`,
    `:symbolic_stoichiometry`, and `:symbolic_exponent`.
    """
    kind::Symbol
    """Human-readable location/context where the issue was found."""
    context::String
    """The expected/reference unit for the check (or first compared unit)."""
    lhs_unit::Union{Nothing, UnitValidationValue}
    """The observed/compared-against unit."""
    rhs_unit::Union{Nothing, UnitValidationValue}
    """Additional explanatory detail."""
    detail::String
end

"""
    UnitValidationReport

Structured output from unit validation.
"""
struct UnitValidationReport
    """`true` if no unit issues were found, otherwise `false`."""
    valid::Bool
    """All collected diagnostics."""
    issues::Vector{UnitValidationIssue}
end

"""
    UnitValidationError

Exception thrown by strict unit-validation entrypoints (e.g.
`assert_valid_units`). Wraps the full validation report.
"""
struct UnitValidationError <: Exception
    """The full set of collected unit issues."""
    report::UnitValidationReport
    """Optional caller-provided context string (e.g. constructor call site)."""
    info::String
end

function Base.showerror(io::IO, err::UnitValidationError)
    nissues = length(err.report.issues)
    print(io, "UnitValidationError: unit validation failed (", nissues, " issue")
    (nissues == 1) || print(io, "s")
    print(io, ").")
    isempty(err.info) || print(io, "\ninfo: ", err.info)
    for (i, issue) in enumerate(err.report.issues)
        print(io, "\n\n[", i, "] ", issue.kind)
        print(io, "\ncontext: ", issue.context)
        (issue.lhs_unit === nothing) || print(io, "\nlhs unit: ", issue.lhs_unit)
        (issue.rhs_unit === nothing) || print(io, "\nrhs unit: ", issue.rhs_unit)
        isempty(issue.detail) || print(io, "\ndetail: ", issue.detail)
    end
end

# Emit @warn messages for a vector of UnitValidationIssues.
# Used by both validate_units(::Reaction) and validate_units(::ReactionSystem).
function _warn_unit_issues(issues::Vector{UnitValidationIssue})
    for issue in issues
        if issue.kind == :species_unit_mismatch
            @warn(string("Species are expected to have units of ", issue.lhs_unit,
                " however, species ", issue.context, " has units ", issue.rhs_unit, "."))
        elseif issue.kind == :reaction_species_unit_mismatch
            @warn(string("In ", issue.context, " the ", issue.detail, "."))
        elseif issue.kind == :reaction_rate_unit_mismatch
            @warn(string(
                "Reaction rate laws are expected to have units of ", issue.lhs_unit,
                " however, ", issue.context, " has units of ", issue.rhs_unit, "."))
        elseif issue.kind == :reaction_side_unit_mismatch
            @warn(string("in ", issue.context,
                " the substrate units are not consistent with the product units."))
        elseif issue.kind == :additive_term_unit_mismatch
            @warn(string(issue.context, ": additive terms have mismatched units [",
                issue.lhs_unit, "] and [", issue.rhs_unit, "]."))
        elseif issue.kind == :equation_unit_mismatch
            @warn(string("Equation unit mismatch in ", issue.context,
                ": lhs has units ", issue.lhs_unit, ", rhs has units ", issue.rhs_unit, "."))
        elseif issue.kind == :comparison_unit_mismatch
            @warn(string(issue.context, ": comparison operands have mismatched units [",
                issue.lhs_unit, "] and [", issue.rhs_unit, "]."))
        elseif issue.kind == :conditional_condition_unit_mismatch
            @warn(string(issue.context, ": ifelse condition must be unitless, got [",
                issue.rhs_unit, "]."))
        elseif issue.kind == :conditional_branch_unit_mismatch
            @warn(string(issue.context, ": ifelse branches have mismatched units [",
                issue.lhs_unit, "] and [", issue.rhs_unit, "]."))
        elseif issue.kind == :exponent_unit_mismatch
            @warn(string(issue.context, ": exponent must be unitless, got [",
                issue.rhs_unit, "]."))
        elseif issue.kind == :symbolic_stoichiometry
            @warn(string(issue.context, ": ", issue.detail, "."))
        elseif issue.kind == :symbolic_exponent
            @warn(string(issue.context, ": ", issue.detail, "."))
        end
    end
end

_validation_push!(issues, issue::UnitValidationIssue) = isnothing(issues) ? nothing : push!(issues, issue)
_units_match(lhs, rhs) = isequal(lhs, rhs)
_is_unitless(unit) = _units_match(unit, SYM_UNITLESS)

"""
    catalyst_get_unit(x, [noise_units])

Get the unit of a symbolic expression, preserving SymbolicDimensions.

Unlike MTKBase's `get_unit`, this does not expand to concrete SI Dimensions,
avoiding floating-point precision loss when using non-SI units like M or μM.
Returns a `DynamicQuantities.Quantity{Float64, SymbolicDimensions}`.

The optional `noise_units` argument is a `Dict` mapping unwrapped noise variable
symbols to their effective units. This is built by `validate_units(rs::ReactionSystem)`
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
        return raw
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
    #    is done separately in unit validation checks.
    if SymbolicUtils.isadd(x)
        return catalyst_get_unit(SymbolicUtils.arguments(x)[1], noise_units)
    end

    # 5. Power: base^exponent
    if SymbolicUtils.ispow(x)
        pargs = SymbolicUtils.arguments(x)
        base_u = catalyst_get_unit(pargs[1], noise_units)
        if _is_unitless(base_u)
            return SYM_UNITLESS
        end
        exp_u = catalyst_get_unit(pargs[2], noise_units)
        _is_unitless(exp_u) || return SYM_UNITLESS
        exp_val = pargs[2]
        if SymbolicUtils.isconst(exp_val)
            return base_u ^ value(exp_val)
        else
            # Symbolic exponent on unitful base — unit is indeterminate.
            # Validation code pre-checks via _has_symbolic_unitful_pow before
            # reaching here, so this is a safe fallback.
            return SYM_UNITLESS
        end
    end

    # 6. Function calls
    if iscall(x)
        op = SymbolicUtils.operation(x)
        args = SymbolicUtils.arguments(x)

        # Dependent variable calls like A(t) — look up metadata on the term
        if SymbolicUtils.issym(op)
            raw2 = getmetadata(x, MT.VariableUnit, nothing)
            return raw2 isa DQ.AbstractQuantity ? raw2 : SYM_UNITLESS
        end

        # Indexed access: getindex(array_var, i...) — unit comes from the array.
        if op === getindex
            return catalyst_get_unit(args[1], noise_units)
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
# (term consistency checked separately in unit validation checks)
function _cgu_op(::Union{typeof(+), typeof(-)}, args, noise_units)
    return catalyst_get_unit(args[1], noise_units)
end

# Comparison operators are dimensionless.
_cgu_op(::Comparison, args, noise_units) = SYM_UNITLESS

# ifelse(condition, x, y): condition must be unitless and branch units must match.
function _cgu_op(::Conditional, args, noise_units)
    (length(args) == 3) || return SYM_UNITLESS
    cond_unit = catalyst_get_unit(args[1], noise_units)
    true_unit = catalyst_get_unit(args[2], noise_units)
    false_unit = catalyst_get_unit(args[3], noise_units)
    return (_is_unitless(cond_unit) && _units_match(true_unit, false_unit)) ? true_unit : SYM_UNITLESS
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
        if result isa DQ.AbstractQuantity
            return oneunit(result)
        end
        if (result isa Bool) || (result isa Number)
            return SYM_UNITLESS
        end
        @warn "Cannot infer unit for $op with units $unit_args: unsupported return type $(typeof(result))."
        return SYM_UNITLESS
    catch e
        @warn "Cannot infer unit for $op with units $unit_args: $e"
        return SYM_UNITLESS
    end
end

"""
    _has_symbolic_unitful_pow(x) -> Bool

Check whether a symbolic expression contains any power node `base^exp` where the
base has units and the exponent is symbolic (not a constant). Such expressions have
indeterminate units at analysis time.
"""
function _has_symbolic_unitful_pow(x)
    x = Symbolics.unwrap(x)
    x isa SymbolicT || return false
    if SymbolicUtils.ispow(x)
        args = SymbolicUtils.arguments(x)
        # Only match unitful base with unitless symbolic exponent (e.g. A^n).
        # Unitful exponents (e.g. V^t) are separately caught as :exponent_unit_mismatch.
        if !_is_unitless(catalyst_get_unit(args[1])) &&
           !SymbolicUtils.isconst(args[2]) &&
           _is_unitless(catalyst_get_unit(args[2]))
            return true
        end
    end
    if iscall(x) || SymbolicUtils.isadd(x)
        return any(_has_symbolic_unitful_pow, SymbolicUtils.arguments(x))
    end
    return false
end

"""
    _validate_unit_expr(expr, label, [noise_units]) -> Bool

Walk a symbolic expression tree and check that all additive terms have consistent
units. Returns `true` if valid, `false` if any addition has terms with mismatched
units. Issues warnings for each mismatch found.
"""
function _validate_unit_expr(expr, label, noise_units=nothing; issues=nothing, warn::Bool = true)
    x = expr isa Symbolics.Num ? Symbolics.unwrap(expr) : expr
    x isa SymbolicT || return true
    return _vue_walk(x, label, noise_units; issues, warn)
end

function _vue_walk(x, label, noise_units; issues=nothing, warn::Bool = true)
    valid = true
    if SymbolicUtils.isadd(x)
        args = SymbolicUtils.arguments(x)
        first_unit = catalyst_get_unit(args[1], noise_units)
        for i in 2:length(args)
            other_unit = catalyst_get_unit(args[i], noise_units)
            if !_units_match(first_unit, other_unit)
                valid = false
                issue = UnitValidationIssue(:additive_term_unit_mismatch, string(label),
                    first_unit, other_unit, "Additive terms have mismatched units.")
                _validation_push!(issues, issue)
                warn && @warn("$label: additive terms have mismatched units [$first_unit] and [$other_unit] in $x.")
            end
        end
        # Also recurse into each term
        for a in args
            a_unwrapped = a isa Symbolics.Num ? Symbolics.unwrap(a) : a
            if a_unwrapped isa SymbolicT
                valid &= _vue_walk(a_unwrapped, label, noise_units; issues, warn)
            end
        end
    elseif SymbolicUtils.ispow(x)
        args = SymbolicUtils.arguments(x)
        (length(args) == 2) || return valid
        exp_unit = catalyst_get_unit(args[2], noise_units)
        if !_is_unitless(exp_unit)
            valid = false
            issue = UnitValidationIssue(:exponent_unit_mismatch, string(label), SYM_UNITLESS,
                exp_unit, "Exponent must be unitless.")
            _validation_push!(issues, issue)
            warn && @warn("$label: exponent has non-unitless unit [$exp_unit] in $x.")
        end
        for a in args
            a_unwrapped = a isa Symbolics.Num ? Symbolics.unwrap(a) : a
            if a_unwrapped isa SymbolicT
                valid &= _vue_walk(a_unwrapped, label, noise_units; issues, warn)
            end
        end
    elseif iscall(x)
        op = SymbolicUtils.operation(x)
        args = SymbolicUtils.arguments(x)

        if op isa Comparison
            (length(args) >= 2) || return valid
            lhs_unit = catalyst_get_unit(args[1], noise_units)
            rhs_unit = catalyst_get_unit(args[2], noise_units)
            if !_units_match(lhs_unit, rhs_unit)
                valid = false
                issue = UnitValidationIssue(:comparison_unit_mismatch, string(label), lhs_unit,
                    rhs_unit, "Comparison operands must have matching units.")
                _validation_push!(issues, issue)
                warn && @warn("$label: comparison operands have mismatched units [$lhs_unit] and [$rhs_unit] in $x.")
            end
        elseif op isa Conditional
            (length(args) == 3) || return valid
            cond_unit = catalyst_get_unit(args[1], noise_units)
            if !_is_unitless(cond_unit)
                valid = false
                issue = UnitValidationIssue(:conditional_condition_unit_mismatch, string(label),
                    SYM_UNITLESS, cond_unit, "ifelse condition must be unitless.")
                _validation_push!(issues, issue)
                warn && @warn("$label: ifelse condition is not unitless [$cond_unit] in $x.")
            end
            true_unit = catalyst_get_unit(args[2], noise_units)
            false_unit = catalyst_get_unit(args[3], noise_units)
            if !_units_match(true_unit, false_unit)
                valid = false
                issue = UnitValidationIssue(:conditional_branch_unit_mismatch, string(label),
                    true_unit, false_unit, "ifelse branches must have matching units.")
                _validation_push!(issues, issue)
                warn && @warn("$label: ifelse branches have mismatched units [$true_unit] and [$false_unit] in $x.")
            end
        end

        for a in args
            a_unwrapped = a isa Symbolics.Num ? Symbolics.unwrap(a) : a
            if a_unwrapped isa SymbolicT
                valid &= _vue_walk(a_unwrapped, label, noise_units; issues, warn)
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
function _validate_equation(eq::Equation; noise_units=nothing, issues=nothing, warn::Bool = true)
    valid = true
    eq_str = string(eq)
    if _has_symbolic_unitful_pow(eq.lhs) || _has_symbolic_unitful_pow(eq.rhs)
        valid = false
        issue = UnitValidationIssue(:symbolic_exponent, eq_str, nothing, nothing,
            "Symbolic exponent on unitful base is not supported for unit validation")
        _validation_push!(issues, issue)
        warn && @warn(string(eq_str, ": symbolic exponent on unitful base, ",
            "cannot validate units."))
        return valid
    end
    valid &= _validate_unit_expr(eq.lhs, string(eq_str, ": lhs"), noise_units; issues, warn)
    valid &= _validate_unit_expr(eq.rhs, string(eq_str, ": rhs"), noise_units; issues, warn)
    lhs_unit = catalyst_get_unit(eq.lhs, noise_units)
    rhs_unit = catalyst_get_unit(eq.rhs, noise_units)
    if !_units_match(lhs_unit, rhs_unit)
        valid = false
        issue = UnitValidationIssue(:equation_unit_mismatch, string(eq), lhs_unit,
            rhs_unit, "LHS and RHS units do not match.")
        _validation_push!(issues, issue)
        warn && @warn(string("Equation unit mismatch in ", eq,
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
