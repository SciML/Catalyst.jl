### Structures & Constants ###

Base.@kwdef mutable struct CatalystLatexParams
    double_linebreak::Bool = false
    starred::Bool = true
    harpoon_arrows::Bool = true
    mathjax::Bool = false
end

const LATEX_DEFS = CatalystLatexParams()

# Wrappers through which `latexify` displays a `ReactionSystem` as its chemical Langevin
# equation (`form = :sde`) or as the random time change representation of its jump process
# (`form = :jump`).
struct LatexCLE{R}
    rs::R
end

struct LatexTimeChange{R}
    rs::R
end

### Latexify Receipt ###

@latexrecipe function f(rs::ReactionSystem; form = :reactions, expand_functions = true)
    expand_functions && (rs = expand_registered_functions(rs))
    if form == :reactions    # Returns chemical reaction network code.
        mult_symbol --> ""
        env --> :chem
        return rs
    elseif form == :ode      # Returns ODE system code (rendered via chemical_arrows).
        mult_symbol --> ""
        env --> :chem
        return system_to_reactionsystem(ode_model(rs); disable_forbidden_symbol_check = true)
    elseif form == :sde      # Returns the chemical Langevin equation.
        mult_symbol --> ""
        return LatexCLE(rs)
    elseif form == :jump     # Returns the random time change representation.
        mult_symbol --> ""
        return LatexTimeChange(rs)
    end
    error("Unrecognised form argument given: $form. This should be either :reactions (default), :ode, :sde, or :jump.")
end

function Latexify.infer_output(env, rs::ReactionSystem, args...)
    if env in (:arrows, :chem, :chemical, :arrow)
        return chemical_arrows
    else
        error("The environment $env is not defined.")
    end
end

Latexify.infer_output(env, ::LatexCLE, args...) = cle_equations
Latexify.infer_output(env, ::LatexTimeChange, args...) = time_change_equations

function processsym(s)
    args = sorted_arguments(s)
    name = MT.getname(s)
    if length(args) <= 1
        var = value(Symbolics.variable(name))
    else
        idxs = args[2:end]
        var = value(Symbolics.variable(MT.getname(s), idxs...))
    end
    # Rendered as Latexify renders the name, with multi-letter names then made upright by
    # `_upright_names`, as for unknowns that keep their time argument.
    var
end

function chemical_arrows(rn::ReactionSystem;
        double_linebreak = LATEX_DEFS.double_linebreak,
        starred = LATEX_DEFS.starred, mathrm = true,
        mathjax = LATEX_DEFS.mathjax, math_delimiters = false,
        show_time_arg = false, kwargs...)
    any_nonrx_subsys(rn) &&
        (@warn "Latexify currently ignores non-ReactionSystem subsystems. Please call `flatsys = flatten(sys)` to obtain a flattened version of your system before trying to Latexify it.")

    rxs = reactions(rn)
    nonrxs = filter(eq -> eq isa Equation, equations(rn))
    if isempty(rxs) && isempty(nonrxs)
        return Latexify.LaTeXString("ReactionSystem $(MT.nameof(rn)) has no reactions or equations.")
    end

    str = starred ? "\\begin{align*}\n" : "\\begin{align}\n"
    eol = double_linebreak ? "\\\\\\\\\n" : "\\\\\n"
    backwards_reaction = false

    rev_arrow = LATEX_DEFS.harpoon_arrows ? "\\xrightleftharpoons" : "\\xleftrightarrow"

    # test if in IJulia since their mathjax is outdated...
    # VSCODE uses Katex and doesn't have this issue.
    if mathjax || (isdefined(Main, :IJulia) && Main.IJulia.inited &&
        !any(s -> occursin("VSCODE", s), collect(keys(ENV))))
        str *= "\\require{mhchem} \n"
    end

    syms = unknowns(rn)
    if isempty(syms)
        subber = identity
    else
        subber = SymbolicUtils.Substituter{true}([s => processsym(s) for s in syms], SymbolicUtils.default_substitute_filter)
    end
    isvar = _unknown_predicate(rn)

    lastidx = length(rxs)
    for (i, r) in enumerate(rxs)
        if backwards_reaction
            backwards_reaction = false
            continue
        end

        ### Expand functions to maths expressions
        rate = (r.rate isa SymbolicT && !show_time_arg) ? subber(r.rate) : r.rate

        ### Generate formatted string of substrates
        substrates = [make_stoich_str(substrate[1], substrate[2], subber; mathrm,
                          kwargs...)
                      for substrate in zip(r.substrates, r.substoich)]
        isempty(substrates) && (substrates = ["\\varnothing"])

        str *= join(substrates, " + ")

        ### Generate reaction arrows
        if i + 1 <= length(rxs) && issetequal(r.products, rxs[i + 1].substrates) &&
           issetequal(r.substrates, rxs[i + 1].products)
            ### Bi-directional arrows
            rate_backwards = (rxs[i + 1].rate isa SymbolicT && !show_time_arg) ? subber(rxs[i + 1].rate) :
                             rxs[i + 1].rate
            str *= " &" * rev_arrow
            str *= "[" * _latex_ordered(rate_backwards, isvar; kwargs...) * "]"
            str *= "{" * _latex_ordered(rate, isvar; kwargs...) * "} "
            backwards_reaction = true
        else
            ### Uni-directional arrows
            str *= " &\\xrightarrow{" * _latex_ordered(rate, isvar; kwargs...) * "} "
        end

        ### Generate formatted string of products
        products = [make_stoich_str(product[1], product[2], subber; mathrm = true,
                        kwargs...)
                    for product in zip(r.products, r.prodstoich)]
        isempty(products) && (products = ["\\varnothing"])
        str *= join(products, " + ")
        if ((i == lastidx) ||
            (((i + 1) == lastidx) && (backwards_reaction == true))) &&
           isempty(nonrxs)
            str *= "  \n "
        else
            str *= " $eol"
        end
    end

    if !isempty(nonrxs)
        # Without time arguments, the unknowns on the right-hand side are replaced by
        # variables without them, as in the rates, so that e.g. `sin(t)` keeps its argument;
        # the left-hand side (e.g. dX(t)/dt) only drops the time argument of its unknown.
        iv_latex = latexraw(get_iv(rn))
        lhs_str(lhs) = show_time_arg ? latexraw(lhs; kwargs...) :
                       replace(latexraw(lhs; kwargs...), "\\left( $iv_latex \\right)" => "")
        rhs_of(eq) = (show_time_arg || !(eq.rhs isa SymbolicT)) ? eq.rhs : subber(eq.rhs)
        eqstrs = [lhs_str(eq.lhs) * " = " * _latex_ordered(rhs_of(eq), isvar; kwargs...)
                  for eq in nonrxs]
        eqstr_list = replace.(eqstrs, "=" => "&=")
        newstr = join(eqstr_list, " $eol")
        str *= newstr
        str *= "  \n "
    end

    str *= starred ? "\\end{align*}\n" : "\\end{align}\n"

    # Multi-letter names upright, single letters italic; no gap before a call's arguments.
    str = _tight_calls(_upright_names(str))

    # Optionally wrap in $$ delimiters for renderers (Documenter, MathJax) that require
    # explicit math delimiters. Not valid standard LaTeX, but needed for web rendering.
    if math_delimiters
        str = "\$\$\n" * str * "\$\$"
    end

    # `latexify` copies the result to the clipboard itself when `copy_to_clipboard(true)`.
    return Latexify.LaTeXString(str)
end

### SDE and Jump Forms ###

# The chemical Langevin equation of a `ReactionSystem`: for each species X,
#   dX = (∑ⱼ νⱼ aⱼ) dt + ∑ⱼ νⱼ √aⱼ dWⱼ,
# with aⱼ the ODE rate law of reaction j (√aⱼ times its noise scaling, if any), νⱼ its net
# stoichiometry of X, and Wⱼ independent Brownian motions, one per reaction, written
# `\mathcal{W}_j` so that they can't be confused with brownians of the model. Non-reaction
# differential equations are shown in the same form, with the terms of any brownians they
# contain; other equations as they are.
function cle_equations(w::LatexCLE; double_linebreak = LATEX_DEFS.double_linebreak,
        starred = LATEX_DEFS.starred, math_delimiters = false, show_time_arg = false,
        mathrm = true, mathjax = LATEX_DEFS.mathjax, kwargs...)
    rs = Catalyst.flatten(w.rs)
    # `MT.poissonians` is not public in ModelingToolkitBase, unlike `MT.brownians` and
    # `MT.jumps`, an oversight (SciML/ModelingToolkit.jl#5389). `sde_model` uses it too.
    (isempty(MT.jumps(rs)) && isempty(MT.poissonians(rs))) ||
        error("The SDE form can not display a ReactionSystem with coupled jumps or poissonians.")
    rxs = reactions(rs)
    nonrxs = nonreactions(rs)
    (isempty(rxs) && isempty(nonrxs)) && return _no_equations_latex(w.rs)

    subber = show_time_arg ? identity : _latex_subber(rs)
    isvar = _unknown_predicate(rs)
    tex(x) = _latex_ordered(subber(x), isvar; kwargs...)
    ivtex = latexraw(get_iv(rs); kwargs...)
    targ = show_time_arg ? "\\left( " * ivtex * " \\right)" : ""

    # The drift and noise terms that the reactions give each species.
    _, ispcs = get_indep_sts(rs)
    crl = get_combinatoric_ratelaws(rs)
    drift = assemble_oderhs(rs, ispcs; combinatoric_ratelaws = crl,
        expand_catalyst_funs = false)
    noise = [String[] for _ in ispcs]
    species_to_idx = Dict(x => i for (i, x) in enumerate(ispcs))
    for (j, rx) in enumerate(rxs)
        amp = sqrt(oderatelaw(rx; combinatoric_ratelaw = crl, expand_catalyst_funs = false))
        hasnoisescaling(rx) && (amp *= getnoisescaling(rx))
        dW = tex(amp) * " \\, \\mathrm{d}\\mathcal{W}_{$j}" * targ
        for (spec, stoich) in rx.netstoich
            drop_dynamics(spec) && continue
            push!(noise[species_to_idx[spec]], _coef_summand(stoich, dW, tex; kwargs...))
        end
    end
    lines = [_sde_line(tex(x), drift[i], noise[i], tex, ivtex)
             for (i, x) in enumerate(ispcs)]

    # Non-reaction equations, with the terms of any brownians in differential equations.
    browns = MT.brownians(rs)
    for eq in nonrxs
        if Symbolics.is_derivative(eq.lhs) &&
           !Symbolics.is_derivative(only(Symbolics.arguments(eq.lhs)))
            rhs = eq.rhs
            eqnoise = String[]
            for B in browns
                coef, rest, islinear = Symbolics.linear_expansion(rhs, B)
                (islinear && !_iszero(coef)) || continue
                dB = "\\mathrm{d}" * latexraw(B; kwargs...) * targ
                summand = _coef_summand(unwrap(coef), dB, tex; sep = " \\, ", kwargs...)
                push!(eqnoise, summand)
                rhs = unwrap(rest)
            end
            x = only(Symbolics.arguments(eq.lhs))
            push!(lines, _sde_line(tex(x), rhs, eqnoise, tex, ivtex))
        else
            push!(lines, tex(eq.lhs) * " &= " * tex(eq.rhs))
        end
    end
    return _latex_align(lines; starred, double_linebreak, math_delimiters)
end

# The random time change representation of the jump process of a `ReactionSystem`: for each
# species X,
#   X(t) = X(0) + ∑ⱼ νⱼ Rⱼ(∫₀ᵗ aⱼ(s′) ds′),
# with aⱼ the jump rate law (propensity) of reaction j, νⱼ its net stoichiometry of X, and Rⱼ
# independent unit-rate Poisson processes, one per reaction, written `\mathcal{R}_j` so that
# they can't be confused with names in the model. The time arguments are always shown,
# since the integrands need them.
function time_change_equations(w::LatexTimeChange;
        double_linebreak = LATEX_DEFS.double_linebreak, starred = LATEX_DEFS.starred,
        math_delimiters = false, show_time_arg = false, mathrm = true,
        mathjax = LATEX_DEFS.mathjax, kwargs...)
    rs = Catalyst.flatten(w.rs)
    (has_nonreactions(rs) || !isempty(MT.brownians(rs)) || !isempty(MT.jumps(rs))) &&
        error("The jump form can only display reactions, but this ReactionSystem has non-reaction equations, brownians, or coupled jumps.")
    rxs = reactions(rs)
    isempty(rxs) && return _no_equations_latex(w.rs)

    # The rate laws are shown at the integration variable `s′`: X(t) becomes X(s′), t s′.
    iv = get_iv(rs)
    s = _integration_variable()
    at_s = Dict{SymbolicT, SymbolicT}(x => _at_time(x, iv, s) for x in species(rs))
    at_s[iv] = s
    isvar = _unknown_predicate(rs)
    tex(x) = _latex_ordered(x, isvar; kwargs...)
    integral_tail = " \\, \\mathrm{d}" * tex(s) * " \\right)"

    _, ispcs = get_indep_sts(rs)
    crl = get_combinatoric_ratelaws(rs)
    terms = [String[] for _ in ispcs]
    species_to_idx = Dict(x => i for (i, x) in enumerate(ispcs))
    for (j, rx) in enumerate(rxs)
        rl = jumpratelaw(rx; combinatoric_ratelaw = crl, expand_catalyst_funs = false)
        Rj = "\\mathcal{R}_{$j}\\left( \\int_{0}^{" * tex(iv) * "} " *
             tex(substitute(rl, at_s)) * integral_tail
        for (spec, stoich) in rx.netstoich
            drop_dynamics(spec) && continue
            push!(terms[species_to_idx[spec]], _coef_summand(stoich, Rj, tex; kwargs...))
        end
    end
    lines = [tex(x) * " &= " * tex(_at_time(x, iv, 0)) * join(terms[i])
             for (i, x) in enumerate(ispcs)]
    return _latex_align(lines; starred, double_linebreak, math_delimiters)
end

# Substitutes each unknown by a variable that displays without its time argument.
function _latex_subber(rs)
    syms = unknowns(rs)
    isempty(syms) && return identity
    return SymbolicUtils.Substituter{true}([s => processsym(s) for s in syms],
        SymbolicUtils.default_substitute_filter)
end

# The jump form's integration variable, `s′`, displayed as `s'`: a name that models are
# unlikely to use.
_integration_variable() = _set_latex_wrapper(unwrap(Symbolics.variable(:s′)), _ -> "s'")

# `x` with the independent variable `iv` replaced by `τ`: `X(t)` becomes `X(τ)`, an indexed
# species `X(t)[i]` becomes `X(τ)[i]`, and `X(t, x)` becomes `X(τ, x)`.
function _at_time(x, iv, τ)
    isequal(x, iv) && return τ
    iscall(x) || return x
    return Symbolics.operation(x)(map(a -> _at_time(a, iv, τ), Symbolics.arguments(x))...)
end

# One SDE: "\mathrm{d}X &= drift \mathrm{d}t + noise terms", with `t` the independent
# variable `ivtex` and the drift in parentheses only if it is a sum, leaving out a zero
# drift, and "0" if there is nothing.
function _sde_line(x, drift, noise, tex, ivtex)
    dt = "\\mathrm{d}" * ivtex
    summands = _iszero(drift) ? noise :
               vcat(_coef_summand(unwrap(drift), dt, tex; sep = " \\, "), noise)
    return "\\mathrm{d}" * x * " &= " * _join_summands(summands)
end

# Joins summands that each start with " + " or " - ", dropping the first one's " + ".
function _join_summands(summands)
    isempty(summands) && return "0"
    str = join(summands)
    return startswith(str, " + ") ? str[4:end] : "-" * str[4:end]
end

# The summand `coef term`, starting with " + " or " - ". A numeric coefficient gives the sign
# and is left out when its magnitude is one; a symbolic sum is put in parentheses, and a
# negated symbolic coefficient gives the sign. A symbolic constant (e.g. the 1 that
# `linear_expansion` gives) counts as a number. `sep` goes between a coefficient and `term`.
function _coef_summand(coef, term, tex; sep = " ", kwargs...)
    (coef isa SymbolicT && SymbolicUtils.isconst(coef)) &&
        (coef = SymbolicUtils.unwrap_const(coef))
    if !(coef isa SymbolicT)
        mag = abs(coef)
        str = isone(mag) ? term : latexraw(mag; kwargs...) * sep * term
        return (coef < zero(coef) ? " - " : " + ") * str
    end
    cstr = tex(coef)
    if iscall(coef) && Symbolics.operation(coef) === (+)
        return " + \\left( " * cstr * " \\right)" * sep * term
    end
    startswith(lstrip(cstr), "-") && return " - " * tex(-coef) * sep * term
    return " + " * cstr * sep * term
end

function _no_equations_latex(rs)
    return Latexify.LaTeXString("ReactionSystem $(MT.nameof(rs)) has no reactions or equations.")
end

# Puts equation lines (each with one `&`) in an align environment, as `chemical_arrows` does.
function _latex_align(lines; starred, double_linebreak, math_delimiters)
    eol = double_linebreak ? "\\\\\\\\\n" : "\\\\\n"
    str = starred ? "\\begin{align*}\n" : "\\begin{align}\n"
    str *= join(lines, " $eol") * "  \n "
    str *= starred ? "\\end{align*}\n" : "\\end{align}\n"
    str = _tight_calls(_upright_names(str))
    math_delimiters && (str = "\$\$\n" * str * "\$\$")
    return Latexify.LaTeXString(str)
end

# Latexify writes multi-character names as `\mathtt{name}` (and a function call's name, such
# as that of X2(t), as `\mathrm{\mathtt{X2}}`). They are made upright, `\mathrm{name}`, so
# that a name like `kp` or `v0` can't be read as a product, while a single (Greek) letter,
# with or without a subscript, is left in italics.
function _upright_names(str)
    out = IOBuffer()
    i = firstindex(str)
    while (r = findnext("\\mathtt{", str, i)) !== nothing
        close = _closing_brace(str, last(r))
        name = str[nextind(str, last(r)):prevind(str, close)]
        print(out, str[i:prevind(str, first(r))])
        print(out, occursin(LATEX_ITALIC_NAME, name) ? name : "\\mathrm{" * name * "}")
        i = nextind(str, close)
    end
    print(out, str[i:end])
    str = String(take!(out))
    return replace(str, r"\\mathrm\{\\mathrm\{([^{}]*)\}\}" => s"\\mathrm{\1}")
end

# A single letter or Greek letter, possibly with subscripts and superscripts (from unicode
# ones, e.g. `X₁₂`, `k₋` or `X²`).
const LATEX_ITALIC_NAME = r"^(\\[A-Za-z]+|[A-Za-z])([_^](\{[^{}]*\}|[^{}\s\\]))*$"

# Latexify writes a call's arguments as `X\left( t \right)`, and TeX puts a thin space
# between a name and a `\left( … \right)` group, giving `X (t)`. `\mathopen{}` removes that
# gap. A `\left(` that starts a factor of a product follows a space, so is left alone.
_tight_calls(str) = replace(str, r"([A-Za-z0-9}'])\\left\(" => s"\1\\mathopen{}\\left(")

# The index of the `}` that closes the `{` at index `open` of `str`.
function _closing_brace(str, open)
    depth = 0
    i = open
    while i <= lastindex(str)
        str[i] == '{' && (depth += 1)
        str[i] == '}' && (depth -= 1) == 0 && return i
        i = nextind(str, i)
    end
    error("Unbalanced braces in LaTeX output: $str")
end

### Factor Order ###

# Renders `x` with the factors of each product in one order: the numeric coefficient, then
# the factors without unknowns (parameters, the independent variable), then those with
# unknowns, so that rate constants come before species. Symbolics orders the factors of a
# product by name instead (with powers first). So the arithmetic (+, *, /, ^) and a few
# elementary functions are turned into an `Expr` here, and the rest (variables, parameters,
# other function calls, numbers) is rendered by Symbolics. `isvar(f)` tells whether the
# factor `f` contains an unknown.
function _latex_ordered(x, isvar; kwargs...)
    x isa SymbolicT || return latexraw(x; kwargs...)
    return latexraw(_ordered_expr(x, isvar; kwargs...); kwargs...)
end

# Whether an expression contains an unknown of `rs`, by name, so that it holds whether the
# unknowns appear as X(t), as the variables `processsym` gives, or at another time.
function _unknown_predicate(rs)
    names = Set(MT.getname(u) for u in unknowns(rs))
    return f -> any(v -> MT.getname(v) in names, Symbolics.get_variables(f))
end

# A real number constant's value (other constants are rendered by Symbolics), or `nothing`.
_latex_number(x) = (v = Symbolics.value(x); v isa Real ? v : nothing)

const LATEX_ORDERED_FUNCTIONS = (sqrt, exp, log, abs, sin, cos, tan, sinh, cosh, tanh)

# A number as Symbolics renders it (e.g. 1//2 as \frac{1}{2}); integers as they are.
function _number_leaf(c; kwargs...)
    c isa Integer && return c
    c < 0 && return Expr(:call, :-, _number_leaf(-c; kwargs...))
    return isone(c) ? 1 : LaTeXString(latexraw(Num(c); kwargs...))
end

# The product of `factors`, with a thin space on either side of a multi-letter name, so that
# e.g. `k1 X2` doesn't run together as `k1X2`.
function _product_expr(factors)
    isempty(factors) && return 1
    length(factors) == 1 && return only(factors)
    ex = Expr(:call, :*, first(factors))
    for i in 2:length(factors)
        if _is_word(factors[i - 1]) || _is_word(factors[i])
            push!(ex.args, LATEX_THIN_SPACE)
        end
        push!(ex.args, factors[i])
    end
    return ex
end

const LATEX_THIN_SPACE = LaTeXString("\\,")

# Whether a rendered factor is a multi-letter name (or a power of one), which
# `_upright_names` will make upright.
function _is_word(f)
    Meta.isexpr(f, :call) && f.args[1] === :^ && return _is_word(f.args[2])
    f isa AbstractString || return false
    return any(m -> !occursin(LATEX_ITALIC_NAME, m.captures[1]),
        eachmatch(r"\\mathtt\{([^{}]*)\}", f))
end

# Whether an ordered term is negated: a negative number or `-(…)`.
function _is_negated(ex)
    (ex isa Real && ex < 0) && return true
    return Meta.isexpr(ex, :call) && ex.args[1] === :- && length(ex.args) == 2
end

function _is_negative_power(x)
    iscall(x) && Symbolics.operation(x) === (^) || return false
    e = _latex_number(Symbolics.arguments(x)[2])
    return e !== nothing && e < 0
end

function _ordered_expr(x, isvar; kwargs...)
    x = unwrap(x)
    c = _latex_number(x)
    c === nothing || return _number_leaf(c; kwargs...)
    leaf(y) = LaTeXString(latexraw(y; kwargs...))
    iscall(x) || return leaf(x)
    op = Symbolics.operation(x)
    args = sorted_arguments(x)
    ordered(y) = _ordered_expr(y, isvar; kwargs...)
    if op === (+)
        terms = Any[ordered(a) for a in args]
        # start with a positive term when there is one: X - 1 rather than -1 + X
        k = findfirst(!_is_negated, terms)
        (k === nothing || k == 1) || pushfirst!(terms, popat!(terms, k))
        return Expr(:call, :+, terms...)
    elseif op === (*)
        coef = 1
        params, vars, var_sums, denom = Any[], Any[], Any[], Any[]
        for a in map(unwrap, args)
            ac = _latex_number(a)
            if ac !== nothing
                coef *= ac
            elseif _is_negative_power(a)
                b, e = Symbolics.arguments(a)
                push!(denom, ordered(b^(-_latex_number(e))))
            elseif isvar(a) && iscall(a) && Symbolics.operation(a) === (+)
                # sums with unknowns after the other factors: X (X - 1), not (X - 1) X
                push!(var_sums, ordered(a))
            else
                group = isvar(a) ? vars : params
                # as Symbolics does, powers first within each group
                if iscall(a) && Symbolics.operation(a) === (^)
                    pushfirst!(group, ordered(a))
                else
                    push!(group, ordered(a))
                end
            end
        end
        coef_factor = isone(abs(coef)) ? Any[] : Any[_number_leaf(abs(coef); kwargs...)]
        num = _product_expr(vcat(coef_factor, params, vars, var_sums))
        ex = isempty(denom) ? num : Expr(:call, :/, num, _product_expr(denom))
        return coef < 0 ? Expr(:call, :-, ex) : ex
    elseif op === (/)
        num, den = ordered(args[1]), ordered(args[2])
        # a negated numerator goes in front of the fraction
        if Meta.isexpr(num, :call) && num.args[1] === :- && length(num.args) == 2
            return Expr(:call, :-, Expr(:call, :/, num.args[2], den))
        end
        return Expr(:call, :/, num, den)
    elseif op === (^)
        base, e = args
        if _is_negative_power(x)
            return Expr(:call, :/, 1, ordered(base^(-_latex_number(e))))
        end
        # parenthesize a call (e.g. X(t) or a function) used as a base, as Symbolics does
        b = (iscall(base) &&
             !(Symbolics.operation(base) in (+, *, /, ^, LATEX_ORDERED_FUNCTIONS...))) ?
            LaTeXString("\\left( " * leaf(base) * " \\right)") : ordered(base)
        return Expr(:call, :^, b, ordered(e))
    elseif op in LATEX_ORDERED_FUNCTIONS
        # functions that Latexify renders itself (as Symbolics would), with their arguments
        # ordered too
        return Expr(:call, nameof(op), map(ordered, args)...)
    end
    return leaf(x)
end

### Utility ###

function any_nonrx_subsys(rn::MT.AbstractSystem)
    !(rn isa ReactionSystem) && (return true)
    for subsys in get_systems(rn)
        any_nonrx_subsys(subsys) && (return true)
    end
    false
end

function make_stoich_str(spec, stoich, subber; mathrm = true, kwargs...)
    if mathrm
        prestr = "\\mathrm{"
        poststr = "}"
    else
        prestr = ""
        poststr = ""
    end

    if isequal(stoich, one(stoich))
        prestr * latexraw(subber(spec); kwargs...) * poststr
    else
        if (stoich isa SymbolicT) && iscall(stoich)
            LaTeXString("(") *
            latexraw(subber(stoich); kwargs...) *
            LaTeXString(")") *
            prestr * latexraw(subber(spec); kwargs...) * poststr
        else
            latexraw(subber(stoich); kwargs...) * LaTeXString(" ") *
            prestr * latexraw(subber(spec); kwargs...) * poststr
        end
    end
end
