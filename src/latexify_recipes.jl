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
    _has_latex_wrapper(var) || (var = _set_latex_wrapper(var, string))
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
        latexstr = Latexify.LaTeXString("ReactionSystem $(MT.nameof(rn)) has no reactions or equations.")
        Latexify.COPY_TO_CLIPBOARD && clipboard(latexstr)
        return latexstr
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
            str *= "[" * latexraw(rate_backwards; kwargs...) * "]"
            str *= "{" * latexraw(rate; kwargs...) * "} "
            backwards_reaction = true
        else
            ### Uni-directional arrows
            str *= " &\\xrightarrow{" * latexraw(rate; kwargs...) * "} "
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
        eqstrs = [latexraw(eq; kwargs...) for eq in nonrxs]
        if !show_time_arg
            iv_latex = latexraw(get_iv(rn))
            eqstrs = [replace(s, "\\left( $iv_latex \\right)" => "") for s in eqstrs]
        end
        eqstr_list = replace.(eqstrs, "=" => "&=")
        newstr = join(eqstr_list, " $eol")
        str *= newstr
        str *= "  \n "
    end

    str *= starred ? "\\end{align*}\n" : "\\end{align}\n"

    # Strip \mathtt{} wrapping so multi-character names render as plain math italic.
    str = replace(str, r"\\mathtt\{([^}]*)\}" => s"\1")

    # Optionally wrap in $$ delimiters for renderers (Documenter, MathJax) that require
    # explicit math delimiters. Not valid standard LaTeX, but needed for web rendering.
    if math_delimiters
        str = "\$\$\n" * str * "\$\$"
    end

    latexstr = Latexify.LaTeXString(str)
    Latexify.COPY_TO_CLIPBOARD && clipboard(latexstr)
    return latexstr
end

### SDE and Jump Forms ###

# The chemical Langevin equation of a `ReactionSystem`: for each species X,
#   dX = (∑ⱼ νⱼ aⱼ) dt + ∑ⱼ νⱼ √aⱼ dWⱼ,
# with aⱼ the ODE rate law of reaction j (√aⱼ times its noise scaling, if any), νⱼ its net
# stoichiometry of X, and Wⱼ independent Brownian motions, one per reaction. Non-reaction
# differential equations are shown in the same form, with the terms of any brownians they
# contain; other equations as they are.
function cle_equations(w::LatexCLE; double_linebreak = LATEX_DEFS.double_linebreak,
        starred = LATEX_DEFS.starred, math_delimiters = false, show_time_arg = false,
        mathrm = true, mathjax = LATEX_DEFS.mathjax, kwargs...)
    rs = Catalyst.flatten(w.rs)
    # `MT.poissonians` is not public in ModelingToolkitBase (`sde_model` uses it the same way).
    (isempty(MT.jumps(rs)) && isempty(MT.poissonians(rs))) ||
        error("The SDE form can not display a ReactionSystem with coupled jumps or poissonians.")
    rxs = reactions(rs)
    nonrxs = nonreactions(rs)
    (isempty(rxs) && isempty(nonrxs)) && return _no_equations_latex(w.rs)

    subber = show_time_arg ? identity : _latex_subber(rs)
    tex(x) = latexraw(subber(x); kwargs...)
    targ = show_time_arg ? "\\left( " * latexraw(get_iv(rs); kwargs...) * " \\right)" : ""

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
        dW = tex(amp) * " \\, \\mathrm{d}W_{$j}" * targ
        for (spec, stoich) in rx.netstoich
            drop_dynamics(spec) && continue
            push!(noise[species_to_idx[spec]], _coef_summand(stoich, dW, tex; kwargs...))
        end
    end
    lines = [_sde_line(tex(x), drift[i], noise[i], tex) for (i, x) in enumerate(ispcs)]

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
                dB = "\\, \\mathrm{d}" * latexraw(B; kwargs...) * targ
                push!(eqnoise, _coef_summand(unwrap(coef), dB, tex; kwargs...))
                rhs = unwrap(rest)
            end
            x = only(Symbolics.arguments(eq.lhs))
            push!(lines, _sde_line(tex(x), rhs, eqnoise, tex))
        else
            push!(lines, tex(eq.lhs) * " &= " * tex(eq.rhs))
        end
    end
    return _latex_align(lines; starred, double_linebreak, math_delimiters)
end

# The random time change representation of the jump process of a `ReactionSystem`: for each
# species X,
#   X(t) = X(0) + ∑ⱼ νⱼ Rⱼ(∫₀ᵗ aⱼ(s) ds),
# with aⱼ the jump rate law (propensity) of reaction j, νⱼ its net stoichiometry of X, and Rⱼ
# independent unit-rate Poisson processes, one per reaction. The time arguments are always
# shown, since the integrands need them.
function time_change_equations(w::LatexTimeChange;
        double_linebreak = LATEX_DEFS.double_linebreak, starred = LATEX_DEFS.starred,
        math_delimiters = false, show_time_arg = false, mathrm = true,
        mathjax = LATEX_DEFS.mathjax, kwargs...)
    rs = Catalyst.flatten(w.rs)
    (has_nonreactions(rs) || !isempty(MT.brownians(rs)) || !isempty(MT.jumps(rs))) &&
        error("The jump form can only display reactions, but this ReactionSystem has non-reaction equations, brownians, or coupled jumps.")
    rxs = reactions(rs)
    isempty(rxs) && return _no_equations_latex(w.rs)

    # The rate laws are shown at the integration variable `s`: X(t) becomes X(s), and t s.
    iv = get_iv(rs)
    s = _integration_variable(rs)
    at_s = Dict{SymbolicT, SymbolicT}(x => Symbolics.operation(x)(s) for x in species(rs))
    at_s[iv] = s
    tex(x) = latexraw(x; kwargs...)
    integral_tail = " \\, \\mathrm{d}" * tex(s) * " \\right)"

    _, ispcs = get_indep_sts(rs)
    crl = get_combinatoric_ratelaws(rs)
    terms = [String[] for _ in ispcs]
    species_to_idx = Dict(x => i for (i, x) in enumerate(ispcs))
    for (j, rx) in enumerate(rxs)
        rl = jumpratelaw(rx; combinatoric_ratelaw = crl, expand_catalyst_funs = false)
        Rj = "R_{$j}\\left( \\int_{0}^{" * tex(iv) * "} " * tex(substitute(rl, at_s)) *
             integral_tail
        for (spec, stoich) in rx.netstoich
            drop_dynamics(spec) && continue
            push!(terms[species_to_idx[spec]], _coef_summand(stoich, Rj, tex; kwargs...))
        end
    end
    lines = [tex(x) * " &= " * tex(Symbolics.operation(x)(0)) * join(terms[i])
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

# `s`, or `τ` if the system already uses the name `s`.
function _integration_variable(rs)
    names = Set(MT.getname(x) for x in vcat(unknowns(rs), parameters(rs)))
    return unwrap(Symbolics.variable((:s in names) ? :τ : :s))
end

# One SDE: "\mathrm{d}X &= \left( drift \right) \mathrm{d}t + noise terms", leaving out a zero
# drift, and "0" if there is nothing.
function _sde_line(x, drift, noise, tex)
    summands = _iszero(drift) ? noise :
               vcat(" + \\left( " * tex(drift) * " \\right) \\, \\mathrm{d}t", noise)
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
# negated symbolic coefficient gives the sign.
function _coef_summand(coef, term, tex; kwargs...)
    if !(coef isa SymbolicT)
        mag = abs(coef)
        str = isone(mag) ? term : latexraw(mag; kwargs...) * " " * term
        return (coef < zero(coef) ? " - " : " + ") * str
    end
    cstr = tex(coef)
    if iscall(coef) && Symbolics.operation(coef) === (+)
        return " + \\left( " * cstr * " \\right) " * term
    end
    startswith(lstrip(cstr), "-") && return " - " * tex(-coef) * " " * term
    return " + " * cstr * " " * term
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
    # Strip \mathtt{} wrapping so multi-character names render as plain math italic.
    str = replace(str, r"\\mathtt\{([^}]*)\}" => s"\1")
    math_delimiters && (str = "\$\$\n" * str * "\$\$")
    return Latexify.LaTeXString(str)
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
