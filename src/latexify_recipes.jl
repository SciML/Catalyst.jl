### Structures & Constants ###

Base.@kwdef mutable struct CatalystLatexParams
    double_linebreak::Bool = false
    starred::Bool = true
    harpoon_arrows::Bool = true
    mathjax::Bool = false
end

const LATEX_DEFS = CatalystLatexParams()

### Latexify Receipt ###

@latexrecipe function f(rs::ReactionSystem; form = :reactions, expand_functions = true)
    expand_functions && (rs = expand_registered_functions(rs))
    if form == :reactions    # Returns chemical reaction network code.
        mult_symbol --> ""
        env --> :chem
        return rs
    elseif form == :ode      # Returns ODE system code.
        mult_symbol --> ""
        return ode_model(rs)
    elseif form == :sde      # Returns SDE system code.
        mult_symbol --> ""
        return sde_model(rs)
    end
    error("Unrecognised form argument given: $form. This should be either reactions (default), :ode, or :sde.")
end

function Latexify.infer_output(env, rs::ReactionSystem, args...)
    if env in (:arrows, :chem, :chemical, :arrow)
        return chemical_arrows
    else
        error("The environment $env is not defined.")
    end
end

function processsym(s)
    args = sorted_arguments(s)
    name = MT.getname(s)
    if length(args) <= 1
        var = value(Symbolics.variable(name))
    else
        idxs = args[2:end]
        var = value(Symbolics.variable(MT.getname(s), idxs...))
    end
    return var
end

function chemical_arrows(
        rn::ReactionSystem;
        double_linebreak = LATEX_DEFS.double_linebreak,
        starred = LATEX_DEFS.starred, mathrm = true,
        mathjax = LATEX_DEFS.mathjax, kwargs...
    )
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
    if mathjax || (
            isdefined(Main, :IJulia) && Main.IJulia.inited &&
                !any(s -> occursin("VSCODE", s), collect(keys(ENV)))
        )
        str *= "\\require{mhchem} \n"
    end

    subber = SymbolicUtils.Substituter{true}([s => processsym(s) for s in species(rn)], SymbolicUtils.default_substitute_filter)

    lastidx = length(rxs)
    for (i, r) in enumerate(rxs)
        if backwards_reaction
            backwards_reaction = false
            continue
        end

        ### Expand functions to maths expressions
        rate = r.rate isa SymbolicT ? subber(r.rate) : r.rate

        ### Generate formatted string of substrates
        substrates = [
            make_stoich_str(
                    substrate[1], substrate[2], subber; mathrm,
                    kwargs...
                )
                for substrate in zip(r.substrates, r.substoich)
        ]
        isempty(substrates) && (substrates = ["\\varnothing"])

        str *= join(substrates, " + ")

        ### Generate reaction arrows
        if i + 1 <= length(rxs) && issetequal(r.products, rxs[i + 1].substrates) &&
                issetequal(r.substrates, rxs[i + 1].products)
            ### Bi-directional arrows
            rate_backwards = rxs[i + 1].rate isa SymbolicT ? subber(rxs[i + 1].rate) :
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
        products = [
            make_stoich_str(
                    product[1], product[2], subber; mathrm = true,
                    kwargs...
                )
                for product in zip(r.products, r.prodstoich)
        ]
        isempty(products) && (products = ["\\varnothing"])
        str *= join(products, " + ")
        if (
                (i == lastidx) ||
                    (((i + 1) == lastidx) && (backwards_reaction == true))
            ) &&
                isempty(nonrxs)
            str *= "  \n "
        else
            str *= " $eol"
        end
    end

    if !isempty(nonrxs)
        eqstrs = latexraw.(nonrxs)
        eqstr_list = replace.(eqstrs, "=" => "&=")
        newstr = join(eqstr_list, " $eol")
        str *= newstr
        str *= "  \n "
    end

    str *= starred ? "\\end{align*}\n" : "\\end{align}\n"

    latexstr = Latexify.LaTeXString(str)
    Latexify.COPY_TO_CLIPBOARD && clipboard(latexstr)
    return latexstr
end

### Utility ###

function any_nonrx_subsys(rn::MT.AbstractSystem)
    !(rn isa ReactionSystem) && (return true)
    for subsys in get_systems(rn)
        any_nonrx_subsys(subsys) && (return true)
    end
    return false
end

function make_stoich_str(spec, stoich, subber; mathrm = true, kwargs...)
    if mathrm
        prestr = "\\mathrm{"
        poststr = "}"
    else
        prestr = ""
        poststr = ""
    end

    return if isequal(stoich, one(stoich))
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
