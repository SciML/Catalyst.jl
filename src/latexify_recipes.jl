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
        cdot --> false
        env --> :chem
        return rs
    elseif form == :ode      # Returns ODE system code.
        cdot --> false
        return convert(ODESystem, rs)
    elseif form == :sde      # Returns SDE system code.
        cdot --> false
        return convert(SDESystem, rs)
    end
    error("Unrecognised form argument given: $form. This should be either reactions (default), :ode, or :sde.")
end

function Latexify.infer_output(env, rs::ReactionSystem, args...)
    env in [:arrows, :chem, :chemical, :arrow] && return chemical_arrows

    error("The environment $env is not defined.")
    latex_function = Latexify.get_latex_function(rs, args...)

    return latex_function
end

function chemical_arrows(rn::ReactionSystem; expand = true,
        double_linebreak = LATEX_DEFS.double_linebreak,
        starred = LATEX_DEFS.starred, mathrm = true,
        mathjax = LATEX_DEFS.mathjax, kwargs...)
    any_nonrx_subsys(rn) &&
        (@warn "Latexify currently ignores non-ReactionSystem subsystems. Please call `flatsys = flatten(sys)` to obtain a flattened version of your system before trying to Latexify it.")

    rxs = reactions(rn)
    nonrxs = filter(eq -> eq isa Equation, equations(rn))
    if isempty(rxs) && isempty(nonrxs)
        latexstr = Latexify.LaTeXString("ReactionSystem $(ModelingToolkit.nameof(rn)) has no reactions or equations.")
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

    subber = ModelingToolkit.substituter([s => value(Symbolics.variable(MT.getname(s)))
                                          for s in species(rn)])

    lastidx = length(rxs)
    for (i, r) in enumerate(rxs)
        if backwards_reaction
            backwards_reaction = false
            continue
        end

        ### Expand functions to maths expressions
        rate = r.rate isa Symbolic ? subber(r.rate) : r.rate
        rate = ModelingToolkit.prettify_expr(toexpr(rate))
        expand && (rate = recursive_clean!(rate))

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
            rate_backwards = ModelingToolkit.prettify_expr(toexpr(rxs[i + 1].rate))
            #rate_backwards = rxs[i+1].rate isa Symbolic ? Expr(subber(rxs[i+1].rate)) : rxs[i+1].rate
            expand && (rate_backwards = recursive_clean!(rate_backwards))
            expand && (rate_backwards = recursive_clean!(rate_backwards))
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
        if (i == lastidx) ||
           (((i + 1) == lastidx) && (backwards_reaction == true)) &&
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
    false
end

# Recursively traverses an expression and removes things like X^1, 1*X. Will not actually have any effect on the expression when used as a function, but will make it much easier to look at it for debugging, as well as if it is transformed to LaTeX code.
function recursive_clean!(expr)
    (expr isa Symbol) && (expr == :no___noise___scaling) && (return 1)
    (typeof(expr) != Expr) && (return expr)
    for i in 1:length(expr.args)
        expr.args[i] = recursive_clean!(expr.args[i])
    end
    (expr.args[1] == :^) && (expr.args[3] == 1) && (return expr.args[2])
    if expr.args[1] == :*
        in(0, expr.args) && (return 0)
        i = 1
        while (i = i + 1) <= length(expr.args)
            if (typeof(expr.args[i]) == Expr) && (expr.args[i].head == :call) &&
               (expr.args[i].args[1] == :*)
                for arg in expr.args[i].args
                    (arg != :*) && push!(expr.args, arg)
                end
            end
        end
        for i in length(expr.args):-1:2
            (typeof(expr.args[i]) == Expr) && (expr.args[i].head == :call) &&
                (expr.args[i].args[1] == :*) && deleteat!(expr.args, i)
            (expr.args[i] == 1) && deleteat!(expr.args, i)
        end
        (length(expr.args) == 2) && (return expr.args[2])                   # We have a multiplication of only one thing, return only that thing.
        (length(expr.args) == 1) && (return 1)                              # We have only * and no real arguments.
        (length(expr.args) == 3) && (expr.args[2] == -1) && return :(-$(expr.args[3]))
        (length(expr.args) == 3) && (expr.args[3] == -1) && return :(-$(expr.args[2]))
    end
    if expr.head == :call
        (expr.args[1] == :/) && (expr.args[3] == 1) && (return expr.args[2])
        (expr.args[1] == :binomial) && (expr.args[3] == 1) && return expr.args[2]
        #@isdefined($(expr.args[1])) || error("Function $(expr.args[1]) not defined.")
    end
    return expr
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
        if (stoich isa Symbolic) && iscall(stoich)
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
