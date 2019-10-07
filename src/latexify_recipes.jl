@latexrecipe function f(
    r::DiffEqBase.AbstractReactionNetwork; 
    bracket=false,
    noise=false,
    noise_only=false,
    noise_var=:W,
    env=:align,
)

    env := env
    ## DiffEqBio's chemical arrow notation support is not yet migrated outside 
    ## of Latexify.jl. If used, simply pass the reaction network object on and 
    ## set the appropriate kwargs.
    if env in [:arrows, :chem, :chemical, :arrow] 
        bracket := bracket
        return r
    end

    ### Get the noise terms.
    if noise || noise_only
        noise_matrix = reshape(r.g_func, length(r.syms), :)
        w = [Symbol("d$(noise_var)_$(i)(t)") for _ in 1:size(noise_matrix, 1), i in 1:size(noise_matrix,2)]
        noise_matrix = Expr.(:call, :*, noise_matrix, w)
        noise_rhs = Meta.parse.([join(noise_matrix[i,:], " + ") for i in 1:size(noise_matrix,1)])

        ### Clean out terms which are either zero or multiplied by zero.
        for i in 1:length(noise_rhs)
            filter!(x -> x != 0, noise_rhs[i].args)
            filter!(x -> hasfield(typeof(x), :args) ? x.args[1:2] != [:*, 0] : true, noise_rhs[i].args)
        end
    end

    ### Define the left-hand side
    if noise || noise_only
        if bracket
            lhs = [Meta.parse("d[$x](t)") for x in r.syms]
        else
            lhs = [Meta.parse("d$x(t)") for x in r.syms]
        end
    else
        if bracket
            lhs = [Latexify.LaTeXString("\\frac{d\\left[ $x \\right](t)}{dt}") for x in r.syms]
        else
            lhs = [Latexify.LaTeXString("\\frac{d$x(t)}{dt}") for x in r.syms]
        end
    end

    ### Define the right-hand side
    if noise_only
        rhs = noise_rhs
        separator --> " ∝& "
    elseif noise
        rhs = Expr.(:call, :+, Expr.(:call, :*, r.f_func, :dt), noise_rhs)
    else
        rhs = r.f_func
    end
    bracket && (rhs = Latexify.add_brackets(rhs, r.syms))

    return lhs, rhs
end
