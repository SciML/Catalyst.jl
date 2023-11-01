# File containing code which I am currently unsure which file in which it belongs. Long-term functions here should be moved elsewhere.


# If u0s are not given while conservation laws are present, throws an error.
# Used in HomotopyContinuation and BifurcationKit extensions.
function conservationlaw_errorcheck(rs, pre_varmap)
    vars_with_vals = Set(p[1] for p in pre_varmap)
    any(s -> s in vars_with_vals, species(rs)) && return
    isempty(conservedequations(rs)) || 
        error("The system has conservation laws but initial conditions were not provided for some species.")
end