#Generates the expression which is then converted into a function which generates polynomials for given parameters. Parameters not given can be given at a later stage, but parameters in the exponential must be given here.
function get_equilibration(params::Vector{Symbol},reactants::OrderedDict{Symbol,Int},f_expr::Vector{Expr})
    func_body = Expr(:block)
    push!(func_body.args,:(@polyvar internal___polyvar___x[1:$(length(reactants))]))
    push!(func_body.args,:(@polyvar internal___polyvar___p[1:$(length(params))]))
    push!(func_body.args,:([]))
    foreach(poly->push!(func_body.args[3].args,recursive_replace!(poly,(reactants,:internal___polyvar___x))), deepcopy(f_expr))
    func_expr = :((;TO___BE___REMOVED=to___be___removed) -> $func_body)
    foreach(i -> push!(func_expr.args[1].args[1].args,Expr(:kw,params[i],:(internal___polyvar___p[$i]))), 1:length(params))
    deleteat!(func_expr.args[1].args[1].args,1)
    return func_expr
end

#Adds information about some fixed concentration to the network. Macro for simplification
macro fixed_concentration(reaction_network,fixed_conc...)
    func_expr = :(add_fixed_concentration($reaction_network))
    foreach(fc -> push!(func_expr.args,fc),fixed_conc)
    return func_expr
end

#Function which does the actual work of adding the fixed concentration information. Can be called directly by inputting polynomials.
function add_fixed_concentration(reaction_network::DiffEqBase.AbstractReactionNetwork,fixed_conc::Polynomial...)
    check_polynomial(reaction_network)
    replaced = keys(reaction_network.fixed_concentrations)
    for fc in fixed_concentrations
        intersection = intersect(setdiff(reaction_network.syms,replaced),Symbil.(variables(fc)))
        (length(intersection)==0) && (@warn "Unable to replace a polynomial"; continue;)
        next_replace = intersection[1]
        push!(replaced,next_replace)
        push!(reaction_network.fixed_concentrations,next_replace=>fc)
    end
    (reaction_network.equilibratium_polynomial==nothing) && return
    foreach(sym -> reaction_network.equilibratium_polynomial[findfirst(reaction_network.syms.==sym)] = reaction_network.fixed_concentrations[sym], keys(reaction_network.fixed_concentrations))
end

#
function fix_parameters(reaction_network::DiffEqBase.AbstractReactionNetwork;kwargs...)
    check_polynomial(reaction_network)
    reaction_network.equilibratium_polynomial = [pol.num for pol in reaction_network.make_polynomial(;kwargs...)]
    foreach(sym -> reaction_network.equilibratium_polynomial[findfirst(reaction_network.syms.==sym)] = reaction_network.fixed_concentrations[sym], keys(reaction_network.fixed_concentrations))
end

#Macro running the HC template function.
macro make_hc_template(reaction_network)
    return Expr(:escape, quote
        internal___var___p___template = randn(ComplexF64, length($(reaction_network).params))
        internal___var___f___template = DynamicPolynomials.subs.($(reaction_network).equilibratium_polynomial, Ref(internal___polyvar___p => internal___var___p___template))
        internal___var___result___template = HomotopyContinuation.solve(internal___var___f___template, report_progress=false)
        $(reaction_network).homotopy_continuation_template = (internal___var___p___template,solutions(internal___var___result___template))
    end)
end

#Solves the system once using ranomd parameters. Saves the solution as a template to be used for further solving.
function make_hc_template(reaction_network::DiffEqBase.AbstractReactionNetwork)
    check_polynomial(reaction_network)
    p_template = randn(ComplexF64, length(reaction_network.params))
    f_template = DynamicPolynomials.subs.(reaction_network.equilibratium_polynomial, Ref(internal___polyvar___p => p_template))
    result_template = HomotopyContinuation.solve(f_template, report_progress=false)
    reaction_network.homotopy_continuation_template = (p_template,solutions(result_template))
end

#
function make_poly_system()
    try
        equilibratium_polynomial =
        return true
    catch
        return false
    end
end

#
function steady_states(reaction_network::DiffEqBase.AbstractReactionNetwork,pars::Vector{Float64})
    (reaction_network.homotopy_continuation_template==nothing) ? make_hc_template(reaction_network) : check_polynomial(reaction_network)
    result = HomotopyContinuation.solve(reaction_network.equilibratium_polynomial, reaction_network.homotopy_continuation_template[2], parameters=internal___polyvar___p, p₁=reaction_network.homotopy_continuation_template[1], p₀=Complex{Float64}.(pars))
    filter(realsolutions(result)) do x
            all(xᵢ -> xᵢ ≥ -0.001, x)
    end
    return result
end

#
function stability(solution::Vector{Float64},reaction_network::DiffEqBase.AbstractReactionNetwork,pars::Vector{Float64})

end

function check_polynomial(reaction_network::DiffEqBase.AbstractReactionNetwork)
    (!reaction_network.is_polynomial_system) && (error("This reaction network does not correspond to a polynomial system. Some of the reaction rate must contain non polynomial terms."))
end
