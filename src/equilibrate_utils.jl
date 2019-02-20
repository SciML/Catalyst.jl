#Collection of functions relating to finding fixed points of reaction networks, as well as the creation of bifurcation diagrams.

### Functions used during the construction of the reaction network. ###

#Generates the expression which is then converted into a function which generates polynomials for given parameters. Parameters not given can be given at a later stage, but parameters in the exponential must be given here. Also checks whenever the system is polynomial.
function get_equilibration(params::Vector{Symbol},reactants::OrderedDict{Symbol,Int},f_expr::Vector{Expr})
    func_body = Expr(:block)
    push!(func_body.args,:([]))
    foreach(poly->push!(func_body.args[1].args,recursive_replace!(poly,(reactants,:internal___polyvar___x))), deepcopy(f_expr))
    func_expr = :((;TO___BE___REMOVED=to___be___removed) -> $(deepcopy(func_body)))
    foreach(i -> push!(func_expr.args[1].args[1].args,Expr(:kw,params[i],:(internal___polyvar___p[$i]))), 1:length(params))
    deleteat!(func_expr.args[1].args[1].args,1)
    push!(func_expr.args[1].args,:internal___polyvar___x)
    is_pol_func = :((internal___polyvar___x,internal___polyvar___dummy=1) -> $(deepcopy(func_body)))
    foreach(p -> push!(is_pol_func.args[1].args,:($p = 1)),params)
    is_pol = quote try
        $(is_pol_func)(internal___polyvar___x)
        true
    catch
        false
    end end
    return (func_expr, is_pol)
end

#Manages post creation modification to a reaction network with respect to utility needed for equilibrium calculations. Attempts to make an equilibrium polynomial.
function manage_equilibrium_functionality!(reaction_network::DiffEqBase.AbstractReactionNetwork)
    reaction_network.is_polynomial_system && (try fix_parameters(bss) catch end)
end



## Functions required for the preparation to use homotopy continuation to find fixed points.

#Calls the internal make_polynomial function to create a polynomial. Used to model creation (will succced unless there is parameters in the exponents). Input include list of parameters to fix. All parameters occuring as exponent must be here (none other is required at this stage). Will insert fixed concentrations if such exists.
function fix_parameters(reaction_network::DiffEqBase.AbstractReactionNetwork;kwargs...)
    check_polynomial(reaction_network)
    reaction_network.equilibratium_polynomial = reaction_network.make_polynomial(reaction_network.polyvars_vars;kwargs...)
    !(typeof(reaction_network.equilibratium_polynomial[1])<:Polynomial) && (reaction_network.equilibratium_polynomial = map(pol->pol.num,reaction_network.equilibratium_polynomial))
    foreach(sym -> reaction_network.equilibratium_polynomial[findfirst(reaction_network.syms.==sym)] = reaction_network.fixed_concentrations[sym], keys(reaction_network.fixed_concentrations))
end

#In some networks some linear combinations of concentrations remain fixed. This supplies this information directly to the network.
macro fixed_concentration(reaction_network,fixed_conc...)
    func_expr = Expr(:escape,:(internal___fix___concentrations($reaction_network)))
    foreach(fc -> push!(func_expr.args[1].args,recursive_replace_vars!(balance_poly(fc), reaction_network)),fixed_conc)
    return func_expr
end
#Changes a polynomial expression of the form p(x) = q(x) to p(x) - q(x).
function balance_poly(poly::Expr)
    (poly.head != :call) && (return :($(poly.args[1])-$(poly.args[2])))
    return poly
end
#Complicated function which replaces the named varriables in the expression supplied to @fix_concentration with the corresponding polyvar stored in the reaction network.
function recursive_replace_vars!(expr::Any, rn::Symbol)
    if (typeof(expr) == Symbol)&&(expr!=:+)&&(expr!=:-)&&(expr!=:*)&&(expr!=:^)&&(expr!=:/)
        return :(in($(QuoteNode(expr)),$(rn).syms) ? $(rn).polyvars_vars[$(rn).syms_to_ints[$(QuoteNode(expr))]] : $(rn).polyvars_params[$(rn).params_to_ints[$(QuoteNode(expr))]])
    elseif typeof(expr) == Expr
        foreach(i -> expr.args[i] = recursive_replace_vars!(expr.args[i], rn), 1:length(expr.args))
    end
    return expr
end
#Function which does the actual work of adding the fixed concentration information. Not meant to be called directly, but is called through the fixed_concentration macro.
function internal___fix___concentrations(reaction_network::DiffEqBase.AbstractReactionNetwork,fixed_concentrations::Polynomial...)
    check_polynomial(reaction_network)
    replaced = Set(keys(reaction_network.fixed_concentrations))
    for fc in fixed_concentrations
        vars_in_fc = []
        foreach(v -> in(v,reaction_network.polyvars_vars) && push!(vars_in_fc,reaction_network.syms[findfirst(v.==reaction_network.polyvars_vars)]), variables(fc))
        intersection = intersect(setdiff(reaction_network.syms,replaced),vars_in_fc)
        (length(intersection)==0) && (@warn "Unable to replace a polynomial"; continue;)
        next_replace = intersection[1]
        push!(replaced,next_replace)
        push!(reaction_network.fixed_concentrations,next_replace=>fc)
    end
    (reaction_network.equilibratium_polynomial==nothing) && return
    foreach(sym -> reaction_network.equilibratium_polynomial[findfirst(reaction_network.syms.==sym)] = reaction_network.fixed_concentrations[sym], keys(reaction_network.fixed_concentrations))
end

#Can be called separatly, but will otherwise be called first time a steady state is to be found. Solves the system once using random parameters. Saves the solution as a template to be used for further solving.
function make_hc_template(reaction_network::DiffEqBase.AbstractReactionNetwork)
    check_polynomial(reaction_network)
    p_template = randn(ComplexF64, length(reaction_network.params))
    f_template = DynamicPolynomials.subs.(reaction_network.equilibratium_polynomial, Ref(reaction_network.polyvars_params => p_template))
    result_template = HomotopyContinuation.solve(f_template, report_progress=false)
    reaction_network.homotopy_continuation_template = (p_template,solutions(result_template))
end
#Macro running the HC template function.
macro make_hc_template(reaction_network)
    return Expr(:escape,:(make_hc_template($reaction_network)))
end

#Checks that the reaction network is a polynomial system.
check_is_polynomial(reaction_network::DiffEqBase.AbstractReactionNetwork) = (!reaction_network.is_polynomial_system) && (error("This reaction network does not correspond to a polynomial system. Some of the reaction rate must contain non polynomial terms."))
check_exists_polynomial(reaction_network::DiffEqBase.AbstractReactionNetwork) = (reaction_network.equilibratium_polynomial==nothing) && (error("No equilibrium polynomial have been created. Please use the @fix_concentration macro to do so."))


### Functions for finding single steady states of fixed points, and for analysing their stability..

#
function steady_states(reaction_network::DiffEqBase.AbstractReactionNetwork,params::Vector{Float64})
    (reaction_network.homotopy_continuation_template==nothing) ? make_hc_template(reaction_network) : check_polynomial(reaction_network)
    result = HomotopyContinuation.solve(reaction_network.equilibratium_polynomial, reaction_network.homotopy_continuation_template[2], parameters=reaction_network.polyvars_params, p₁=reaction_network.homotopy_continuation_template[1], p₀=params)
    filter(realsolutions(result)) do x
            all(xᵢ -> xᵢ ≥ -0.001, x)
    end
end

#
function stability(solution::Vector{Float64},reaction_network::DiffEqBase.AbstractReactionNetwork,pars::Vector{Float64})

end

### Structures holding information about bifurcation diagrams. ###

### Bifurcation Diagrams ###

#
struct bifur_path
    param::Symbol
    p_vals::Vector{Float64}
    vals::Vector{Vector{Float64}}
    jac_eigenvals::Vector{Vector{ComplexF64}}
    stability_types::Vector{Int64}
    leng::Int64
end

#
struct bifur_grid
end

### Functions called to generate bifurcation diagrams, and for plotting them. ###


function bifurcations(reaction_network::DiffEqBase.AbstractReactionNetwork,params::Vector{Float64},param::Symbol,range::Tuple{Float64,Float64})
    (reaction_network.homotopy_continuation_template==nothing) ? make_hc_template(reaction_network) : check_polynomial(reaction_network)
    p1 = copy(params); p1[reaction_network.params_to_ints[param]] = range[1];
    p2 = copy(params); p2[reaction_network.params_to_ints[param]] = range[2];
    result1 = solutions(HomotopyContinuation.solve(reaction_network.equilibratium_polynomial, reaction_network.homotopy_continuation_template[2], parameters=reaction_network.polyvars_params, p₁=reaction_network.homotopy_continuation_template[1], p₀=p1))
    result2 = solutions(HomotopyContinuation.solve(reaction_network.equilibratium_polynomial, reaction_network.homotopy_continuation_template[2], parameters=reaction_network.polyvars_params, p₁=reaction_network.homotopy_continuation_template[1], p₀=p2))
    tracker1 = pathtracker_startsolutions(reaction_network.equilibratium_polynomial, parameters=reaction_network.polyvars_params, p₁=p1, p₀=p2)[1]
    tracker2 = pathtracker_startsolutions(reaction_network.equilibratium_polynomial, parameters=reaction_network.polyvars_params, p₁=p2, p₀=p1)[1]
    paths_complete = Vector{Any}()
    paths_incomplete = Vector{Any}()
    for result in result1
        path = track_solution(tracker1,result)
        if (currstatus(tracker1) == PathTrackerStatus.success)
            remove_sol!(result2,path[2][end])
            push!(paths_complete,(1. .- path[1],path[2]))
        else
            push!(paths_incomplete,(1. .-path[1],path[2]))
        end
    end
    for result in result2
        path = track_solution(tracker2,result)
        (currstatus(tracker2) == PathTrackerStatus.success)&&remove_path!(paths_incomplete,path[2][end])
        push!(paths_complete,path)
    end
    append!(paths_complete,paths_incomplete)
    return bifur_paths(positive_real_projection.(paths_complete),param,range[1],range[2],reaction_network,params)
end

function plot_bifs(bps)
    plot();
    plot_bifs!(bps)
end
function plot_bifs!(bps,val=1)
    for bp in bps
        color = stab_color(stability_type(bp.jac_eigenvals[1]))
        plot!(bp.p_vals,getindex.(bp.vals,val),color=stab_color.(bp.stability_types),label="")
    end
    plot!()
end

### Functions required for generating bifurcation diagrams. ###

function bifur_paths(paths,param,r1,r2,reaction_network,params)
    bps = Vector{bifur_path}()
    for path in paths
        (length(path[1])==0) && continue
        (abs(1-path[1][1]/path[1][end])<0.0001) && continue
        jac_eigenvals = stabilities(path[2],param,r1 .+ ((r2-r1) .* path[1]),reaction_network,params)
        push!(bps,bifur_path(param, r1 .+ ((r2-r1) .* path[1]), path[2], jac_eigenvals, stability_type.(jac_eigenvals), length(path[1])))
    end
    return bps
end

function split_bifur_path!(bp,pos)
    bp1 = bifur_path(bp.param,bp.p_vals[1:pos],bp.vals[1:pos],bp.jac_eigenvals[1:pos],pos)
    bp2 = bifur_path(bp.param,bp.p_vals[pos:end],bp.vals[pos:end],bp.jac_eigenvals[pos:end],bp.leng-pos+1)
    return (bp1,bp2)
end

function stability_type(eigenvalues)
    stab_type = 0
    (maximum(real(eigenvalues))<1e-6)&&(stab_type+=1)
    any(imag(eigenvalues).>1e-6)&&(stab_type+=2)
    return stab_type
end

function stab_color(stab_type)
    (stab_type==0) && (return :red)
    (stab_type==1) && (return :blue)
    (stab_type==2) && (return :orange)
    (stab_type==3) && (return :cyan)
end

function stabilities(vals,param,param_vals,reaction_network,params)
    stabs = Vector{Vector{Any}}()
    Jac_temp = zeros(length(vals[1]),length(vals[1]))
    for i = 1:length(vals)
        params_i = copy(params)
        params_i[reaction_network.params_to_ints[param]] = param_vals[i]
        push!(stabs,eigen(reaction_network.jac(Jac_temp,vals[i],params_i,0.)).values)
    end
    return stabs
end

function remove_sol!(results,path_fin)
    for i = length(results):-1:1
        if maximum(abs.([imag.(path_fin.-results[i])..., real.(path_fin.-results[i])...]))<0.0000001
            deleteat!(results,i)
            return
        end
    end
end

function remove_path!(paths,path_fin)
    for i = length(paths):-1:1
        if maximum(abs.([imag.(path_fin.-paths[i][2][1])..., real.(path_fin.-paths[i][2][1])...]))<0.0000001
            deleteat!(paths,i)
            return
        end
    end
end

function track_solution(tracker,sol)
    T = []; X = [];
    for (x,t) in iterator(tracker, sol)
        push!(T,t); push!(X,x);
    end
    return (T,X)
end

function positive_real_projection(track_result)
    T = []; X = [];
    for i = 1:length(track_result[1])
        if (minimum(real.(track_result[2][i]))>-0.0001)&&(maximum(abs.(imag.(track_result[2][i])))<0.0001)
            push!(T,track_result[1][i]); push!(X,real.(track_result[2][i]))
        end
    end
    return (T,X)
end
