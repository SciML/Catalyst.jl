#Collection of functions relating to finding fixed points of reaction networks, as well as the creation of bifurcation diagrams.

### Functions used during the construction of the reaction network. ###

#Equilibrate content structure, with all the stuff required for homotopy continuation base equilibration calculations.
mutable struct EquilibrateContent
    make_polynomial::Function
    constraints::Vector{Polynomial{true,Float64}}
    homotopy_continuation_templates::Vector{NamedTuple{(:p, :sol),Tuple{Vector{Complex{Float64}},Vector{Vector{Complex{Float64}}}}}}
    equilibrium_polynomial::Union{Vector{Polynomial{true,Float64}},Nothing}
    is_polynomial_system::Bool
    polyvars::NamedTuple{(:x, :t, :p),Tuple{Vector{PolyVar{true}},PolyVar{true},Vector{PolyVar{true}}}}
    function EquilibrateContent(make_polynomial,pvxs,pvt,pvps)
        is_polynomial_system = try
            make_polynomial(pvxs,pvt,pvps,internal___var___paramvals=fill(1,length(pvps)))
            true
        catch; false; end
        equilibrium_polynomial = try
            tmp_pol = make_polynomial(pvxs,pvt,pvps);
            (typeof(tmp_pol[1])<:Polynomial) ? tmp_pol : map(pol->pol.num,tmp_pol)
        catch; nothing; end
        new(make_polynomial,Vector{Polynomial{true,Float64}}(),Vector{NamedTuple{(:p, :sol),Tuple{Vector{Complex{Float64}},Vector{Vector{Complex{Float64}}}}}}(),equilibrium_polynomial,is_polynomial_system,(x=pvxs,t=pvt,p=pvps))
    end
    function EquilibrateContent(mp,cs,hcts,ep,ips,pvs) #only for use in addequi1!()
        return new(mp,cs,hcts,ep,ips,pvs)
    end
end
#Various Functions for accessing and updating an EquilibrateContent structure.
get_equi_poly(rn::DiffEqBase.AbstractReactionNetwork) = (return rn.equilibrate_content.equilibrium_polynomial)
has_equi_poly(rn::DiffEqBase.AbstractReactionNetwork) = (return rn.equilibrate_content.equilibrium_polynomial != nothing)
is_poly(rn::DiffEqBase.AbstractReactionNetwork) = (return rn.equilibrate_content.is_polynomial_system)
get_polyvars(rn::DiffEqBase.AbstractReactionNetwork) = (return rn.equilibrate_content.polyvars)
get_hc_templates(rn::DiffEqBase.AbstractReactionNetwork) = (return rn.equilibrate_content.homotopy_continuation_templates)
has_hc_templates(rn::DiffEqBase.AbstractReactionNetwork) = (return length(rn.equilibrate_content.homotopy_continuation_templates)!=0)

push_hc_template!(rn::DiffEqBase.AbstractReactionNetwork,p_template,sol_template) = (push!(rn.equilibrate_content.homotopy_continuation_templates,(p=p_template,sol=sol_template)))
push_constraint!(rn::DiffEqBase.AbstractReactionNetwork,constraint) = (push!(rn.equilibrate_content.constraints,constraint); has_equi_poly(rn) && push!(rn.equilibrate_content.equilibrium_polynomial,constraint);)
set_equi_poly!(rn::DiffEqBase.AbstractReactionNetwork,t;kwargs...) = (rn.equilibrate_content.equilibrium_polynomial = derationalise_polys(rn.equilibrate_content.make_polynomial(rn.equilibrate_content.polyvars.x,t,rn.equilibrate_content.polyvars.p;kwargs...)))
derationalise_polys(polys) = (typeof(polys[1])<:Polynomial) ? (return polys) : (return map(pol->pol.num,polys))
add_constraints!(rn::DiffEqBase.AbstractReactionNetwork) = (has_equi_poly(rn) && foreach(constraint -> push!(rn.equilibrate_content.equilibrium_polynomial,constraint), rn.equilibrate_content.constraints))
derationalise_equi_poly!(rn::DiffEqBase.AbstractReactionNetwork) = (!(typeof(rn.equilibrate_content.equilibrium_polynomial[1])<:Polynomial) && (rn.equilibrate_content.equilibrium_polynomial = map(pol->pol.num,rn.equilibrate_content.equilibrium_polynomial)))

reset_equi_poly!(rn::DiffEqBase.AbstractReactionNetwork) = (rn.equilibrate_content.equilibrium_polynomial = nothing)
reset_hc_templates!(rn::DiffEqBase.AbstractReactionNetwork) = (rn.equilibrate_content.homotopy_continuation_templates = Vector{NamedTuple{(:p, :sol),Tuple{Vector{Complex{Float64}},Vector{Vector{Complex{Float64}}}}}}())


#Generates the expression which is then converted into a function which generates polynomials for given parameters. Parameters not given can be given at a later stage, but parameters in the exponential must be given here. Also checks whenever the system is polynomial.
function get_equilibration(params::Vector{Symbol}, reactants::OrderedDict{Symbol,Int}, f_expr::Vector{ExprValues})
    func_body = Expr(:block)
    push!(func_body.args,Expr(:if,:(length(internal___var___paramvals) != 0),Expr(:block)))
    for i = 1:length(params)
        push!(func_body.args[1].args[2].args, :((internal___var___paramvals___exemption!=$i) && ($(params[i])=(isinteger(internal___var___paramvals[$i]) ? Int64(internal___var___paramvals[$i]) : internal___var___paramvals[$i]))))
    end
    push!(func_body.args,:([]))
    foreach(poly -> push!(func_body.args[2].args, recursive_replace!(poly,(reactants,:internal___polyvar___x),(OrderedDict(:t=>1),:internal___polyvar___t))), deepcopy(f_expr))
    func_expr = :((;TO___BE___REMOVED=to___be___removed) -> $(deepcopy(func_body)))
    foreach(i -> push!(func_expr.args[1].args[1].args,Expr(:kw,params[i],:(internal___polyvar___p[$i]))), 1:length(params))
    push!(func_expr.args[1].args[1].args,Expr(:kw,:internal___var___paramvals,:([])))
    push!(func_expr.args[1].args[1].args,Expr(:kw,:internal___var___paramvals___exemption,-1))
    deleteat!(func_expr.args[1].args[1].args,1)
    push!(func_expr.args[1].args,:internal___polyvar___x)
    push!(func_expr.args[1].args,:internal___polyvar___t)
    push!(func_expr.args[1].args,:internal___polyvar___p)
    return func_expr
end


### Functions required for the preparation to use homotopy continuation to find fixed points. ###

#--- Fixes parameters, required when some exponents contain parameters. --#
#Calls the internal make_polynomial function to create a polynomial. Used to model creation (will succced unless there is parameters in the exponents). Input include list of parameters to fix. All parameters occuring as exponent must be here (none other is required at this stage). Will insert fixed concentrations if such exists.
function fix_parameters(rn::DiffEqBase.AbstractReactionNetwork, params=Vector{Number}()::Vector{Number}; t=get_polyvars(rn).t, full_vector_exemption=-1, kwargs...)
    check_is_polynomial(rn)
    set_equi_poly!(rn,t;internal___var___paramvals=params,internal___var___paramvals___exemption=full_vector_exemption,kwargs...)
    add_constraints!(rn)
    #derationalise_equi_poly!(rn)
    return nothing
end

#--- Adds additional constraints to the system, typically a fixed concentration such as in the system X <--> Y.
#In some networks some linear combinations of concentrations remain fixed. This supplies this information directly to the network.
macro add_constraint(rn::Symbol, constraints::Expr...)
    func_expr = Expr(:escape,:(internal___add___constraint!($rn)))
    foreach(constraint -> push!(func_expr.args[1].args,recursive_replace_vars!(balance_poly(constraint), rn)),constraints)
    return func_expr
end
#In some networks some linear combinations of concentrations remain fixed. This supplies this information directly to the network.
macro add_constraints(rn::Symbol, constraints::Expr...)
    func_expr = Expr(:escape,:(internal___add___constraint!($rn)))
    foreach(constraint -> push!(func_expr.args[1].args,recursive_replace_vars!(balance_poly(constraint), rn)),MacroTools.striplines(constraints...).args)
    return func_expr
    output = Expr(:block)
    for constraint in MacroTools.striplines(constraints...).args
        push!(output.args,Expr(:escape,:(internal___add___constraint!($rn))))
        push!(output.args[end].args[1].args,recursive_replace_vars!(balance_poly(constraint), rn))
    end
    return output
end
#Changes a polynomial expression of the form p(x) = q(x) to p(x) - q(x). Used by add_constraint(s).
function balance_poly(poly::Expr)
    (poly.head != :call) && (return :($(poly.args[1])-$(poly.args[2])))
    return poly
end
#Complicated function which replaces the named varriables in the expression supplied to @fix_concentration with the corresponding polyvar stored in the reaction network.
function recursive_replace_vars!(expr::Union{ExprValues,LineNumberNode}, rn::Symbol)
    if (typeof(expr) == Symbol)&&(expr!=:+)&&(expr!=:-)&&(expr!=:*)&&(expr!=:^)&&(expr!=:/)
        return :(in($(QuoteNode(expr)),$(rn).syms) ? $(rn).equilibrate_content.polyvars.x[$(rn).syms_to_ints[$(QuoteNode(expr))]] : $(rn).equilibrate_content.polyvars.p[$(rn).params_to_ints[$(QuoteNode(expr))]])
    elseif typeof(expr) == Expr
        foreach(i -> expr.args[i] = recursive_replace_vars!(expr.args[i], rn), 1:length(expr.args))
    end
    return expr
end
#Function which does the actual work of adding the fixed concentration information. Not meant to be called directly, but is called through the fixed_concentration macro.
function internal___add___constraint!(rn::DiffEqBase.AbstractReactionNetwork, constraints::Polynomial...)
    check_is_polynomial(rn)
    foreach(constraint -> push_constraint!(rn,constraint),constraints)
    remove_replicate_polys!(rn)
end
#Replicate polynomials might be bad for HomotopyContinuation, this function removes them when they might occur.
function remove_replicate_polys!(rn::DiffEqBase.AbstractReactionNetwork)
    for i = length(get_equi_poly(rn)):-1:2, j = (i-1):-1:1
        (get_equi_poly(rn)[i]==get_equi_poly(rn)[j])&&deleteat!(get_equi_poly(rn),i)
        (get_equi_poly(rn)[i]==-get_equi_poly(rn)[j])&&deleteat!(get_equi_poly(rn),i)
    end
end

#--- Makes homotopy continuation templates to be used at a later stage when solving for fixed points ---#
#Can be called separatly, but will otherwise be called first time a steady state is to be found. Solves the system once using random parameters. Saves the solution as a template to be used for further solving.
function make_hc_template(reaction_network::DiffEqBase.AbstractReactionNetwork)
    check_is_polynomial(reaction_network)
    check_exists_polynomial(reaction_network)
    p_template = randn(ComplexF64, length(reaction_network.params))
    f_template = DynamicPolynomials.subs.([reaction_network.equilibratium_polynomial...,reaction_network.fixed_concentrations...], Ref(reaction_network.polyvars_params => p_template))
    result_template = HomotopyContinuation.solve(f_template, show_progress=false)
    reaction_network.homotopy_continuation_template = (p_template,solutions(result_template))
end
#Similar to make_hc_template, but if a template is already presence this will add another one (as opposed to make which creates one and discard previous ones). Having several templates might add some robustness.
function add_hc_template(rn::DiffEqBase.AbstractReactionNetwork)
    p_template = randn(ComplexF64, length(rn.params))
    f_template = DynamicPolynomials.subs.(get_equi_poly(rn), Ref(get_polyvars(rn).p => p_template))
    solution_template = solutions(HomotopyContinuation.solve(f_template, show_progress=false))
    if (isempty(get_hc_templates(rn))) || (length(solution_template) == length(get_hc_templates(rn)[1]))
        push_hc_template!(rn,p_template,solution_template)
    elseif length(solution_template) > length(get_hc_templates(rn)[1])
        reset_hc_templates!(rn)
        push_hc_template!(rn,p_template,solution_template)
        @warn "Tried to add another homotopy continuation template to the set. However the new template contained more solutions than the previous one(s). Previous template(s) discarded."
    else
        @warn "Tried to add another homotopy continuation template to the set. However the new template contained less solutions than the previous one(s). No new template added."
    end
end
#Macro running the HC template function.
macro make_hc_template(reaction_network::Symbol)
    return Expr(:escape,:(make_hc_template($reaction_network)))
end
#Macro running the HC template function.
macro add_hc_template(reaction_network::Symbol)
    return Expr(:escape,:(add_hc_template($reaction_network)))
end

#--- Check whenever there currently is a polynomial, ot whenever such a polynomial is even possible. Returns corresponding error messages. ---#
#Checks that the reaction network is a polynomial system.
check_is_polynomial(rn::DiffEqBase.AbstractReactionNetwork) = (!is_poly(rn)) && (error("This reaction network does not correspond to a polynomial system. Some of the reaction rate must contain non polynomial terms."))
#Checks that a polynomial have been created for the reaction network.
check_exists_polynomial(reaction_network::DiffEqBase.AbstractReactionNetwork) = (reaction_network.equilibratium_polynomial==nothing) && (error("No equilibrium polynomial have been created. Please use the @fix_concentration macro to do so."))



### Functions for finding single steady states of fixed points, and for analysing their stability. ###

#Finds the steady states of a reaction network
function steady_states(rn::DiffEqBase.AbstractReactionNetwork, p=Vector{Float64}()::Vector{Float64}; solver=HcSteadyStateSolver::Function)
    return solver(rn,p)
end

#Finds steady states of a system using homotopy continuation.
function HcSteadyStateSolver(rn::DiffEqBase.AbstractReactionNetwork,p::Vector{Float64})
    using_temp_poly  = initialise_solver!(rn,p)
    (length(p)==0) && (return positive_real_solutions(solutions(HomotopyContinuation.solve(get_equi_poly(rn),show_progress=false))))
    result = hc_solve_at(rn,p)
    using_temp_poly && finalise_solver!(rn)
    return positive_real_solutions(result)
end
# Solves the reaction networks at the specific parameter values using Homotopy Continuation.
function hc_solve_at(rn::DiffEqBase.AbstractReactionNetwork,p::Vector{Float64})
    for hc_template in get_hc_templates(rn)
        result = solutions(HomotopyContinuation.solve(get_equi_poly(rn), hc_template.sol, parameters=get_polyvars(rn).p, p₁=hc_template.p, p₀=p, show_progress=false))
        (length(result) == length(hc_template.sol)) && (return result)
    end
    @warn "While solving the system using homotopy continuation some solutions were lost."
    best_length = 0; best_solution = [];
    for hc_template in get_hc_templates(rn)
        result = solutions(HomotopyContinuation.solve(get_equi_poly(rn), hc_template.sol, parameters=get_polyvars(rn).p, p₁=hc_template.p, p₀=p, show_progress=false))
        (length(result) > best_length) && (best_length = length(result); best_solution = result;)
    end
    return best_solution
end

#Prepares the reaction network for being solved with Homotopy Continuation, in case some initialising operations have not been done yet.
function initialise_solver!(rn::DiffEqBase.AbstractReactionNetwork, p::Vector{Float64}, bifurcation_exception_parameter::Int64=-1)
    check_is_polynomial(rn)
    using_temp_poly = !has_equi_poly(rn)
    using_temp_poly && fix_parameters(rn, p, full_vector_exemption=bifurcation_exception_parameter)
    !has_hc_templates(rn) && add_hc_template(rn)
    return using_temp_poly
end
# In case a temporary equilibrium polynomial were used, this one resets it.
function finalise_solver!(rn::DiffEqBase.AbstractReactionNetwork)
    reset_equi_poly!(rn); reset_hc_templates!(rn);
end
# Extracts the biological plausible solutions (real and non negative).
function positive_real_solutions(result)
    tmp_sol = filter(result) do y
        all(yᵢ -> abs(imag(yᵢ)) < 0.001, y)&&all(yᵢ -> real(yᵢ) ≥ -0.001, y)
    end
    return map(x->real.(x),tmp_sol)
end

#Returns a boolean which will be true if the given fixed point is stable.
function stability(solution::Vector{Float64}, params::Vector{Float64}, rn::DiffEqBase.AbstractReactionNetwork, t=0.::Float64)
    jac = zeros(length(rn.syms),length(rn.syms))
    rn.jac(jac,solution,params,t)
    return maximum(real.(eigen(jac).values))<0.
end



### Structures holding information about bifurcation diagrams. ###

### Bifurcation Diagrams Structures###

#Contains information from a single bifurcation line. A bifurcation plot would typically contain several such lines. These are the results from a method which tracks a solution through parameter space.
struct bifurcation_path
    p_vals::Vector{Float64}
    vals::Vector{Vector{Float64}}
    jac_eigenvals::Vector{Vector{ComplexF64}}
    stability_types::Vector{Int8}
    length::Int64
end
#A bifurcation diagram, contains a set of bifurcation paths.
struct bifurcation_diagram
    param::Symbol
    range::Tuple{Float64,Float64}
    paths::Vector{bifurcation_path}
end

#This is a single point in a grid like bifurcation scheme. It contains a parameter value and corresponding fixed points, as well as stability information.
struct bifurcation_point
    vals::Vector{Vector{Float64}}
    jac_eigenvals::Vector{Vector{ComplexF64}}
    stability_types::Vector{Int8}
end
#Contains a grid like bifurgation diagrams with a vector of bifurcation points.
struct bifurcation_grid
    param::Symbol
    range::AbstractRange
    grid_points::Vector{bifurcation_point}
    length::Int64
end
#Contains a grid like bifurgation diagrams over 2 parameters. The bifurcation points are contained in a matrix.
struct bifurcation_grid_2d
    param1::Symbol
    range1::AbstractRange
    param2::Symbol
    range2::AbstractRange
    grid_points::Matrix{bifurcation_point}
    size::Tuple{Int64,Int64}
end

#A 2d hybrid of grid and diagram. Contains a 1d grid of bifurcation diargam.
struct bifurcation_diagram_grid
    param1::Symbol
    range1::AbstractRange
    param2::Symbol
    range2::Tuple{Float64,Float64}
    bifurcation_diagrams::Vector{bifurcation_diagram}
    length::Int64
end


### Functions called to generate bifurcation diagrams. ###

#Generates a bifurcation diagram for a given reaction network and between two given parameter values. Does so by tracking paths

#Generates a bifurcation diagram.
function bifurcations(rn::DiffEqBase.AbstractReactionNetwork,p::Vector{Float64},param::Symbol,range::Tuple{Float64,Float64};solver=HcBifurcationSolver::Function,dp=(range[2]-range[1])/200.::Float64)
    return bifurcation_diagram(param,range,solver(rn,p,param,range,dp=dp))
end

#Generates a grid of bifurcation points, using a given steady state method.
function bifurcations_grid(rn::DiffEqBase.AbstractReactionNetwork,p::Vector{Float64},param::Symbol,range::AbstractRange;solver=HcSteadyStateSolver::Function)
    grid_points = Vector{Union{bifurcation_point,Nothing}}(fill(nothing,length(range)))
    for i = 1:length(range)
        p_i=copy(p); p_i[rn.params_to_ints[param]]=range[i];
        sol = solver(rn,p_i)
        jac_eigenvals = get_jac_eigenvals(map(s->ComplexF64.(s),sol),param,fill(range[i],length(sol)),rn,p)
        grid_points[i] = bifurcation_point(sol,jac_eigenvals,stability_type.(jac_eigenvals))
    end
    return bifurcation_grid(param,range,Vector{bifurcation_point}(grid_points),length(range))
end

#Generates a 2d grid of bifurcation points, using a given steady state method.
function bifurcations_grid_2d(rn::DiffEqBase.AbstractReactionNetwork,p::Vector{Float64},param1::Symbol,range1::AbstractRange,param2::Symbol,range2::AbstractRange;solver=HcSteadyStateSolver::Function)
    grid_points = Matrix{Union{bifurcation_point,Nothing}}(fill(nothing,length(range1),length(range2)))
    for i = 1:length(range1), j = 1:length(range2)
        p_ij = copy(p);
        p_ij[rn.params_to_ints[param1]] = range1[i];
        p_ij[rn.params_to_ints[param2]] = range2[j];
        sol = solver(rn,p_ij)
        jac_eigenvals = get_jac_eigenvals(map(s->ComplexF64.(s),sol),param1,fill(range1[i],length(sol)),rn,p)
        grid_points[i,j] = bifurcation_point(sol,jac_eigenvals,stability_type.(jac_eigenvals))
    end
    return bifurcation_grid_2d(param1,range1,param2,range2,Matrix{bifurcation_point}(grid_points),(length(range1),length(range2)))
end

#Generates a grid of bifurcation diagram.
function bifurcations_diagram_grid(rn::DiffEqBase.AbstractReactionNetwork,p::Vector{Float64},param1::Symbol,range1::AbstractRange,param2::Symbol,range2::Tuple{Float64,Float64};solver=HcBifurcationSolver::Function,dp=(range2[2]-range2[1])/200.::Float64)
    diagram_grid = Vector{Union{bifurcation_diagram,Nothing}}(fill(nothing,length(range1)))
    for i = 1:length(range1)
        p_i=copy(p); p_i[rn.params_to_ints[param1]]=range1[i];
        diagram_grid[i] = bifurcations(rn,p_i,param2,range2,solver=solver,dp=dp)
    end
    return bifurcation_diagram_grid(param1,range1,param2,range2,Vector{bifurcation_diagram}(diagram_grid),length(range1))
end


### Functions required for generating bifurcation diagrams. These should return a vector of bifurcation paths. ###

### Homotopy continuation based method for making bifurcation diagrams. Only works for polynomial systems and contains some heuretics for dealing with certain bifurcation points. ###

#--- Function for the simple HC bifurcation solver ---#
#Simple bifurcation solver which tracks the path from both sides. Might have problem with some diagrams with several bifurcation points.
function HcBifurcationSolverSimple(rn::DiffEqBase.AbstractReactionNetwork,p::Vector{Float64},param::Symbol,range::Tuple{Float64,Float64};dp=(range[2]-range[1])/200.::Float64,d_sol=0.0000001)
    using_temp_poly = initialise_solver!(rn,p,rn.params_to_ints[param])
    p1 = copy(p); p1[rn.params_to_ints[param]] = range[1];
    p2 = copy(p); p2[rn.params_to_ints[param]] = range[2];
    sol1 = hc_solve_at(rn, p1)
    sol2 = hc_solve_at(rn, p2)
    tracker1 = make_coretracker(rn,sol1,p1,p2,dp/(range[2]-range[1]))
    tracker2 = make_coretracker(rn,sol2,p2,p1,dp/(range[2]-range[1]))
    paths_complete   = Vector{NamedTuple{(:p, :x),Tuple{Vector{Float64},Vector{Vector{Complex{Float64}}}}}}()
    paths_incomplete = Vector{NamedTuple{(:p, :x),Tuple{Vector{Float64},Vector{Vector{Complex{Float64}}}}}}()
    for sol in sol1
        path = track_path(sol,tracker1,range...)
        if (currstatus(tracker1) == CoreTrackerStatus.success)
            remove_sol!(sol2,path.x[end],d_sol)
            push!(paths_complete,path)
        else
            push!(paths_incomplete,path)
        end
    end
    for sol in sol2
        path = track_path(sol,tracker2,reverse(range)...)
        (currstatus(tracker2) == CoreTrackerStatus.success) && remove_path!(paths_incomplete,path.x[end],d_sol)
        push!(paths_complete,path)
    end
    append!(paths_complete,paths_incomplete)
    using_temp_poly && finalise_solver!(rn)
    return bifurcation_paths(positive_real_projection.(paths_complete),param,range[1],range[2],rn,p)
end
#Used in the homotopy continuation bifurcation traces heuristics. If the path was traced succesfully to its end, checks if it corresponds to a solution in that end and removes that solution (it does not need to be attempted in the other direction).
function remove_sol!(results::Vector{Vector{ComplexF64}},path_fin::Vector{ComplexF64},d_sol::Float64)
    for i = length(results):-1:1
        if maximum(abs.([imag.(path_fin.-results[i])..., real.(path_fin.-results[i])...]))<d_sol
            deleteat!(results,i)
            return
        end
    end
end
#Used in the homotopy continuation bifurcation traces heuristics. Similar to the previous but removes an unfinished path.
function remove_path!(paths::Vector{NamedTuple{(:p, :x),Tuple{Vector{Float64},Vector{Vector{Complex{Float64}}}}}},path_fin::Vector{ComplexF64},d_sol::Float64)
    for i = length(paths):-1:1
        if maximum(abs.([imag.(path_fin.-paths[i][2][1])..., real.(path_fin.-paths[i][2][1])...]))<d_sol
            deleteat!(paths,i)
            return
        end
    end
end

#--- Function for the main HC bifurcation solver ---#
#Main bifurcation solver. Track paths from the first to the last parameter values. Detects whenever a solution is lost due to a complicated bifurcation. In this case it resumes the path tracking just after that bifurcation.
#Δp = how long step in parameter value to take after a bifrucation is encountered (should be smaller than the minimum expected distance between two bifurcations).
#Δx = distance between two solutions for them to be considered identical. Should larger than the distance between two solutions.
function HcBifurcationSolver(rn::DiffEqBase.AbstractReactionNetwork,p::Vector{Float64},param::Symbol,range::Tuple{Float64,Float64};dp=(range[2]-range[1])/200.::Float64,Δp=0.01::Float64,Δx=0.01::Float64)
    using_temp_poly = initialise_solver!(rn,p,rn.params_to_ints[param])
    parameters(p_val) = (p = copy(p); p[rn.params_to_ints[param]] = p_val; return p;)
    p_cur = range[1]; paths = Vector{NamedTuple{(:p, :x),Tuple{Vector{Float64},Vector{Vector{Complex{Float64}}}}}}();
    while p_cur < range[2]
        sol = hc_solve_at(rn, parameters(p_cur))
        substract_sols!(sol,p_cur,paths,Δx)
        ct1 = make_coretracker(rn,sol,parameters(p_cur),parameters(range[1]),dp/(range[2]-range[1]))
        ct2 = make_coretracker(rn,sol,parameters(p_cur),parameters(range[2]),dp/(range[2]-range[1]))
        new_paths = track_path_two_ways(sol,ct1,ct2,p_cur,range)
        paths = combine_paths(paths,new_paths,p_cur,Δx)
        p_cur = minimum(map(path->path.p[end],new_paths)) + Δp
    end
    using_temp_poly && finalise_solver!(rn)
    return bifurcation_paths(positive_real_projection.(paths),param,range[1],range[2],rn,p)
end
#Tracks a given solution in both directions, and combines the paths together.
function track_path_two_ways(start_points,coretracker1,coretracker2,p_cur,range)
    paths = Vector{NamedTuple{(:p, :x),Tuple{Array{Float64,1},Array{Array{Complex{Float64},1},1}}}}();
    for sp in start_points
        (P1,X1) = track_path(sp,coretracker1,range[1],p_cur)
        (P2,X2) = track_path(sp,coretracker2,p_cur,range[2])
        push!(paths,(p=[P1...,P2...],x=[X1...,X2...]))
    end
    return paths
end
#For a given solution at a certain parameter value, remove those solutions already present in the given paths.
function substract_sols!(solutions,p,paths,Δx)
    for path in filter(path -> path.p[end]≥p, paths)
        idx = findfirst(path.p .> p)
        p_pre = path.p[idx-1]; p_post = path.p[idx];
        f1 = (p_post-p)/(p_post-p_pre); f2 = (p_post-p)/(p_post-p_pre);
        path_value_at_p = f1*path.x[idx-1]+f2*path.x[idx]
        closest_sol = argmin(map(sol->norm(sol-path_value_at_p),solutions))
        (norm(solutions[closest_sol]-path_value_at_p) < Δx) && deleteat!(solutions,closest_sol)
    end
end
#For two sets of paths, removes the paths in the old path vector which seems to be similar to those in the new (starts in the same point).
function combine_paths(paths_old,paths_new,p,Δx)
    output_paths = copy(paths_new)
    first_vals = map(path -> path.x[1],paths_new)
    midpoint_vals = map(path -> path.x[findfirst(path.p .> p/2)],paths_new)
    for path in paths_old
        any(norm.(map(fv->fv-path.x[1],first_vals)) .< Δx) || any(norm.(map(mv->mv-path.x[findfirst(path.p .> p/2)],midpoint_vals)) .< Δx) && continue
        push!(output_paths,path)
    end
    return output_paths
end

#--- Bifurcation functions used by both solvers ---#
#Makes a coretracker from the given information
function make_coretracker(rn,sol,p1,p2,Δt)
    return coretracker(get_equi_poly(rn), sol, parameters=get_polyvars(rn).p, p₁=p1, p₀=p2, max_step_size=Δt)
end
#Tracks a path for the given coretracker, returns the path.
function track_path(solution,coretracker,p_start,p_end)
    P = Vector{Float64}(); X = Vector{Vector{ComplexF64}}();
    for (x,t) in iterator(coretracker, solution)
        push!(P,t2p(t,p_start,p_end)); push!(X,x);
    end
    return (p=P,x=X)
end
#Coverts a given t value (as used by homotopy continuation, 0<=t<=1), to the corresponding parameter value.
function t2p(t,p_start,p_end)
    return p_start + (1-t)*(p_end-p_start)
end
#For a given path, filters away all biologically unplausible values (negative and imaginary).
function positive_real_projection(track_result::NamedTuple{(:p, :x),Tuple{Array{Float64,1},Array{Array{Complex{Float64},1},1}}})
    T = Vector{Float64}(); X = Vector{Vector{ComplexF64}}();
    for i = 1:length(track_result.p)
        if (minimum(real.(track_result.x[i]))>-0.0001)&&(maximum(abs.(imag.(track_result.x[i])))<0.0001)
            push!(T,track_result.p[i]); push!(X,real.(track_result.x[i]))
        end
    end
    return (p=T,x=X)
end
#Takes a set of paths as tracked by homotopy continuation and turns them into a vector of bifurcation paths.
function bifurcation_paths(paths::Vector{NamedTuple{(:p, :x),Tuple{Array{Float64,1},Array{Array{Complex{Float64},1},1}}}},param::Symbol,r1::Number,r2::Number,reaction_network::DiffEqBase.AbstractReactionNetwork,params::Vector{Float64})
    bps = Vector{bifurcation_path}()
    for path in paths
        (length(path.p)==0) && continue
        (abs(1-path.p[1]/path.p[end])<0.0001) && continue
        jac_eigenvals = get_jac_eigenvals(path.x,param,path.p,reaction_network,params)
        push!(bps,bifurcation_path(path.p, path.x, jac_eigenvals, stability_type.(jac_eigenvals), length(path.p)))
    end
    return bps
end


#---Functions relating to the detection of stability in bifurcation diagrams -###

#Generates a number between 0 and 3 corresponing to some stability type (0=unstable, 1=stable, 2=unstable with imaginary eigenvalues, 3=stable with imaginary eigenvalues).
function stability_type(eigenvalues::Vector{Any})
    stab_type = Int8(0)
    (maximum(real(eigenvalues))<1e-6)&&(stab_type+=1)
    any(imag(eigenvalues).>1e-6)&&(stab_type+=2)
    return stab_type
end
#For a specific stability types, returns a color (to be used to distinguish the types in plots).
function stab_color(stab_type::Int8)
    (stab_type==0) && (return :red)
    (stab_type==1) && (return :blue)
    (stab_type==2) && (return :orange)
    (stab_type==3) && (return :cyan)
end
#For a set of fixed points, return a corresponding set with the eigen values of the jacobian at the specified points.
function get_jac_eigenvals(vals::Vector{Vector{ComplexF64}},param::Symbol,param_vals::Vector{Float64},reaction_network::DiffEqBase.AbstractReactionNetwork,params::Vector{Float64})
    stabs = Vector{Vector{Any}}()
    (length(vals)==0)&&(return stabs)
    Jac_temp = zeros(length(vals[1]),length(vals[1]))
    for i = 1:length(vals)
        params_i = copy(params)
        params_i[reaction_network.params_to_ints[param]] = param_vals[i]
        reaction_network.jac(Jac_temp,vals[i],params_i,0.)
        push!(stabs,eigen(Jac_temp).values)
    end
    return stabs
end
