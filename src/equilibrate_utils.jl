#Collection of functions relating to finding fixed points of reaction networks, as well as the creation of bifurcation diagrams.

### Functions used during the construction of the reaction network. ###

#Generates the expression which is then converted into a function which generates polynomials for given parameters. Parameters not given can be given at a later stage, but parameters in the exponential must be given here. Also checks whenever the system is polynomial.
function get_equilibration(params::Vector{Symbol}, reactants::OrderedDict{Symbol,Int}, f_expr::Vector{ExprValues})
    func_body = Expr(:block)
    push!(func_body.args,Expr(:if,:(length(internal___var___polyvector) != 0),Expr(:block)))
    foreach(i -> push!(func_body.args[1].args[2].args, :($(params[i])=internal___var___polyvector[$i])), 1:length(params))
    push!(func_body.args,:([]))
    foreach(poly -> push!(func_body.args[2].args, recursive_replace!(poly,(reactants,:internal___polyvar___x),(OrderedDict(:t=>1),:internal___polyvar___t))), deepcopy(f_expr))
    func_expr = :((;TO___BE___REMOVED=to___be___removed) -> $(deepcopy(func_body)))
    foreach(i -> push!(func_expr.args[1].args[1].args,Expr(:kw,params[i],:(internal___reaction___network.polyvars_params[$i]))), 1:length(params))
    push!(func_expr.args[1].args[1].args,Expr(:kw,:internal___var___polyvector,:([])))
    deleteat!(func_expr.args[1].args[1].args,1)
    push!(func_expr.args[1].args,:internal___reaction___network)
    push!(func_expr.args[1].args,:internal___polyvar___t)
    push!(func_expr.args[1].args,Expr(:kw,:internal___polyvar___x,:(internal___reaction___network.polyvars_vars)))
    push!(func_expr.args[1].args,Expr(:kw,:print_msg,1))
    return func_expr
end

#Manages post creation modification to a reaction network with respect to utility needed for equilibrium calculations. Attempts to make an equilibrium polynomial.
function manage_equilibrium_functionality!(reaction_network::DiffEqBase.AbstractReactionNetwork)
    reaction_network.is_polynomial_system = try
        fix_parameters(reaction_network,fill(1,length(reaction_network.params)),t=1)
        true
    catch
        false
    end
    reaction_network.is_polynomial_system && try fix_parameters(reaction_network)
    catch
        (reaction_network.equilibratium_polynomial = nothing)
    end
end

## Functions required for the preparation to use homotopy continuation to find fixed points.

#Calls the internal make_polynomial function to create a polynomial. Used to model creation (will succced unless there is parameters in the exponents). Input include list of parameters to fix. All parameters occuring as exponent must be here (none other is required at this stage). Will insert fixed concentrations if such exists.
function fix_parameters(reaction_network::DiffEqBase.AbstractReactionNetwork, params=Vector{Number}()::Vector{Number}; t=reaction_network.polyvars_t, kwargs...)
    check_is_polynomial(reaction_network)
    reaction_network.equilibratium_polynomial = reaction_network.make_polynomial(reaction_network,[t];internal___var___polyvector=params,kwargs...)
    !(typeof(reaction_network.equilibratium_polynomial[1])<:Polynomial) && (reaction_network.equilibratium_polynomial = map(pol->pol.num,reaction_network.equilibratium_polynomial))
    foreach(sym -> reaction_network.equilibratium_polynomial[findfirst(reaction_network.syms.==sym)] = reaction_network.fixed_concentrations[sym], keys(reaction_network.fixed_concentrations))
end

#In some networks some linear combinations of concentrations remain fixed. This supplies this information directly to the network.
macro add_constraint(reaction_network::Symbol, fixed_conc::Expr...)
    func_expr = Expr(:escape,:(internal___add___constraint($reaction_network)))
    foreach(fc -> push!(func_expr.args[1].args,recursive_replace_vars!(balance_poly(fc), reaction_network)),fixed_conc)
    return func_expr
end
#Changes a polynomial expression of the form p(x) = q(x) to p(x) - q(x).
function balance_poly(poly::Expr)
    (poly.head != :call) && (return :($(poly.args[1])-$(poly.args[2])))
    return poly
end
#Complicated function which replaces the named varriables in the expression supplied to @fix_concentration with the corresponding polyvar stored in the reaction network.
function recursive_replace_vars!(expr::Union{ExprValues,LineNumberNode}, rn::Symbol)
    if (typeof(expr) == Symbol)&&(expr!=:+)&&(expr!=:-)&&(expr!=:*)&&(expr!=:^)&&(expr!=:/)
        return :(in($(QuoteNode(expr)),$(rn).syms) ? $(rn).polyvars_vars[$(rn).syms_to_ints[$(QuoteNode(expr))]] : $(rn).polyvars_params[$(rn).params_to_ints[$(QuoteNode(expr))]])
    elseif typeof(expr) == Expr
        foreach(i -> expr.args[i] = recursive_replace_vars!(expr.args[i], rn), 1:length(expr.args))
    end
    return expr
end
#Function which does the actual work of adding the fixed concentration information. Not meant to be called directly, but is called through the fixed_concentration macro.
function internal___add___constraint(reaction_network::DiffEqBase.AbstractReactionNetwork, fixed_concentrations::Polynomial...)
    check_is_polynomial(reaction_network)
    for fc in fixed_concentrations
        reaction_network.fixed_concentrations = [reaction_network.fixed_concentrations...,fc]
    end
    #replaced = Set(keys(reaction_network.fixed_concentrations))
    #for fc in fixed_concentrations
    #    vars_in_fc = []
    #    foreach(v -> in(v,reaction_network.polyvars_vars) && push!(vars_in_fc,reaction_network.syms[findfirst(v.==reaction_network.polyvars_vars)]), variables(fc))
    #    intersection = intersect(setdiff(reaction_network.syms,replaced),vars_in_fc)
    #    (length(intersection)==0) && (@warn "Unable to replace a polynomial"; continue;)
    #    next_replace = intersection[1]
    #    push!(replaced,next_replace)
    #    push!(reaction_network.fixed_concentrations,next_replace=>fc)
    #end
    #(reaction_network.equilibratium_polynomial==nothing) && return
    #foreach(sym -> reaction_network.equilibratium_polynomial[findfirst(reaction_network.syms.==sym)] = reaction_network.fixed_concentrations[sym], keys(reaction_network.fixed_concentrations))
end
#In some networks some linear combinations of concentrations remain fixed. This supplies this information directly to the network.
macro add_constraints(reaction_network::Symbol, constraints::Expr...)
    output = Expr(:block)
    for constraint in MacroTools.striplines(constraints...).args
        push!(output.args,Expr(:escape,:(internal___add___constraint($reaction_network))))
        push!(output.args[end].args[1].args,recursive_replace_vars!(balance_poly(constraint), reaction_network))
    end
    return output
end

#Can be called separatly, but will otherwise be called first time a steady state is to be found. Solves the system once using random parameters. Saves the solution as a template to be used for further solving.
function make_hc_template(reaction_network::DiffEqBase.AbstractReactionNetwork)
    check_is_polynomial(reaction_network)
    check_exists_polynomial(reaction_network)
    p_template = randn(ComplexF64, length(reaction_network.params))
    f_template = DynamicPolynomials.subs.([reaction_network.equilibratium_polynomial...,reaction_network.fixed_concentrations...], Ref(reaction_network.polyvars_params => p_template))
    result_template = HomotopyContinuation.solve(f_template, show_progress=false)
    reaction_network.homotopy_continuation_template = (p_template,solutions(result_template))
end
#Macro running the HC template function.
macro make_hc_template(reaction_network::Symbol)
    return Expr(:escape,:(make_hc_template($reaction_network)))
end

#Checks that the reaction network is a polynomial system.
check_is_polynomial(reaction_network::DiffEqBase.AbstractReactionNetwork) = (!reaction_network.is_polynomial_system) && (error("This reaction network does not correspond to a polynomial system. Some of the reaction rate must contain non polynomial terms."))
#Checks that a polynomial have been created for the reaction network.
check_exists_polynomial(reaction_network::DiffEqBase.AbstractReactionNetwork) = (reaction_network.equilibratium_polynomial==nothing) && (error("No equilibrium polynomial have been created. Please use the @fix_concentration macro to do so."))
#If not avaiable, ensures a temporary polynomial is set.
function set_temporary_polynomial!(reaction_network::DiffEqBase.AbstractReactionNetwork, params::Vector{Float64})
    (!reaction_network.is_polynomial_system || reaction_network.equilibratium_polynomial!=nothing) && return
    fix_parameters(reaction_network, params)
    make_hc_template(reaction_network)
    reaction_network.has_temporary_polynomial = true
end
#If a temporary polynomial was used, removes it.
function reset_temporary_polynomial!(reaction_network::DiffEqBase.AbstractReactionNetwork)
    (reaction_network.has_temporary_polynomial) && return
    reaction_network.equilibratium_polynomial = nothing
    reaction_network.homotopy_continuation_template = nothing
end

### Functions for finding single steady states of fixed points, and for analysing their stability..

#Finds the steady states of a reaction network
function steady_states(reaction_network::DiffEqBase.AbstractReactionNetwork, params=Vector{Float64}()::Vector{Float64}, solver=HcSteadyStateSolver::Function)
    return solver(reaction_network,params)::Vector{Vector{Float64}}
end

#Finds steady states of a system using homotopy continuation.
function HcSteadyStateSolver(reaction_network::DiffEqBase.AbstractReactionNetwork,params::Vector{Float64})
    set_temporary_polynomial!(reaction_network, params)
    if length(params)==0
        result = HomotopyContinuation.solve([reaction_network.equilibratium_polynomial...,reaction_network.fixed_concentrations...], show_progress=false)
    else
        (reaction_network.homotopy_continuation_template==nothing) && make_hc_template(reaction_network)
        result = HomotopyContinuation.solve([reaction_network.equilibratium_polynomial...,reaction_network.fixed_concentrations...], reaction_network.homotopy_continuation_template[2], parameters=reaction_network.polyvars_params, p₁=reaction_network.homotopy_continuation_template[1], p₀=params, show_progress=false)
    end
    reset_temporary_polynomial!(reaction_network)
    filter(realsolutions(result)) do x
            all(xᵢ -> xᵢ ≥ -0.001, x)
    end
end

#Returns a boolean which will be true if the given fixed point is stable.
function stability(solution::Vector{Float64}, params::Vector{Float64}, reaction_network::DiffEqBase.AbstractReactionNetwork, t=0.::Float64)
    return maximum(real.(eigen(reaction_network.jac(zeros(length(reaction_network.syms),length(reaction_network.syms)),solution,params,t)).values))<0.
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
function bifurcations(reaction_network::DiffEqBase.AbstractReactionNetwork,params::Vector{Float64},param::Symbol,range::Tuple{Float64,Float64};solver=HcBifurcationSolver1::Function,stepsize=0.01::Float64)
    return bifurcation_diagram(param,range,solver(reaction_network,params,param,range,stepsize=stepsize))
end

#Generates a grid of bifurcation points, using a given steady state method.
function bifurcations_grid(reaction_network::DiffEqBase.AbstractReactionNetwork,params::Vector{Float64},param::Symbol,range1::AbstractRange;solver=HcSteadyStateSolver::Function)
    grid_points = Vector{Union{bifurcation_point,Nothing}}(fill(nothing,length(range1)))
    for i = 1:length(range1)
        params_i=copy(params); params_i[reaction_network.params_to_ints[param]]=range1[i];
        ss = solver(reaction_network,params_i)
        jac_eigenvals = get_jac_eigenvals(map(s->ComplexF64.(s),ss),param,fill(range1[i],length(ss)),reaction_network,params)
        grid_points[i] = bifurcation_point(ss,jac_eigenvals,stability_type.(jac_eigenvals))
    end
    return bifurcation_grid(param,range1,Vector{bifurcation_point}(grid_points),length(range1))
end

#Generates a 2d grid of bifurcation points, using a given steady state method.
function bifurcations_grid_2d(reaction_network::DiffEqBase.AbstractReactionNetwork,params::Vector{Float64},param1::Symbol,range1::AbstractRange,param2::Symbol,range2::AbstractRange;solver=HcSteadyStateSolver::Function)
    grid_points = Matrix{Union{bifurcation_point,Nothing}}(fill(nothing,length(range1),length(range2)))
    for i = 1:length(range1), j = 1:length(range2)
        params_ij = copy(params);
        params_ij[reaction_network.params_to_ints[param1]] = range1[i];
        params_ij[reaction_network.params_to_ints[param2]] = range2[j];
        ss = solver(reaction_network,params_ij)
        jac_eigenvals = get_jac_eigenvals(map(s->ComplexF64.(s),ss),param1,fill(range1[i],length(ss)),reaction_network,params)
        grid_points[i,j] = bifurcation_point(ss,jac_eigenvals,stability_type.(jac_eigenvals))
    end
    return bifurcation_grid_2d(param1,range1,param2,range2,Matrix{bifurcation_point}(grid_points),(length(range1),length(range2)))
end

#Generates a grid of bifurcation diagram.
function bifurcations_diagram_grid(reaction_network::DiffEqBase.AbstractReactionNetwork,params::Vector{Float64},param1::Symbol,range1::AbstractRange,param2::Symbol,range2::Tuple{Float64,Float64};solver=HcBifurcationSolver1::Function,stepsize=0.01::Float64)
    diagram_grid = Vector{Union{bifurcation_diagram,Nothing}}(fill(nothing,length(range1)))
    for i = 1:length(range1)
        params_i=copy(params); params_i[reaction_network.params_to_ints[param1]]=range1[i];
        diagram_grid[i] = bifurcations(reaction_network,params_i,param2,range2,solver=solver,stepsize=stepsize)
    end
    return bifurcation_diagram_grid(param1,range1,param2,range2,Vector{bifurcation_diagram}(diagram_grid),length(range1))
end


### Functions required for generating bifurcation diagrams. These shoudl return a vector of bifurcation paths. ###

### Homotopy continuation based method for making bifurcation diagrams. Only works for polynomial systems and contains some heuretics for dealing with certain bifurcation points. ###

#Creates a bifurcation diagram by tracking using homotopy continuation. Might have problems with certain bifurcations (which might be fixed by introducing a few heuristics).
function HcBifurcationSolver1(reaction_network::DiffEqBase.AbstractReactionNetwork,params::Vector{Float64},param::Symbol,range::Tuple{Float64,Float64};stepsize=0.01::Float64)
    (reaction_network.homotopy_continuation_template==nothing) ? make_hc_template(reaction_network) : check_is_polynomial(reaction_network)
    p1 = copy(params); p1[reaction_network.params_to_ints[param]] = range[1];
    p2 = copy(params); p2[reaction_network.params_to_ints[param]] = range[2];
    result1 = solutions(HomotopyContinuation.solve([reaction_network.equilibratium_polynomial...,reaction_network.fixed_concentrations...], reaction_network.homotopy_continuation_template[2], parameters=reaction_network.polyvars_params, p₁=reaction_network.homotopy_continuation_template[1], p₀=p1, show_progress=false))
    result2 = solutions(HomotopyContinuation.solve([reaction_network.equilibratium_polynomial...,reaction_network.fixed_concentrations...], reaction_network.homotopy_continuation_template[2], parameters=reaction_network.polyvars_params, p₁=reaction_network.homotopy_continuation_template[1], p₀=p2, show_progress=false))
    tracker1 = coretracker(reaction_network.equilibratium_polynomial, parameters=reaction_network.polyvars_params, p₁=p1, p₀=p2, max_step_size=stepsize)
    tracker2 = coretracker(reaction_network.equilibratium_polynomial, parameters=reaction_network.polyvars_params, p₁=p2, p₀=p1, max_step_size=stepsize)
    paths_complete = Vector{Tuple{Vector{Float64},Vector{Vector{ComplexF64}}}}()
    paths_incomplete = Vector{Tuple{Vector{Float64},Vector{Vector{ComplexF64}}}}()
    for result in result1
        path = track_solution(tracker1,result)
        if (currstatus(tracker1) == CoreTrackerStatus.success)
            remove_sol!(result2,path[2][end])
            push!(paths_complete,(1. .- path[1],path[2]))
        else
            push!(paths_incomplete,(1. .-path[1],path[2]))
        end
    end
    for result in result2
        path = track_solution(tracker2,result)
        (currstatus(tracker2) == CoreTrackerStatus.success)&&remove_path!(paths_incomplete,path[2][end])
        push!(paths_complete,path)
    end
    append!(paths_complete,paths_incomplete)
    reset_temporary_polynomial!(reaction_network)
    return bifurcation_paths(positive_real_projection.(paths_complete),param,range[1],range[2],reaction_network,params)
end

#Used in the homotopy continuation bifurcation traces heuristics. If the path was traced succesfully to its end, checks if it corresponds to a solution in that end and removes that solution (it does not need to be attempted in the other direction).
function remove_sol!(results::Vector{Vector{ComplexF64}},path_fin::Vector{ComplexF64})
    for i = length(results):-1:1
        if maximum(abs.([imag.(path_fin.-results[i])..., real.(path_fin.-results[i])...]))<0.0000001
            deleteat!(results,i)
            return
        end
    end
end
#Used in the homotopy continuation bifurcation traces heuristics. Similar to the previous but removes an unfinished path.
function remove_path!(paths::Vector{Tuple{Vector{Float64},Vector{Vector{ComplexF64}}}},path_fin::Vector{ComplexF64})
    for i = length(paths):-1:1
        if maximum(abs.([imag.(path_fin.-paths[i][2][1])..., real.(path_fin.-paths[i][2][1])...]))<0.0000001
            deleteat!(paths,i)
            return
        end
    end
end

#Second iteration of the HC solver. Whenever a solution is lost, the parameter is incremented slightly (to pass the singularity), after which the path is taken up again.
function HcBifurcationSolver2(reaction_network::DiffEqBase.AbstractReactionNetwork,params::Vector{Float64},param::Symbol,range::Tuple{Float64,Float64};stepsize=0.01::Float64,p_skip=0.001::Float64,thres1=0.0001::Float64,thres2=0.2::Float64)
    (reaction_network.homotopy_continuation_template==nothing) ? make_hc_template(reaction_network) : check_is_polynomial(reaction_network)
    p_end = copy(params); p_end[reaction_network.params_to_ints[param]] = range[2];
    function make_coretracker(p_start)
        pars = copy(params); pars[reaction_network.params_to_ints[param]] = p_start;
        coretracker(reaction_network.equilibratium_polynomial, parameters=reaction_network.polyvars_params, p₁=pars, p₀=p_end, max_step_size=stepsize)
    end
    function solve_system(p)
        pars = copy(params); pars[reaction_network.params_to_ints[param]] = p;
        solutions(HomotopyContinuation.solve([reaction_network.equilibratium_polynomial...,reaction_network.fixed_concentrations...], reaction_network.homotopy_continuation_template[2], parameters=reaction_network.polyvars_params, p₁=reaction_network.homotopy_continuation_template[1], p₀=pars, show_progress=false))
    end
    results = solve_system(range[1])
    paths_complete = Vector{Tuple{Vector{Float64},Vector{Vector{ComplexF64}}}}()
    breakpoints = Vector{Vector{ComplexF64}}()
    for result in results
        recursive_track_solutions!(result,range[1],paths_complete,breakpoints,make_coretracker,solve_system,range[2],p_skip=0.0001,thres1=0.0001,thres2=0.1)
    end
    reset_temporary_polynomial!(reaction_network)
    return bifurcation_paths(positive_real_projection.(paths_complete),param,range[1],range[2],reaction_network,params)
end
#Recursively creates new paths. required when one path splits into several.
function recursive_track_solutions!(start_point::Vector{ComplexF64},start_p::Float64,paths_complete::Vector{Tuple{Vector{Float64},Vector{Vector{ComplexF64}}}},breakpoints::Vector{Vector{ComplexF64}},make_coretracker::Function,solve_system::Function,p_end::Float64;p_skip=0.0001::Float64,thres1=0.0001::Float64,thres2=0.1::Float64)
    any(map(bp->maximum(abs.([start_p,start_point...].-bp))<thres1, breakpoints)) && return
    path = track_solution(make_coretracker(start_p),start_point)
    push!(paths_complete,(1. .- path[1],path[2]))
    (path[1][end]==0.0) && return
    push!(breakpoints,[path[1][end],path[2][end]...])
    p_skip = start_p + (1-path[1][end]+p_skip)*(p_end-start_p)
    skip_sols = solve_system(p_skip)
    for sol in skip_sols
        (maximum(abs.(path[2][end].-sol))>thres2) && continue
        recursive_track_solutions!(sol,p_skip,paths_complete,breakpoints,make_coretracker,solve_system,p_end,p_skip=0.0001,thres1=0.0001,thres2=0.1)
    end
end

#For a given tracker and solution, tries to track the solution from t=1 to t=0.
function track_solution(tracker::CoreTracker,sol::Vector{ComplexF64})
    T = Vector{Float64}(); X = Vector{Vector{ComplexF64}}();
    for (x,t) in iterator(tracker, sol)
        push!(T,t); push!(X,x);
    end
    return (T,X)
end
#For a given path, filters away all biologically unplausible values (negative and imaginary).
function positive_real_projection(track_result::Tuple{Vector{Float64},Vector{Vector{ComplexF64}}})
    T = Vector{Float64}(); X = Vector{Vector{ComplexF64}}();
    for i = 1:length(track_result[1])
        if (minimum(real.(track_result[2][i]))>-0.0001)&&(maximum(abs.(imag.(track_result[2][i])))<0.0001)
            push!(T,track_result[1][i]); push!(X,real.(track_result[2][i]))
        end
    end
    return (T,X)
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
function bifurcation_paths(paths::Vector{Tuple{Vector{Float64},Vector{Vector{ComplexF64}}}},param::Symbol,r1::Number,r2::Number,reaction_network::DiffEqBase.AbstractReactionNetwork,params::Vector{Float64})
    bps = Vector{bifurcation_path}()
    for path in paths
        (length(path[1])==0) && continue
        (abs(1-path[1][1]/path[1][end])<0.0001) && continue
        jac_eigenvals = get_jac_eigenvals(path[2],param,r1 .+ ((r2-r1) .* path[1]),reaction_network,params)
        push!(bps,bifurcation_path(r1 .+ ((r2-r1) .* path[1]), path[2], jac_eigenvals, stability_type.(jac_eigenvals), length(path[1])))
    end
    return bps
end
#Takes a set of paths as tracked by homotopy continuation and turns them into a vector of bifurcation paths.
function bifurcation_paths(paths::Vector{NamedTuple{(:p, :x),Tuple{Array{Float64,1},Array{Array{Complex{Float64},1},1}}}},param::Symbol,r1::Number,r2::Number,reaction_network::DiffEqBase.AbstractReactionNetwork,params::Vector{Float64})
    bps = Vector{bifurcation_path}()
    for path in paths
        (length(path.p)==0) && continue
        (abs(1-path.p[1]/path.p[end])<0.0001) && continue
        jac_eigenvals = get_jac_eigenvals(path.x,param,r1 .+ ((r2-r1) .* path.p),reaction_network,params)
        push!(bps,bifurcation_path(path.p, path.x, jac_eigenvals, stability_type.(jac_eigenvals), length(path.p)))
    end
    return bps
end

function HcBifurcationSolver3(reaction_network::DiffEqBase.AbstractReactionNetwork,params::Vector{Float64},param::Symbol,range::Tuple{Float64,Float64};stepsize=0.01::Float64,Δp=0.01::Float64,Δx=0.01::Float64)
    uses_temp_polynomial = temp_initialise_solver(reaction_network,params)
    parameters(p_val) = (p = copy(params); p[reaction_network.params_to_ints[param]] = p_val; return p;)
    p_cur = range[1]; paths = Vector{NamedTuple{(:p, :x),Tuple{Array{Float64,1},Array{Array{Complex{Float64},1},1}}}}();
    while p_cur < range[2]
        sol = solve_at(reaction_network,parameters(p_cur))
        substract_sols!(sol,p_cur,paths,Δx)
        ct1 = make_coretracker(reaction_network,sol,parameters(p_cur),parameters(range[1]),stepsize)
        ct2 = make_coretracker(reaction_network,sol,parameters(p_cur),parameters(range[2]),stepsize)
        new_paths = track_path_two_ways(sol,ct1,ct2,p_cur,range)
        paths = combine_paths(paths,new_paths,p_cur,Δx)
        println("\n",p_cur)
        for path in new_paths
            println(path.p[1]," ",path.p[end],"\t",path.x[1][1]," ",path.x[end][1])
        end
        p_cur = minimum(map(path->path.p[end],paths)) + Δp
        println(p_cur)
    end
    temp_finalise_solver(uses_temp_polynomial,reaction_network)
    return bifurcation_paths(positive_real_projection.(paths),param,range[1],range[2],reaction_network,params)
end
function solve_at(reaction_network,p)
    return solutions(HomotopyContinuation.solve([reaction_network.equilibratium_polynomial...,reaction_network.fixed_concentrations...], reaction_network.homotopy_continuation_template[2], parameters=reaction_network.polyvars_params, p₁=reaction_network.homotopy_continuation_template[1], p₀=p, show_progress=false))
end
function t2p(t::AbstractFloat,range::Tuple{Float64,Float64})
    return range[1] + t*(range[2]-range[1])
end
function make_coretracker(reaction_network,solution,p1,p2,max_step_size)
    return coretracker([reaction_network.equilibratium_polynomial...,reaction_network.fixed_concentrations...], solution, parameters=reaction_network.polyvars_params, p₁=p1, p₀=p2, max_step_size=max_step_size)
end
function track_path_two_ways(start_points,coretracker1,coretracker2,p_cur,range)
    paths = Vector{NamedTuple{(:p, :x),Tuple{Array{Float64,1},Array{Array{Complex{Float64},1},1}}}}();
    for sp in start_points
        (P1,X1) = track_path(sp,coretracker1,p_cur,range[1])
        (P2,X2) = track_path(sp,coretracker2,p_cur,range[2])
        push!(paths,(p=[reverse(P1)...,reverse(P2)...],x=[X1...,X2...]))
    end
    return compress_path!.(paths)
end
function track_path(solution,coretracker,p_start,p_end)
    P = Vector{Float64}(); X = Vector{Vector{ComplexF64}}();
    for (x,t) in iterator(coretracker, solution)
        push!(P,t2p(t,p_start,p_end)); push!(X,x);
    end
    return (p=P,x=X)
end
function compress_path!(path)
    for i = length(path.p):-1:2
        if (path.p[i]==path.p[i-1]) && (path.x[i]==path.x[i-1])
            deleteat!(path.p,i); deleteat!(path.x,i);
        end
    end
    return path
end
function t2p(t,p_start,p_end)
    return p_start + t*(p_end-p_start)
end
function substract_sols!(solutions,p,paths,Δx)
    for path in paths
        idx = findfirst(path.p .> p)
        p_pre = path.p[idx-1]; p_post = path.p[idx];
        f1 = (p_post-p)/(p_post-p_pre); f2 = (p_post-p)/(p_post-p_pre);
        path_value_at_p = f1*path.x[idx-1]+f2*path.x[idx]
        closest_sol = argmin(map(sol->norm(sol-path_value_at_p),solutions))
        (norm(solutions[closest_sol]-path_value_at_p) < Δx) && deleteat!(solutions,closest_sol)
    end
end
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
function temp_initialise_solver(reaction_network::DiffEqBase.AbstractReactionNetwork, params::Vector{Float64})
    check_is_polynomial(reaction_network)
    using_temp_poly = (reaction_network.equilibratium_polynomial==nothing)
    using_temp_poly && fix_parameters(reaction_network, params)
    (reaction_network.homotopy_continuation_template==nothing) && make_hc_template(reaction_network)
    return using_temp_poly
end
function temp_finalise_solver(uses_temporary_polynomial::Bool,reaction_network::DiffEqBase.AbstractReactionNetwork)
    !uses_temporary_polynomial && return
    reaction_network.equilibratium_polynomial = nothing
    reaction_network.homotopy_continuation_template = nothing
end



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
        push!(stabs,eigen(reaction_network.jac(Jac_temp,vals[i],params_i,0.)).values)
    end
    return stabs
end

### Functions for plotting various types of bifurcation diagrams ###
#- SHOULD BECOME A PROPER OVERLOAD OF PLOT AT SOME POINT, RATHER THAN A SEPARATE FUNCTION -#

#Normal plot overlaods into the plot!, same for all.
function bif_plot(bif_dia::Union{bifurcation_diagram,bifurcation_grid,bifurcation_grid_2d,bifurcation_diagram_grid}, var=1::Integer, args...; kwargs...)
    plot(); bif_plot!(bif_dia, var, args...; kwargs...)
end

#Plots a bifurcation diagram.
function bif_plot!(bif_dia::bifurcation_diagram, var=1::Integer, args...; kwargs...)
    foreach(bp -> plot!(bp.p_vals,getindex.(bp.vals,var), args...; label="", color=stab_color.(bp.stability_types), kwargs...), bif_dia.paths)
    plot!(args...; xlabel=bif_dia.param, kwargs...)
end

#Plots a bifurcation grid.
function bif_plot!(bif_grid::bifurcation_grid, var=1::Integer, args...; kwargs...)
    for i = 1:bif_grid.length
        for k = 1:length(bif_grid.grid_points[i].vals)
            scatter!((bif_grid.range[i],bif_grid.grid_points[i].vals[k][var]), args...; label="", color=stab_color(bif_grid.grid_points[i].stability_types[k]), kwargs...)
        end
    end
    plot!(args...; xlabel=bif_grid.param, kwargs...)
end

#Plots a two-dimensional bifurcation grid.
function bif_plot!(bif_grid_2d::bifurcation_grid_2d, var=1::Integer, args...; kwargs...)
    for i = 1:bif_grid_2d.size[1], j = 1:bif_grid_2d.size[2]
        for k = 1:length(bif_grid_2d.grid_points[i,j].vals)
            scatter!((bif_grid_2d.range1[i],bif_grid_2d.range2[j],bif_grid_2d.grid_points[i,j].vals[k][var]), args...; label="", color=stab_color(bif_grid_2d.grid_points[i,j].stability_types[k]), kwargs...)
        end
    end
    plot!(args...; xlabel=bif_grid_2d.param1, zlabel=bif_grid_2d.param2, kwargs...)
end

#Plots a grid of bifurcation diagrams.
function bif_plot!(bif_dia_grid::bifurcation_diagram_grid, var=1::Integer, args...; kwargs...)
    for i = 1:bif_dia_grid.length
        foreach(bp -> plot!(fill(bif_dia_grid.range1[i],bp.length),bp.p_vals,getindex.(bp.vals,var), args...; label="", color=stab_color.(bp.stability_types), kwargs...), bif_dia_grid.bifurcation_diagrams[i].paths)
    end
    plot!(args...; xlabel=String(bif_dia_grid.param1), zlabel=String(bif_dia_grid.param2), kwargs...)
end
#Makes a scatter type of plot versio for the bifurcation grid diagrams.
function bif_scatter(bif_dia_grid::bifurcation_diagram_grid, var=1::Integer, args...; kwargs...)
    scatter(); bif_scatter!(bif_dia_grid, var, args...; kwargs...)
end
function bif_scatter!(bif_dia_grid::bifurcation_diagram_grid, var=1::Integer, args...; kwargs...)
    for i = 1:bif_dia_grid.length
        for bp in bif_dia_grid.bifurcation_diagrams[i].paths
            for j = 1:bp.length
                scatter!((bif_dia_grid.range1[i],bp.p_vals[j],bp.vals[j][var]), args...; label="", color=stab_color(bp.stability_types[j]), kwargs...)
            end
        end
    end
    scatter!(args...; xlabel=bif_dia_grid.param1, zlabel=bif_dia_grid.param2, kwargs...)
end
