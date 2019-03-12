### ODEProblem from AbstractReactionNetwork ###
function DiffEqBase.ODEProblem(rn::DiffEqBase.AbstractReactionNetwork, u0::Union{AbstractArray, Number}, args...; kwargs...) 
    isa(rn.odefun, ODEFunction) || error("Call addodes! before constructing ODEProblems")
    ODEProblem(rn.odefun, u0::Union{AbstractArray, Number}, args...; kwargs...)
end

### SDEProblem from AbstractReactionNetwork ###
function DiffEqBase.SDEProblem(rn::DiffEqBase.AbstractReactionNetwork, u0::Union{AbstractArray, Number}, args...; kwargs...) 
    isa(rn.sdefun, SDEFunction) || error("Call addsdes! before constructing SDEProblems")
    SDEProblem(rn.sdefun, rn.g::Function, u0, args...;noise_rate_prototype=rn.p_matrix, kwargs...)
end 

### DiscreteProblem, passing through syms
function DiffEqBase.DiscreteProblem(rn::DiffEqBase.AbstractReactionNetwork, u0, tspan::Tuple, p=nothing; kwargs...)
    f = (du,u,p,t) -> du.=u    # identity function to make syms works
    df = DiscreteFunction(f,syms=rn.syms)
    DiscreteProblem(df, u0, tspan, p; kwargs...)
end

### JumpProblem ###
function build_jump_problem(prob, aggregator, rn, kwargs...)
    if typeof(prob)<:DiscreteProblem && any(x->typeof(x) <: VariableRateJump, rn.jumps)
        error("When using time dependant reaction rates a DiscreteProblem should not be used (try an ODEProblem). Also, use a continious solver.")
    end

    # get a JumpSet of the possible jumps
    jset = network_to_jumpset(rn, prob.p)

    # construct map from species index to indices of reactions that depend on it
    if needs_vartojumps_map(aggregator) || needs_depgraph(aggregator)
        rxidxs_to_jidxs   = rxidxs_to_jidxs_map(rn, get_num_majumps(jset))
        spec_to_jumps_set = spec_to_dep_jumps_map(rn, rxidxs_to_jidxs)
        spec_to_jumps_vec = [sort!(collect(specset)) for specset in spec_to_jumps_set]
        jump_to_specs_vec = jump_to_dep_specs_map(rn, rxidxs_to_jidxs)
    else
        rxidxs_to_jidxs   = nothing
        spec_to_jumps_set = nothing
        spec_to_jumps_vec = nothing
        jump_to_specs_vec = nothing
    end

    # construct reaction dependency graph
    if needs_depgraph(aggregator)
        dep_graph = depgraph_from_network(rn, jset, rxidxs_to_jidxs, spec_to_jumps_set)
    else
        dep_graph = nothing
    end

    JumpProblem(prob, aggregator, jset; dep_graph=dep_graph,
                                        vartojumps_map=spec_to_jumps_vec,
                                        jumptovars_map=jump_to_specs_vec,
                                        kwargs...)
end

### JumpProblem from AbstractReactionNetwork ###
function DiffEqJump.JumpProblem(prob, aggregator, rn::DiffEqBase.AbstractReactionNetwork; kwargs...)
    (rn.jumps !== nothing) || error("Call addjumps! before constructing JumpProblems")
    build_jump_problem(prob, aggregator, rn, kwargs...)
end

### SteadyStateProblems from AbstractReactionNetwork ###
function DiffEqBase.SteadyStateProblem(rn::DiffEqBase.AbstractReactionNetwork, args...; kwargs...) 
    isa(rn.odefun, ODEFunction) || error("Call addodes! before constructing SteadyStateProblems")
    SteadyStateProblem(rn.odefun, args...; kwargs...)
end

function DiffEqBase.SteadyStateProblem{isinplace}(rn::DiffEqBase.AbstractReactionNetwork, args...; kwargs...) where isinplace
    isa(rn.odefun, ODEFunction) || error("Call addodes! before constructing SteadyStateProblems")
    SteadyStateProblem{isinplace}(rn.odefun, args...; kwargs...)
end
