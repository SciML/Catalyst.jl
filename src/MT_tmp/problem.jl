### ODEProblem from AbstractReactionNetwork ###
function MT_ODEProblem(rn::DiffEqBase.AbstractReactionNetwork, u0::Union{AbstractArray, Number}, tspan, p, args...; kwargs...)
    (rn.ode_system == nothing) && (rn.ode_system = convert(ODESystem,rn.reaction_system))
    return ODEProblem(rn.ode_system,Pair.(rn.reaction_system.states,u0),tspan,Pair.(rn.reaction_system.ps,p), args...; kwargs...)
end

### SDEProblem from AbstractReactionNetwork ###
function MT_SDEProblem(rn::DiffEqBase.AbstractReactionNetwork, u0::Union{AbstractArray, Number}, tspan, p, args...; kwargs...)
    (rn.sde_system == nothing) && (rn.sde_system = convert(SDESystem,rn.reaction_system))
    p_matrix = zeros(length(rn.reaction_system.states), length(rn.reaction_system.eqs))
    return SDEProblem(rn.sde_system,Pair.(rn.reaction_system.states,u0),tspan,Pair.(rn.reaction_system.ps,p),args...; noise_rate_prototype=p_matrix,kwargs...)
end
