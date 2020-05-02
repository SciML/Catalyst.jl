### ODEProblem from AbstractReactionNetwork ###
function MT_ODEProblem(rn::DiffEqBase.AbstractReactionNetwork, u0::Union{AbstractArray, Number}, tspan, p, args...; kwargs...)
    (rn.ode_system == nothing) && (rn.ode_system = convert(ODESystem,rn.reaction_system))
    return ODEProblem(rn.ode_system,Pair.(rn.reaction_system.states,u0),tspan,Pair.(rn.reaction_system.ps,p))
end
