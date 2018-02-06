GillespieProblem(prob,aggregator::AbstractAggregatorAlgorithm,
                 rs::AbstractReaction...;kwargs...) =
                 JumpProblem(prob,aggregator,build_jumps_from_reaction(rs...);kwargs...)

GillespieProblem(prob,aggregator::AbstractAggregatorAlgorithm,
                 r::ReactionSet;kwargs...) =
                 JumpProblem(prob,aggregator,build_jumps_from_reaction(r);kwargs...)

### Added by Me ###
### ODEProblem ###
newODEProblem(rn::ReactionNetwork,u0,tspan) =
                 ODEProblem(rn.f,u0,tspan)

### SDEProblem ###
newSDEProblem(rn::ReactionNetwork,u0,tspan) =
                 SDEProblem(rn.f,rn.g,u0,tspan,noise_rate_prototype=rn.p_matrix)

### JumpProblem ###
newJumpProblem(prob,aggregator::Direct,rn::ReactionNetwork) =
                 JumpProblem(prob,aggregator::Direct,rn.jumps...)
