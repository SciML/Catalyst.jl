GillespieProblem(prob,aggregator::AbstractAggregatorAlgorithm,
                 rs::AbstractReaction...;kwargs...) =
                 JumpProblem(prob,aggregator,build_jumps_from_reaction(rs...);kwargs...)

GillespieProblem(prob,aggregator::AbstractAggregatorAlgorithm,
                 r::ReactionSet;kwargs...) =
                 JumpProblem(prob,aggregator,build_jumps_from_reaction(r);kwargs...)

### Added by Me ###
### ODEProblem ###
newODEProblem(rn::ReactionNetwork, args...; kwargs...) =
                 ODEProblem(rn.f, args...; kwargs...)

### SDEProblem ###
newSDEProblem(rn::ReactionNetwork, args...; kwargs...) =
                 SDEProblem(rn.f,rn.g, args...;noise_rate_prototype=rn.p_matrix, kwargs...)

### JumpProblem ###
newJumpProblem(prob,aggregator::Direct,rn::ReactionNetwork) =
                 JumpProblem(prob,aggregator::Direct,rn.jumps...)
