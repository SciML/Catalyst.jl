GillespieProblem(prob,aggregator::AbstractAggregatorAlgorithm,
                 rs::AbstractReaction...;kwargs...) =
                 JumpProblem(prob,aggregator,build_jumps_from_reaction(rs...);kwargs...)
                 
GillespieProblem(prob,aggregator::AbstractAggregatorAlgorithm,
                 r::ReactionSet;kwargs...) =
                 JumpProblem(prob,aggregator,build_jumps_from_reaction(r);kwargs...)
