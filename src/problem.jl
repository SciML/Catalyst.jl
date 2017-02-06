GillespieProblem(prob,aggregator::AbstractAggregatorAlgorithm,rs::Reaction...;kwargs...) =
                      JumpProblem(prob,aggregator,build_jumps_from_reaction(rs...);kwargs...)
