using DiffEqBiological, DiffEqJump, Test

# only should agree for pure MassAction or pure ConstantRateJumps
# mixed systems will be reordered when generating jump problems!

# pure mass action
dnanetwork = @reaction_network begin
    c1, G --> G + M
    c2, M --> M + P
    c3, M --> 0
    c4, P --> 0
    c5, 2P --> P2
    c6, P2 --> 2P
    c7, P2 + G --> P2G
    c8, P2G --> P2 + G
end c1 c2 c3 c4 c5 c6 c7 c8
rnpar = [.09, .05, .001, .0009, .00001, .0005, .005, .9]
varlabels = ["G", "M", "P", "P2","P2G"]
u0 = [1000, 0, 0, 0,0]
tf = 4000.
prob = DiscreteProblem(dnanetwork, u0, (0.0, tf), rnpar)
jprob = JumpProblem(prob, RSSA(), dnanetwork)
@test all(rxtospecies_depgraph(dnanetwork) .== jprob.discrete_jump_aggregation.jumptovars_map)
@test all(speciestorx_depgraph(dnanetwork) .== jprob.discrete_jump_aggregation.vartojumps_map)

jprob = JumpProblem(prob, NRM(), dnanetwork)
@test all(DiffEqBiological.rxtorx_depgraph(dnanetwork) .== jprob.discrete_jump_aggregation.dep_gr)

# pure ConstantRateJump
hillnetwork = @reaction_network begin
    hillr(m₃,α,K,n), ∅ --> m₁
    hillr(m₁,α,K,n), ∅ --> m₂
    hill(m₂,α,K,n), ∅ --> m₃
end α K η
p = (1.,1.,1.)
u0 = [10,10,10]
tf = 10.
prob = DiscreteProblem(hillnetwork, u0, (0.0, tf), rnpar)
jprob = JumpProblem(prob, RSSA(), hillnetwork)
@test all(rxtospecies_depgraph(hillnetwork) .== jprob.discrete_jump_aggregation.jumptovars_map)
@test all(speciestorx_depgraph(hillnetwork) .== jprob.discrete_jump_aggregation.vartojumps_map)

jprob = JumpProblem(prob, NRM(), hillnetwork)
@test all(DiffEqBiological.rxtorx_depgraph(hillnetwork) .== jprob.discrete_jump_aggregation.dep_gr)
