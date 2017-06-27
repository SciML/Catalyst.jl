using DiffEqBiological

# Birth-Death-Immigration Process
bdip = @reaction_network begin
  2.0, X --> X + X
  1.0, X --> 0
  0.5, 0 --> X
end

bdip1 = Reaction(2.0, [1], [(1, 1)])
bdip2 = Reaction(1.0, [1], [(1, -1)])
bdip3 = Reaction(0.5, [],  [(1, 1)])

@test bdip.reactions[1] == bdip1
@test bdip.reactions[2] == bdip2
@test bdip.reactions[3] == bdip3
