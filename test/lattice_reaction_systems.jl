### Very simple tests to have during development ###

using Catalyst, OrdinaryDiffEq, Graphs

rs = @reaction_network begin
    A, ∅ → X
    1, 2X + Y → 3X
    B, X → Y
    1, X → ∅
end A B
srs = [SpatialReaction(:D, ([:X], []), ([], [:X]), ([1], []), ([], [1]))]
lattice = grid([20, 20])
lrs = LatticeReactionSystem(rs, srs, lattice);

u0_in = [:X => 10 * rand(nv(lattice)), :Y => 10 * rand(nv(lattice))]
tspan = (0.0, 100.0)
p_in = [:A => 1.0, :B => 4.0, :D => 0.2]

oprob = ODEProblem(lrs, u0_in, tspan, p_in)
@test solve(oprob, Tsit5()).retcode == :Success