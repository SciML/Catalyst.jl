# Fetch packages.
using Catalyst, JumpProcesses, OrdinaryDiffEqTsit5, Statistics, Test

# Sets stable rng number.
using StableRNGs
rng = StableRNG(12345)
seed = rand(rng, 1:100)

# test JumpInputs function auto problem selection
let
    rn = @reaction_network begin
        k*(1 + sin(t)), 0 --> A
    end
    jinput = JumpInputs(rn, [:A => 0], (0.0, 10.0), [:k => .5])
    @test jinput.prob isa ODEProblem
    jprob = JumpProblem(jinput; rng)
    sol = solve(jprob, Tsit5())
    @test sol(10.0; idxs = :A) > 0

    rn = @reaction_network begin
        k, 0 --> A
    end
    jinput = JumpInputs(rn, [:A => 0], (0.0, 10.0), [:k => .5])
    @test jinput.prob isa DiscreteProblem
    jprob = JumpProblem(jinput; rng)
    sol = solve(jprob)
    @test sol(10.0; idxs = :A) > 0

    rn = @reaction_network begin
        @parameters λ
        k*V, 0 --> A
        @equations D(V) ~ λ*V
        @continuous_events begin
            [V ~ 2.0] => [V ~ V/2, A ~ A/2]
        end
    end        
    jinput = JumpInputs(rn, [:A => 0, :V => 1.0], (0.0, 10.0), [:k => 1.0, :λ => .4])
    @test jinput.prob isa ODEProblem
    jprob = JumpProblem(jinput; rng)
    sol = solve(jprob, Tsit5()) 
end

# this test requires a decision on how to handle when a user gives 
# an ODE for a chemical species in addition to reactions 
# let
#     seed = 1111
#     Random.seed!(rng, seed)
#     rn = @reaction_network begin
#         @parameters α β
#         β, X --> 0
#         β, Y --> 0  
#         α, 0 --> Y
#         @equations begin
#             D(X) ~ α * (1 + Y)
#         end
#     end    
#     p = (α = 6.0, β = 2.0, X₀ = 2.0, Y₀ = 1.0)
#     u0map = [:X => p.X₀, :Y => p.Y₀]
#     pmap = [:α => p.α, :β => p.β]
#     tspan = (0.0, 20.0)
#     jinputs = JumpInputs(rn, u0map, tspan, pmap)
#     jprob = JumpProblem(jinputs; rng, save_positions = (false, false))
#     times = range(0.0, tspan[2], length = 100)
#     Nsims = 4000
#     Xv = zeros(length(times))
#     Yv = zeros(length(times))
#     for n in 1:Nsims
#         sol = solve(jprob, Tsit5(); saveat = times, seed)
#         Xv .+= sol[1, :]
#         Yv .+= sol[2, :]
#         seed += 1
#     end
#     Xv ./= Nsims
#     Yv ./= Nsims

#     function Yf(t, p)
#         local α, β, X₀, Y₀ = p
#         return (α / β) + (Y₀ - α / β) * exp(-β * t)
#     end
#     function Xf(t, p)
#         local α, β, X₀, Y₀ = p
#         return (α / β) + (α^2 / β^2) + α * (Y₀ - α / β) * t * exp(-β * t) +
#                (X₀ - α / β - α^2 / β^2) * exp(-β * t)
#     end
#     Xact = [Xf(t, p) for t in times]
#     Yact = [Yf(t, p) for t in times]
#     @test all(abs.(Xv .- Xact) .<= 0.05 .* Xv)
#     @test all(abs.(Yv .- Yact) .<= 0.05 .* Yv)

#     function affect!(integ, u, p, ctx)
#         savevalues!(integ, true)
#         terminate!(integ)
#         nothing
#     end
#     cevents = [t ~ 0.2] => (affect!, [], [], [], nothing)
#     @named jsys = JumpSystem([maj, crj, vrj, eqs[1]], t, [X, Y], [α, β];
#         continuous_events = cevents)
#     jsys = complete(jsys)
#     tspan = (0.0, 200.0)
#     oprob = ODEProblem(jsys, u0map, tspan, pmap)
#     jprob = JumpProblem(jsys, oprob; rng, save_positions = (false, false))
#     Xsamp = 0.0
#     Nsims = 4000
#     for n in 1:Nsims
#         sol = solve(jprob, Tsit5(); saveat = tspan[2], seed)
#         @test sol.retcode == ReturnCode.Terminated
#         Xsamp += sol[1, end]
#         seed += 1
#     end
#     Xsamp /= Nsims
#     @test abs(Xsamp - Xf(0.2, p) < 0.05 * Xf(0.2, p))
# end