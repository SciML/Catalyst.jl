#Tests steady state finders for minimal reaction network.
for rn in min_reaction_networks_standard
    addequi!(rn)
    manage_reaction_network!(rn)
    p_vals = map(amp->amp*rand(length(rn.params)),[1.,10.,100.])
    for p in p_vals
        fps = steady_states(rn,p)
        stab = map(fp->stability(fp,p,rn), fps)
        for i=1:length(fps)
            stab[i] ? (tend = 10) : (tend = -10)
            end_point = DifferentialEquations.solve(ODEProblem(rn,fps[i],(0.,tend),p),Rosenbrock23())[end]
            @test maximum(end_point-fps[i])<0.0001
        end
    end
end
