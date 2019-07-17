#Tests on the networks containing hill functions.
for rn in reaction_networks_hill
    p_vals = map(amp->rand(1:amp,length(rn.params)),[2,3,6])
    for p in p_vals
        fix_parameters(rn,p)
        fps = steady_states(rn,Vector{Float64}())
        stab = map(fp->stability(fp,rn,Float64.(p)), fps)
        for i=1:length(fps)
            stab[i] ? (tend = 10) : (tend = -10)
            end_point = OrdinaryDiffEq.solve(ODEProblem(rn,fps[i],(0.,tend),p),Rosenbrock23())[end]
            @test maximum(abs.(end_point-fps[i]))<0.0001
        end
        #start_points = map(amp->amp*rand(length(rn.syms)),[1.,10.,100.])
        #for sp in start_points
        #    end_point = OrdinaryDiffEq.solve(ODEProblem(rn,sp,(0.,1000.),p),Rosenbrock23())[end]
        #    @test any(map(fp->maximum(abs.(fp-end_point)),fps) .< 0.001)
        #end
    end
end
