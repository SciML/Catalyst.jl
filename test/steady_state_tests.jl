#Tests on the standard networks.
for rn in reaction_networks_standard
    p_vals = map(amp->amp*rand(length(rn.params)),[1.,10.,100.])
    for p in p_vals
        fps = steady_states(rn,p)
        stab = map(fp->stability(fp,p,rn), fps)
        for i=1:length(fps)
            stab[i] ? (tend = 10) : (tend = -10)
            end_point = OrdinaryDiffEq.solve(ODEProblem(rn,fps[i],(0.,tend),p),Rosenbrock23())[end]
            @test maximum(end_point-fps[i])<0.0001
        end
    end
end

#Tests on the networks containing hill functions.
for rn in reaction_networks_hill
    p_vals = map(amp->rand(1:amp,length(rn.params)),[2,3,6])
    for p in p_vals
        fix_parameters(rn,p)
        fps = steady_states(rn,Vector{Float64}())
        stab = map(fp->stability(fp,Float64.(p),rn), fps)
        for i=1:length(fps)
            stab[i] ? (tend = 10) : (tend = -10)
            end_point = OrdinaryDiffEq.solve(ODEProblem(rn,fps[i],(0.,tend),p),Rosenbrock23())[end]
            @test maximum(end_point-fps[i])<0.0001
        end
    end
end

#Tests on the networks requiring fixed concentrations.
for rn in reaction_networks_fixed_conc
    p_vals = map(amp->amp*rand(length(rn.params)),[1.,10.,100.])
    for p in p_vals
        fps = steady_states(rn,p)
        stab = map(fp->stability(fp,p,rn), fps)
        for i=1:length(fps)
            stab[i] ? (tend = 10) : (tend = -10)
            end_point = OrdinaryDiffEq.solve(ODEProblem(rn,fps[i],(0.,tend),p),Rosenbrock23())[end]
            @test maximum(end_point-fps[i])<0.0001
        end
    end
end


#Tests various bifurcation capabilities
for rn in reaction_networks_standard
    p_vals = map(amp->amp*rand(length(rn.params)),[1.,10.,100.])
    param1 = rand(rn.params)
    param2 = rand(rn.params)
    range_tupple = (1.,rand(5.:10.))
    range1 = 1.: rand(5.:10.)
    range2 = 1.: rand(5.:10.)
    for p in p_vals
        bif_plot(bifurcations(rn,p,param1,range_tupple))
        bif_plot(bifurcations_grid(rn,p,param1,range1))
        bif_plot(bifurcations_grid_2d(rn,p,param1,range1,param2,range2))
        bdg = bifurcations_diagram_grid(rn,p,param1,range1,param2,range_tupple)
        bif_plot!(bdg); bif_plot(bdg); bif_scatter(bdg); bif_scatter!(bdg);
    end
end
