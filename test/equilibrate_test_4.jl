#Tests various bifurcation capabilities.
for rn in reaction_networks_standard
    p_vals = map(amp->amp*rand(length(rn.params)),[1.,10.,100.])
    param1 = rand(rn.params)
    param2 = rand(rn.params)
    range_tupple = (1.,rand(5.:10.))
    range1 = 1.: rand(5.:10.)
    range2 = 1.: rand(5.:10.)
    for p in p_vals
        bif1 = bifurcations(rn,p,param1,range_tupple)
        bif2 = bifurcations(rn,p,param1,range_tupple,solver=HcBifurcationSolverSimple)
        bif_grid = bifurcations_grid(rn,p,param1,range1)
        bif_grid_2d = bifurcations_grid_2d(rn,p,param1,range1,param2,range2)
        bif_grid_dia = bifurcations_diagram_grid(rn,p,param1,range1,param2,range_tupple)
        plot(bif1); plot(bif2); plot(bif_grid); plot(bif_grid_2d); plot(bif_grid_dia);
    end
end
