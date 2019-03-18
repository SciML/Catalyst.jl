#Tests various bifurcation capabilities.
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
