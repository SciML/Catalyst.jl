import Base: length
length(p::BifurcationPath) = length(p.p_vals)
length(p::BifurcationGridDiagram) = p.length


@recipe function f(path::BifurcationPath, var=1)
    linewidth --> 4
    color --> stab_color.(path.stability_types)
    path.p_vals, getindex.(path.vals, var)
end

@recipe function f(path::BifurcationPath, var, z_val)
    linewidth --> 4
    color --> stab_color.(path.stability_types)
    fill(z_val, length(path)), path.p_vals, getindex.(path.vals, var)
end

@recipe function f(bif::BifurcationDiagram, var, z_val)
    xlabel --> string(bif.param)
        label --> ""
    for (i, path) in enumerate(bif.paths)
        @series begin
            primary --> i==1
            path, var, z_val
        end
    end
end

@recipe function f(bif::BifurcationDiagram, var=1)
    xlabel --> string(bif.param)
        label --> ""
    for (i, path) in enumerate(bif.paths)
        @series begin
            primary --> i==1
            path, var
        end
    end
end

@recipe function f(grid::BifurcationGridDiagram, var=1)
    label --> ""
    zlabel --> "Concentration"
    ylabel --> string(grid.param2)
    xlabel --> string(grid.param1)
    for (i, bif) in enumerate(grid.bifurcation_diagrams)
        @series begin
            bif, var, grid.range1[i]
        end
    end
end

@recipe function f(u, point::BifurcationPoint, var)
    seriestype := :scatter
    color --> stab_color.(point.stability_types)
    fill(u, length(point.vals)), getindex.(point.vals,var)
end

@recipe function f(u, y, point::BifurcationPoint, var)
    seriestype := :scatter
    color --> stab_color.(point.stability_types)
    fill(u, length(point.vals)), fill(y, length(point.vals)), getindex.(point.vals,var)
end

@recipe function f(grid::BifurcationGrid, var=1)
    label --> ""
    ylabel --> "Concentration"
    xlabel --> string(grid.param)
    for (i, point) in enumerate(grid.grid_points)
        @series begin
            primary --> i==1
            grid.range[i], point, var
        end
    end
end

@recipe function f(grid::BifurcationGrid2D, var=1)
    label --> ""
    zlabel --> "Concentration"
    ylabel --> string(grid.param2)
    xlabel --> string(grid.param1)
    for i = 1:length(grid.range1), j = 1:length(grid.range2)
        @series begin
            primary --> i==1
            grid.range1[i], grid.range2[j], grid.grid_points[i,j], var
        end
    end
end

#For a specific stability types, returns a color (to be used to distinguish the types in plots).
function stab_color(stab_type::Int8)
    (stab_type==0) && (return :red)
    (stab_type==1) && (return :blue)
    (stab_type==2) && (return :orange)
    (stab_type==3) && (return :cyan)
end
