import Base: length
length(p::bifurcation_path) = length(p.p_vals)
length(p::bifurcation_diagram_grid) = p.length


@recipe function f(path::bifurcation_path, var=1)
    linewidth --> 4
    color --> stab_color.(path.stability_types)
    path.p_vals, getindex.(path.vals, var)
end

@recipe function f(path::bifurcation_path, var, z_val)
    linewidth --> 4
    color --> stab_color.(path.stability_types)
    fill(z_val, length(path)), path.p_vals, getindex.(path.vals, var)
end

@recipe function f(bif::bifurcation_diagram, var, z_val)
    xlabel --> string(bif.param)
    for (i, path) in enumerate(bif.paths)
        @series begin
            primary --> i==1
            path, var, z_val
        end
    end
end

@recipe function f(bif::bifurcation_diagram, var=1)
    xlabel --> string(bif.param)
    for (i, path) in enumerate(bif.paths)
        @series begin
            primary --> i==1
            path, var
        end
    end
end

@recipe function f(grid::bifurcation_diagram_grid, var=1)
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

@recipe function f(x, point::bifurcation_point)
    seriestype := :scatter
    color --> stab_color.(point.stability_types)
    fill(x, length(point.vals)), point.vals
end

@recipe function f(x, y, point::bifurcation_point)
    seriestype := :scatter
    color --> stab_color.(point.stability_types)
    fill(x, length(point.vals)), fill(y, length(point.vals)), point.vals
end

@recipe function f(grid::bifurcation_grid)
    label --> ""
    ylabel --> "Concentration"
    xlabel --> string(grid.param)
    for (i, point) in enumerate(grid.grid_points)
        @series begin
            primary --> i==1
            grid.range[i], point
        end
    end
end

@recipe function f(grid::bifurcation_grid_2d)
    label --> ""
    zlabel --> "Concentration"
    ylabel --> string(grid.param1)
    xlabel --> string(grid.param2)
    for i = 1:length(grid.range1), j = 1:length(grid.range2)
        @series begin
            primary --> i==1
            grid.range1[i], grid.range2[j], grid.grid_points[i,j]
        end
    end
end
