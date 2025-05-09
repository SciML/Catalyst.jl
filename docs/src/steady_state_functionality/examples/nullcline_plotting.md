# [Plotting Nullclines and Steady States in Phase Space](@id nullcline_plotting)
In this tutorial we will show how to extract a system's steady states and [nullclines](https://en.wikipedia.org/wiki/Nullcline), and how to plot these in [phase space](https://en.wikipedia.org/wiki/Phase_space). Generally, while nullclines are not directly needed for much analysis, plotting these can give some understanding of a system's steady state and stability properties.

For an ordinary differential equation
```math
\begin{aligned}
\frac{dx_1}{dt} &= f_1(x_1, x_2, ..., x_n) \\
\frac{dx_2}{dt} &= f_2(x_1, x_2, ..., x_n) \\
                &\vdots\\
\frac{dx_n}{dt} &= f_n(x_1, x_2, ..., x_n) \\
\end{aligned}
```
the $i$'th nullcline is the surface along which $\frac{dx_i}{dt} = 0$, i.e. the implicit surface given by $f_i(x_1,\dots,x_n) = 0$. Nullclines are frequently used when visualizing the phase-planes of two-dimensional models (as these can be easily plotted).

## [Computing nullclines and steady states for a bistable switch](@id nullcline_plotting_computation)
For our example we will use a simple bistable switch model, consisting of two species ($X$ and $Y$) which mutually inhibit each other through repressive Hill functions. 
```@example nullcline_plotting
using Catalyst
bs_switch = @reaction_network begin
    hillr(Y, v, K, n), 0 --> X
    hillr(X, v, K, n), 0 --> Y
    1.0, (X,Y) --> 0
end
```

Next, we compute the steady states [using homotopy continuation](@ref homotopy_continuation).
```@example nullcline_plotting
import HomotopyContinuation
ps = [v => 1.0, K => 0.6, n => 4.0]
sss = hc_steady_states(bs_switch, ps; show_progress = false)
```

Finally, we will compute the nullclines. First we create a function which, for species values $(X,Y)$, returns the evaluation of the model's ODE's right-hand side.
```@example nullcline_plotting
nlprob = NonlinearProblem(complete(nlsys), [X => 0.0, Y => 0.0], ps)
function get_XY(Xval, Yval)
    prob = Catalyst.remake(nlprob; u0 = [X => Xval, Y => Yval])
    return prob[equations(nlsys)[2].rhs], prob[equations(nlsys)[1].rhs]
end
```
Next, we plot our nullclines by using Plot.jl's [`contour` function](https://docs.juliaplots.org/latest/series_types/contour/). Here, we plot the $0$ contour lines (which corresponds to the nullclines) for both the $X$ and $Y$ nullclines. We will plot the steady states using `scatter`. We use the `steady_state_stability` function to [compute steady state stabilities](@ref steady_state_stability) (and use this to determine how to plot the steady state markers).
```@example nullcline_plotting
using Plots
# Plot the nullclines. Line labels added in separate `plot` commands (due to how the `contour` works).
plt_mesh = 0:0.02:1.25
contour(plt_mesh, plt_mesh, (x,y) -> get_XY(x,y)[1]; levels = [0.0], lw = 7, la = 0.7, color = 1, cbar=false)
contour!(plt_mesh, plt_mesh, (x,y) -> get_XY(x,y)[2]; levels = [0.0], lw = 7, la = 0.7, color = 2, cbar=false)
plot!([]; label = "dX/dt = 0", lw = 7, la = 0.7, color = 1)
plot!([]; label = "dY/dt = 0", lw = 7, la = 0.7, color = 2)

# Plot the steady states.
for ss in sss
    color, markershape = if steady_state_stability(ss, bs_switch, ps)
        :blue, :circle
    else
        :red, :star4
    end
    scatter!((ss[2], ss[1]); color, markershape, label = "", markersize = 10)
end
scatter!([], []; color = :blue, markershape = :circle, label = "Stable stead state")
scatter!([], []; color = :red, markershape = :star4, label = "Unstable stead state")

# Finishing touches.
plot!(xlimit = span, ylimit = span, xguide = "X", yguide = "Y", legendfontsize = 10, size = (600,600))
```
Here we can see how the steady states occur at the nullclines intersections.

!!! note
    Here we use an inherent Plots function to plot the nullclines. However, there are also specialised packages for these kinds of plots, such as [ImplicitPlots.jl](https://github.com/saschatimme/ImplicitPlots.jl).

## [Plotting system directions in phase space](@id nullcline_plotting_directions)
One useful property of nullclines is that the sign of $dX/dt$ will only switch whenever the solution crosses the $dX/dt=0$ nullcline. This means that, within each region defined by the nullclines, the direction of the solution remains constant. Below we use this to, for each such region, plot arrows showing the solution's direction.
```@example nullcline_plotting
# Creates a function for plotting the ODE's direction at a point in phase space.
function plot_xy_arrow!(Xval, Yval)
    dX, dY = get_XY(Xval, Yval)
    dX = dX > 0 ? 0.05 : -0.05
    dY = dY > 0 ? 0.05 : -0.05
    plot!([Xval, Xval + dX], [Yval, Yval]; color = :black, lw = 2, arrow = :arrow, label = "")
    plot!([Xval, Xval], [Yval, Yval + dY]; color = :black, lw = 2, arrow = :arrow, label = "")
end

# Plots the ODE's direction in phase space at selected positions.
arrow_positions = [(0.25, 0.25), (0.75, 0.75), (0.35, 0.8), (0.8, 0.35), (0.02, 1.1), (1.1, 0.02)]
foreach(pos -> plot_xy_arrow!(pos...), arrow_positions)
plot!()
```
This also works as a form of simple stability analysis, where we can see how the solution moves *away* from the unstable steady state, and *to* the stable ones.
