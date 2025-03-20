# [Plotting Nullclines and Steady States in Phase Space](@id nullcline_plotting)
In this tutorial we will show how to extract a system's steady states and [nullclines](https://en.wikipedia.org/wiki/Nullcline), and how to plot these in [phase space](https://en.wikipedia.org/wiki/Phase_space). Generally, while nullclines are not directly needed for much analysis, plotting these can give some understanding of a systems steady state and stability properties. While nullclines can be "brute forced", we will here use [BifurcationKit.jl](https://github.com/bifurcationkit/BifurcationKit.jl) to rack their path exactly.

For a ordinary differential equation
```math
\begin{aligned}
\frac{dx_1}{dt} &= f_1(x_1, x_2, ..., x_3) \\
\frac{dx_2}{dt} &= f_2(x_1, x_2, ..., x_3) \\
                &\vdots\\
\frac{dx_n}{dt} &= f_3(x_1, x_2, ..., x_3) \\
\end{aligned}
```
The nullclines are the curves
```math
\begin{aligned}
0 &= f_1(x_1, x_2, ..., x_3) \\
0 &= f_2(x_1, x_2, ..., x_3) \\
                &\vdots\\
0 &= f_3(x_1, x_2, ..., x_3) \\
\end{aligned}
```
where the $i$'th nullclines is the curve along which $\frac{dx_i}{dt} = 0$. Generally, nullclines are primarily computed for models with 2 variable (as only here can they be easily plotted).

## [Computing nucclines and steady states for a bistable switch](@id nullcline_plotting_computation)
For our example we will use a simple bistable switch model, consisting of two species ($X$ and $Y$) which mutually inhibits each other through repressive hill functions. We will create our model [programmatically](@ref programmatic_CRN_construction).
```@example nullcline_plotting
using Catalyst
t = default_t()
@parameters v K n
@species X(t) Y(t)
rxs = [
    Reaction(hillr(Y, v, K, n), [], [X]),
    Reaction(hillr(X, v, K, n), [], [Y]),
    Reaction(1.0, [X], []),
    Reaction(1.0, [Y], []),
]
@named bs_switch = ReactionSystem(rxs, t)
bs_switch = complete(bs_switch)
```

Next, we compute the steady states [using homotopy continuation](@ref homotopy_continuation).
```@example nullcline_plotting
import HomotopyContinuation
ps = [v => 1.0, K => 0.6, n => 4.0]
sss = hc_steady_states(bs_switch, ps; show_progress = false)
```

Finally, we will compute the nullclines. We will first extract our model's steady state equations (which we do by creating a `NonlinearSystem`).
```@example nullcline_plotting
nlsys = convert(NonlinearSystem, bs_switch)
eqs = equations(nlsys)
```
Here, each equation is an expression of two variables ($X$ and $Y$). We wish to plot $Y$ as a function of $X$. To do this, we will substitute the species variable $X$ with the parameter $Xpar$. This will enable us to carry out bifurcation analysis of the equations' solutions as $Xpar$ is varied
```@example nullcline_plotting
@parameters Xpar
nc_eq_X = substitute(equations(nlsys)[1], Dict([X => Xpar]))
nc_eq_Y = substitute(equations(nlsys)[2], Dict([X => Xpar]))
```
To input these into BifurcationKit we need to convert these into `NonlinearSystem`s.
```@example nullcline_plotting
@named nc_X_sys = NonlinearSystem([nc_eq_X])
@named nc_Y_sys = NonlinearSystem([nc_eq_Y])
nc_X_sys = complete(nc_X_sys)
nc_Y_sys = complete(nc_Y_sys)
nothing # hide
```
Finally, for out nullcline equations, we can use BifurcationKit's `continuation` function to track their solutions across all values of $Xpar$ (the workflow is similar to when we used `bifurcationdiagram` [here](@ref bifurcation_diagrams)).
```@example nullcline_plotting
using BifurcationKit
span = (0.0, 1.2)
function compute_nullcline(nc_sys)
    bprob = BifurcationProblem(nc_sys, [Y => 1.0], [ps; Xpar => 0.1], Xpar)
    opts_br = ContinuationPar(p_min = span[1], p_max = span[2], dsmax = 0.01)
    return continuation(bprob, PALC(), opts_br; bothside = true)
end
nc_X = compute_nullcline(nc_X_sys)
nc_Y = compute_nullcline(nc_Y_sys)
nothing # hide
```

We are ready to create our plot. We will plot the steady states using `scatter`. We use the `steady_state_stability` function to [compute steady state stabilities](@ref steady_state_stability) (and use this to determine how to plot the steady state markers).
```@example nullcline_plotting
using Plots
plot(nc_X.x, nc_X.param; label = "dX/dt = 0", lw = 7, la = 0.7)
plot!(nc_Y.x, nc_Y.param; label = "dY/dt = 0", lw = 7, la = 0.7)
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
plot!(xlimit = span, ylimit = span, xguide = "X", yguide = "Y", legendfontsize = 10, size = (500,500))
```
Here we can see how the steady states occur at the nullclines intersections.

!!! warn
    BifurcationKit's `continuation` function will not detect disjoint branches, and above we can only use it because we know that each nullcline consists of a single branch. To handle disjoint nullclines, consider using [deflated continuation](@ref bifurcation_diagrams_disjoint_branches) (or possibly a package like [Contour.jl](https://github.com/JuliaGeometry/Contour.jl)).

## [Plotting system directions in phase space](@id nullcline_plotting_directions)
One useful property of nullclines is that the sign of $dX/dt$ will only switch whenever the solution crosses the $dX/dt=0$ nullcline. This mean that, within each region defined by the nullclines, the direction of the solution remains constant. Below we use this to, for each such region, plot arrows showing the solution's direction.

```@example nullcline_plotting
nlprob = NonlinearProblem(complete(nlsys), [X => 0.0, Y => 0.0], ps)
function get_XY(Xval, Yval)
    prob = Catalyst.remake(nlprob; u0 = [X => Xval, Y => Yval])
    return prob[equations(nlsys)[2].rhs], prob[equations(nlsys)[1].rhs]
end
function plot_xy_arrow!(Xval, Yval)
    dX, dY = get_XY(Xval, Yval)
    dX = dX > 0 ? 0.05 : -0.05
    dY = dY > 0 ? 0.05 : -0.05
    plot!([Xval, Xval + dX], [Yval, Yval]; color = :black, lw = 2, arrow = :arrow, label = "")
    plot!([Xval, Xval], [Yval, Yval + dY]; color = :black, lw = 2, arrow = :arrow, label = "")
end
arrow_positions = [(0.25, 0.25), (0.75, 0.75), (0.35, 0.8), (0.8, 0.35), (0.02, 1.1), (1.1, 0.02)]
foreach(pos -> plot_xy_arrow!(pos...), arrow_positions)
plot!()
```
This also works as a form of simple stability analysis, where we can see how the solution moves *away* from the unstable steady state, and *to* the stable ones.
