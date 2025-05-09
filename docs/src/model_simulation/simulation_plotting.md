# [Simulation Plotting](@id simulation_plotting)
Catalyst uses the [Plots.jl](https://github.com/JuliaPlots/Plots.jl) package for performing all plots. This section provides a brief summary of some useful plotting options, while [Plots.jl's documentation](https://docs.juliaplots.org/stable/) provides a more throughout description of how to tune your plots.

!!! note
    [Makie.jl](https://github.com/MakieOrg/Makie.jl) is a popular alternative to the Plots.jl package. While it is not used within Catalyst's documentation, it is worth considering (especially for users interested in interactivity, or increased control over their plots).

## [Common plotting options](@id simulation_plotting_options)
Let us consider the oscillating [Brusselator](@ref basic_CRN_library_brusselator) model. We have previously shown how model simulation solutions can be plotted using the `plot` function. Here we plot an ODE simulation from the Brusselator:
```@example simulation_plotting
using Catalyst, OrdinaryDiffEqDefault, Plots

brusselator = @reaction_network begin
    A, ∅ --> X
    1, 2X + Y --> 3X
    B, X --> Y
    1, X --> ∅
end
u0 = [:X => 1.0, :Y => 0.0]
tspan = (0.0, 50.0)
ps = [:A => 1.0, :B => 4.0]

oprob = ODEProblem(brusselator, u0, tspan, ps)
sol = solve(oprob)
plot(sol)
```

Various plotting options can be provided as optional arguments to the `plot` command. Common options include:
- `lw`: Determine plot line widths.
- `la`: Determine plot line's transparency (at `la = 0.0` lines are fully transparent, i.e. not visible).
- `linestyle`: Determines plot line style.
- `color`: Determines the line colours.
- `legend`: Determines the position of the legend/labels.
- `label`: Determines label texts.
- `xguide`, `yguide`: Determines x and y axis labels.
- `title`: Determines plot title.
- `legendfontsize`, `guidefontsize`, `titlefontsize`: Determines the font size of the labels, x and y guides, and title, respectively.

Here, we re-plot our simulations, utilising some of these options (`legend = :none` is used to disable the legends).
```@example simulation_plotting
plot(sol; lw = 4, linestyle = :dash, color = :green, xguide = "Time", yguide = "Concentration", guidefontsize = 14)
```
Note that, by default, Catalyst uses `xguide = "t"`. Here, however, we modify this to `xguide = "Time"`. We also note that the `color = :green` change both lines' colours to green. To set different colours for each line, we provide these as *a vector without `,` in-between elements* (in Julia interpreted as a matrix with its first dimension equal to `1`):
```@example simulation_plotting
plot(sol; lw = 4, color = [:green :purple])
```
A full list of available colours can be found [here](https://juliagraphics.github.io/Colors.jl/stable/namedcolors/). A full list of possible plotting options can be found [here](https://docs.juliaplots.org/stable/attributes/) (look at the list of various plot attributes, e.g. "Series Attributes"). if there is some option(s) you intend to use multiple times, you can call the `default` function using these, in which case they will be used for all subsequent plots. E.g. here:
```@example simulation_plotting
default(framestyle = :box, grid = false)
```
we designate a box-style frame, and remove the faint background grid, for all subsequent plots in this tutorial.

A useful option unique to Catalyst (and other DifferentialEquations.jl-based) plots is `idxs`. Its input is a vector, listing all the species (or quantities) that should be plotted. I.e.
```@example simulation_plotting
plot(sol; idxs = [:X])
```
can be used to plot `X` only. When only a single argument is given, the vector form is unnecessary (e.g. `idxs = :X` could have been used instead). If symbolic species representation is used, this can be used to designate any algebraic expression(s) that should be plotted. E.g. here we plot the total concentration of $X + Y$ throughout the simulation:
```@example simulation_plotting
plot(sol; idxs = brusselator.X + brusselator.Y)
```

## [Multi-plot plots](@id simulation_plotting_options_subplots)
It is possible to save plots in variables. These can then be used as input to the `plot` command. Here, the plot command can be used to create plots containing multiple plots (by providing multiple inputs). E.g. here we plot the concentration of `X` and `Y` in separate subplots:
```@example simulation_plotting
plt_X = plot(sol; idxs = [:X])
plt_Y = plot(sol; idxs = [:Y])
plot(plt_X, plt_Y)
```

When working with subplots, the [`layout`](https://docs.juliaplots.org/latest/layouts/) and [`size`](https://docs.juliaplots.org/latest/generated/attributes_plot/) options are typically useful. Here we use `layout` to put the first plot above the second one, and `size` to reshape the plot size:
```@example simulation_plotting
plot(plt_X, plt_Y; layout = (2,1), size = (700,500))
```

## [Saving plots](@id simulation_plotting_options_saving)
Once a plot has been saved to a variable, the `savefig` function can be used to save it to a file. Here we save our Brusselator plot simulation (the first argument) to a file called "saved_plot.png" (the second argument):
```@example simulation_plotting
plt = plot(sol)
savefig(plt, "saved_plot.png")
rm("saved_plot.png") # hide
```
The plot file type is automatically determined from the extension (if none is given, a .png file is created).

## [Phase-space plots](@id simulation_plotting_options_phasespace)
By default, simulations are plotted as species concentrations over time. However, [phase space](https://en.wikipedia.org/wiki/Phase_space#:~:text=In%20dynamical%20systems%20theory%20and,point%20in%20the%20phase%20space.) plots are also possible. This is done by designating the axis arguments using the `idxs` option, but providing them as a tuple. E.g. here we plot our simulation in $X-Y$ space:
```@example simulation_plotting
plot(sol; idxs = (:X, :Y))
```
