# [Interactive Visualization of the Brusselator Model](@id interactive_brusselator)

```@contents
Pages = ["interactive_brusselator.md"]
Depth = 3
```

Catalyst can utilize the [GLMakie.jl](https://github.com/JuliaPlots/GLMakie.jl) package for creating interactive visualizations. This tutorial provides a step-by-step guide to creating an interactive visualization of the Brusselator model, building upon the basic [Brusselator](@ref basic_CRN_library_brusselator) example.

!!! note
    This tutorial assumes you have GLMakie.jl installed. If not, you can install it by running `using Pkg; Pkg.add("GLMakie")` in your Julia REPL.

## [Background: The Brusselator Model](@id brusselator_background)

The Brusselator is a theoretical model for autocatalytic reactions, proposed by Ilya Prigogine and his collaborators. It is described by the following set of differential equations:

```math
\begin{aligned}
\frac{dX}{dt} &= A + X^2Y - BX - X \\
\frac{dY}{dt} &= BX - X^2Y
\end{aligned}
```

Where X and Y are chemical species, and A and B are parameters. This system can exhibit oscillatory behavior for certain parameter values, making it an interesting subject for study in nonlinear dynamics and chemical kinetics.

## [Setting up the Brusselator model](@id setup_brusselator)

First, let's import the necessary packages and define our Brusselator model:

```julia
using Catalyst
using OrdinaryDiffEq
using GLMakie

# Define the Brusselator model
brusselator = @reaction_network begin
    A, ∅ → X
    1, 2X + Y → 3X
    B, X → Y
    1, X → ∅
end

# Initial parameter values and conditions
p = [:A => 1.0, :B => 3.0]
u0 = [:X => 1.0, :Y => 1.0]
tspan = (0.0, 50.0)

oprob = ODEProblem(brusselator, u0, tspan, p)

# Function to solve the ODE
function solve_brusselator(A, B, X0, Y0)
    p = [:A => A, :B => B]
    u0 = [:X => X0, :Y => Y0]
    prob = remake(oprob, p=p, u0=u0)
    solve(prob, Tsit5(), saveat = 0.1)
end
```

This code sets up our Brusselator model using Catalyst.jl's `@reaction_network` macro. We also define initial parameters, initial conditions, create an `ODEProblem`, and define a function to solve the ODE with given parameters.

## [Basic static plotting](@id basic_static_plotting)

Let's start by creating a basic plot of our Brusselator model:

```julia
# Create the main figure
fig = Figure(size = (800, 600), fontsize = 18)

# Create an axis for the plot
ax = Axis(fig[1, 1], 
    title = "Brusselator Model", 
    xlabel = "Time", 
    ylabel = "Concentration")

# Solve the ODE
sol = solve_brusselator(1.0, 3.0, 1.0, 1.0)

# Plot the solution
lines!(ax, sol.t, sol[:X], label = "X", color = :blue, linewidth = 3)
lines!(ax, sol.t, sol[:Y], label = "Y", color = :red, linewidth = 3)

# Add a legend
axislegend(ax, position = :rt)

# Display the figure
display(fig)
```

This will produce a basic time series plot of the Brusselator model:

![Basic Brusselator Plot](../assets/brusselator_basic_plot.svg)

The plot shows the concentrations of species X and Y over time. Notice the oscillatory behavior characteristic of the Brusselator model.

## [Adding interactivity](@id adding_interactivity)

Now, let's add interactivity to our plot using Observables and sliders. We'll build this up step by step.

### [Creating Observables](@id creating_observables)

Observables are a key concept in reactive programming and are central to how Makie.jl creates interactive visualizations.

```julia
# Create observables for parameters and initial conditions
A = Observable(1.0)
B = Observable(3.0)
X0 = Observable(1.0)
Y0 = Observable(1.0)
```

An Observable is a container for a value that can change over time. When the value changes, any dependent computations are automatically updated.

### [Adding sliders and connecting to Observables](@id adding_sliders)

Let's add sliders that will control our Observables:

```julia
# Create the main figure
fig = Figure(size = (800, 600), fontsize = 18)

# Create layout for plot and sliders
plot_layout = fig[1, 1] = GridLayout()
slider_layout = fig[2, 1] = GridLayout()

# Create sliders
slider_A = Slider(slider_layout[1, 1], range = 0.1:0.1:5.0, startvalue = 1.0)
slider_B = Slider(slider_layout[2, 1], range = 0.1:0.1:5.0, startvalue = 3.0)
slider_X0 = Slider(slider_layout[3, 1], range = 0.1:0.1:5.0, startvalue = 1.0)
slider_Y0 = Slider(slider_layout[4, 1], range = 0.1:0.1:5.0, startvalue = 1.0)

# Add labels for sliders
Label(slider_layout[1, 1, Left()], "A")
Label(slider_layout[2, 1, Left()], "B")
Label(slider_layout[3, 1, Left()], "X₀")
Label(slider_layout[4, 1, Left()], "Y₀")

# Connect sliders to observables
connect!(A, slider_A.value)
connect!(B, slider_B.value)
connect!(X0, slider_X0.value)
connect!(Y0, slider_Y0.value)
```

These sliders allow us to interactively change the parameters A and B, as well as the initial conditions X₀ and Y₀.

### [Creating a reactive plot](@id reactive_plot)

Now, let's create a plot that reacts to changes in our sliders:

```julia
# Create an axis for the plot
ax = Axis(plot_layout[1, 1], 
    title = "Brusselator Model", 
    xlabel = "Time", 
    ylabel = "Concentration")

# Create an observable for the solution
solution = @lift(solve_brusselator($A, $B, $X0, $Y0))

# Plot the solution
lines!(ax, lift(sol -> sol.t, solution), lift(sol -> sol[:X], solution), label = "X", color = :blue, linewidth = 3)
lines!(ax, lift(sol -> sol.t, solution), lift(sol -> sol[:Y], solution), label = "Y", color = :red, linewidth = 3)

# Add a legend
axislegend(ax, position = :rt)

# Display the figure
display(fig)
```

The resulting figure should look like this:

![Interactive Brusselator Plot](../assets/brusselator_interactive_plot.svg)

This plot will now update in real-time as you move the sliders, allowing for interactive exploration of the Brusselator's behavior under different conditions.

## [Enhanced visualization with phase plot](@id adding_phase_plot)

To gain more insight into the system's behavior, let's enhance our visualization by adding a phase plot:

```julia
# Create the main figure
fig = Figure(size = (1200, 800), fontsize = 18)

# Create main layout: plots on top, sliders at bottom
plot_grid = fig[1, 1] = GridLayout()
slider_grid = fig[2, 1] = GridLayout()

# Create sub-grids for plots
time_plot = plot_grid[1, 1] = GridLayout()
phase_plot = plot_grid[1, 2] = GridLayout()

# Create axes for the time series plot and phase plot
ax_time = Axis(time_plot[1, 1], 
    title = "Brusselator Model - Time Series", 
    xlabel = "Time", 
    ylabel = "Concentration")

ax_phase = Axis(phase_plot[1, 1], 
    title = "Brusselator Model - Phase Plot", 
    xlabel = "X", 
    ylabel = "Y")

# Create sub-grids for sliders
param_grid = slider_grid[1, 1] = GridLayout()
ic_grid = slider_grid[1, 2] = GridLayout()

# Create sliders with labels and group titles
Label(param_grid[1, 1:2], "Parameters", fontsize = 22)
slider_A = Slider(param_grid[2, 2], range = 0.1:0.1:5.0, startvalue = 1.0)
slider_B = Slider(param_grid[3, 2], range = 0.1:0.1:5.0, startvalue = 3.0)
Label(param_grid[2, 1], "A")
Label(param_grid[3, 1], "B")

Label(ic_grid[1, 1:2], "Initial Conditions", fontsize = 22)
slider_X0 = Slider(ic_grid[2, 2], range = 0.1:0.1:5.0, startvalue = 1.0)
slider_Y0 = Slider(ic_grid[3, 2], range = 0.1:0.1:5.0, startvalue = 1.0)
Label(ic_grid[2, 1], "X₀")
Label(ic_grid[3, 1], "Y₀")

# Connect sliders to observables
A = Observable(1.0)
B = Observable(3.0)
X0 = Observable(1.0)
Y0 = Observable(1.0)
connect!(A, slider_A.value)
connect!(B, slider_B.value)
connect!(X0, slider_X0.value)
connect!(Y0, slider_Y0.value)

# Create an observable for the solution
solution = @lift(solve_brusselator($A, $B, $X0, $Y0))

# Plot the time series
lines!(ax_time, lift(sol -> sol.t, solution), lift(sol -> sol[:X], solution), label = "X", color = :blue, linewidth = 3)
lines!(ax_time, lift(sol -> sol.t, solution), lift(sol -> sol[:Y], solution), label = "Y", color = :red, linewidth = 3)

# Plot the phase plot
phase_plot_obj = lines!(ax_phase, lift(sol -> sol[:X], solution), lift(sol -> sol[:Y], solution), 
                        color = lift(sol -> sol.t, solution), colormap = :viridis)

# Add a colorbar for the phase plot
Colorbar(phase_plot[1, 2], phase_plot_obj, label = "Time")

# Add legends
axislegend(ax_time, position = :rt)

# Adjust layout
colgap!(plot_grid, 20)
rowgap!(fig.layout, 20)
colgap!(param_grid, 10)
colgap!(ic_grid, 10)

# Display the figure
display(fig)
```

This will create a visualization with both time series and phase plots:

![Interactive Brusselator Plot with Time Series and Phase Plot](../assets/brusselator_interactive_plot_full.svg)


## [Performance considerations](@id interactive_performance_considerations)

When working with interactive visualizations, especially for more complex models, performance can become an issue. Here are some tips to improve performance:

1. Use appropriate ODE solvers. For stiff systems, consider using `Rodas4()` instead of `Tsit5()`.
2. Adjust the `saveat` parameter in the `solve` function to reduce the number of points plotted.
3. For very large systems, consider using `GLMakie.jl` instead of `Makie.jl` for better performance.

## [Common plotting options](@id common_makie_plotting_options)

Various plotting options can be provided as optional arguments to the `lines!` command. Common options include:
- `linewidth` or `lw`: Determine plot line widths.
- `linestyle`: Determines plot line style.
- `color`: Determines the line colours.
- `label`: Determines label texts.

For example:

```julia
lines!(ax_time, lift(sol -> sol.t, solution), lift(sol -> sol[:X], solution), 
       label = "X", color = :blue, linewidth = 2, linestyle = :dash)
```

## [Extending the interactive visualization](@id extending_interactive_visualization)

You can further extend this visualization by:
- Adding other interactive elements, such as buttons or dropdowns to control different aspects of the simulation or visualization.
- Adding additonal axes to the plot.


