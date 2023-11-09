# [Fitting Parameters for an Oscillatory System](@id parameter_estimation)
In this examples we will use [Optimization.jl](https://github.com/SciML/Optimization.jl) to fit the parameters of an oscillatory system (the Brusselator). Here, special consideration is taken to avoid reaching a local minimum. Instead of fitting the entire time series directly, we will start with fitting parameter values for the first period, and then use those as an initial guess for fitting the next. As we complete this procedure (which can be advantageous for oscillatory systems) we reach a global optimum.

First, we fetch the required packages.
```@example pe1
using Catalyst
using OrdinaryDiffEq
using Optimization
using OptimizationOptimisers # Required for the ADAM optimizer.
using SciMLSensitivity # Required for `Optimization.AutoZygote()` automatic differentiation option.
```

Next, we declare our model, the Brusselator oscillator.
```@example pe1
brusselator = @reaction_network begin
    A, ∅ --> X
    1, 2X + Y --> 3X
    B, X --> Y
    1, X --> ∅
end
p_real = [:A => 1., :B => 2.]
nothing # hide
```

We simulate our model, and from the simulation generate sampled data points
(to which we add noise). We will use this data to fit the parameters of our model.
```@example pe1
u0 = [:X => 1.0, :Y => 1.0]
tspan = (0.0, 30.0)

sample_times = range(tspan[1]; stop = tspan[2], length = 100)
prob = ODEProblem(brusselator, u0, tspan, p_real)
sol_real = solve(prob, Rosenbrock23(); tstops = sample_times)
sample_vals = Array(sol_real(sample_times))
sample_vals .*= (1 .+ .1 * rand(Float64, size(sample_vals)) .- .05)
nothing   # hide
```

We can plot the real solution, as well as the noisy samples.
```@example pe1
using Plots
default(; lw = 3, framestyle = :box, size = (800, 400))

plot(sol_real; legend = nothing, color = [:darkblue :darkred])
scatter!(sample_times, sample_vals'; color = [:blue :red], legend = nothing)
```

Next, we create a function to fit the parameters using the `ADAM` optimizer. For
a given initial estimate of the parameter values, `pinit`, this function will
fit parameter values, `p`, to our data samples. We use `tend` to indicate the
time interval over which we fit the model.
```@example pe1
function optimise_p(pinit, tend)
    function loss(p, _)
        newtimes = filter(<=(tend), sample_times)
        newprob = remake(prob; tspan = (0.0, tend), p = p)
        sol = Array(solve(newprob, Rosenbrock23(); saveat = newtimes))
        loss = sum(abs2, sol .- sample_vals[:, 1:size(sol,2)])
        return loss, sol
    end

    # optimize for the parameters that minimize the loss
    optf = OptimizationFunction(loss, Optimization.AutoZygote())
    optprob = OptimizationProblem(optf, pinit)
    sol = solve(optprob, ADAM(0.1); maxiters = 100)

    # return the parameters we found
    return sol.u
end
nothing # hide
```

Next, we will fit a parameter set to the data on the interval `(0, 10)`.
```@example pe1
p_estimate = optimise_p([5.0, 5.0], 10.0)
```

We can compare this to the real solution, as well as the sample data
```@example pe1
newprob = remake(prob; tspan = (0., 10.), p = p_estimate)
sol_estimate = solve(newprob, Rosenbrock23())
plot(sol_real; color = [:blue :red], label = ["X real" "Y real"], linealpha = 0.2)
scatter!(sample_times, sample_vals'; color = [:blue :red],
         label = ["Samples of X" "Samples of Y"], alpha = 0.4)
plot!(sol_estimate; color = [:darkblue :darkred], linestyle = :dash,
                    label = ["X estimated" "Y estimated"], xlimit = tspan)
```

Next, we use this parameter estimate as the input to the next iteration of our
fitting process, this time on the interval `(0, 20)`.
```@example pe1
p_estimate = optimise_p(p_estimate, 20.)
newprob = remake(prob; tspan = (0., 20.), p = p_estimate)
sol_estimate = solve(newprob, Rosenbrock23())
plot(sol_real; color = [:blue :red], label = ["X real" "Y real"], linealpha = 0.2)
scatter!(sample_times, sample_vals'; color = [:blue :red],
         label = ["Samples of X" "Samples of Y"], alpha = 0.4)
plot!(sol_estimate; color = [:darkblue :darkred], linestyle = :dash,
                    label = ["X estimated" "Y estimated"], xlimit = tspan)
```

Finally, we use this estimate as the input to fit a parameter set on the full
time interval of the sampled data.
```@example pe1
p_estimate = optimise_p(p_estimate, 30.0)

newprob = remake(prob; tspan = (0., 30.0), p = p_estimate)
sol_estimate = solve(newprob, Rosenbrock23())
plot(sol_real; color = [:blue :red], label = ["X real" "Y real"], linealpha = 0.2)
scatter!(sample_times, sample_vals'; color = [:blue :red],
        label = ["Samples of X" "Samples of Y"], alpha = 0.4)
plot!(sol_estimate; color = [:darkblue :darkred], linestyle = :dash,
                    label = ["X estimated" "Y estimated"], xlimit = tspan)
```

The final parameter estimate is then
```@example pe1
p_estimate
```
which is close to the actual parameter set of `[1.0, 2.0]`.

## Why we fit the parameters in iterations
As previously mentioned, the reason we chose to fit the model on a smaller interval to begin with, and
then extend the interval, is to avoid getting stuck in a local minimum. Here
specifically, we chose our initial interval to be smaller than a full cycle of
the oscillation. If we had chosen to fit a parameter set on the full interval
immediately we would have received an inferior solution.
```@example pe1
p_estimate = optimise_p([5.0,5.0], 30.0)

newprob = remake(prob; tspan = (0.0,30.0), p = p_estimate)
sol_estimate = solve(newprob, Rosenbrock23())
plot(sol_real; color = [:blue :red], label = ["X real" "Y real"], linealpha = 0.2)
scatter!(sample_times,sample_vals'; color = [:blue :red],
         label = ["Samples of X" "Samples of Y"], alpha = 0.4)
plot!(sol_estimate; color = [:darkblue :darkred], linestyle = :dash,
                    label = ["X estimated" "Y estimated"], xlimit = tspan)
```
