# [Introduction to Catalyst and Julia for New Julia users](@id catalyst_for_new_julia_users)
The Catalyst tool for the modelling of chemical reaction networks is based in the Julia programming language. While experience in Julia programming is advantageous for using Catalyst, it is not necessary for accessing some of its basic features. This tutorial serves as an introduction to Catalyst for those unfamiliar with Julia, also introducing some basic Julia concepts. Anyone who plans on using Catalyst extensively is recommended to familiarise oneself more thoroughly with the Julia programming language. A collection of resources for learning Julia can be found [here](https://julialang.org/learning/), and a full documentation is available [here](https://docs.julialang.org/en/v1/).

Julia can be downloaded [here](https://julialang.org/downloads/).

*Users who are already familiar with Julia can skip to the [Using Catalyst](@ref using_catalyst) tutorial.*

## Basic Julia usage
On the surface, Julia has many similarities to languages like MATLAB, Python, and R.

*Values* can be assigned to *variables* through the use of a `=` sign. Values (possibly stored in variables) can be used for most basic computations.
```@example ex1
length = 2.0
width = 4.0
area = length*width
```

*Functions* take one or more inputs (enclosed by `()`) and return some output. E.g. the `min` function returns the minimum of two values
```@example ex1
min(1.0, 3.0)
```

A line of Julia code is not required to end with `;`, however, if it does, the output of that line is not displayed.
```julia
min(1.0, 3.0);
```

Each Julia variable has a specific *type*, designating what type of value it is. While not directly required to use Catalyst, this is useful to be aware of. To learn the type of a specific variable, use the `typeof` function. More information about types can be [found here](https://docs.julialang.org/en/v1/manual/types/).
```@example ex1
typeof(1.0)
```

Here, `Float64` denotes decimal-valued numbers. Integer-valued numbers instead are of the `Int64` type.
```@example ex1
typeof(1)
```

Finally, we note that the first time some code is run in Julia, it has to be *compiled*. However, this is only required once per Julia session. Hence, the second time the same code is run, it runs much faster. E.g. try running this line of code first one time, and then one additional time. You will note that the second run is much faster.
```@example ex1
rand(100, 100)^3.5;
```
(This code creates a random 100x100 matrix, and take it to teh power of 3.5)

This is useful to know when you e.g. declare, simulate, or plot, a Catalyst model. The first time you run a command there might be a slight delay. However, subsequent runs will execute much quicker. This holds even if you do minor adjustments before the second run (such as changing simulation initial conditions).

## Installing and activating packages
Except for some base Julia packages (such as `Pkg`, the package manager) that are available by default, Julia packages must be installed locally before they can be used. Most packages are registered with Julia, and can be added through the `Pkg.add("desired_package")` command (where `desired_package` is the name of the package you wish to install). We can thus install Catalyst:
```julia
using Pkg
Pkg.add("Catalyst")
```

Here, the command `using Pkg` is required to activate the package manager.

Next, we also wish to add the `DifferentialEquations` and `Plots` packages (for numeric simulation of models, and plotting, respectively).
```julia
Pkg.add("DifferentialEquations")
Pkg.add("Plots")
```
Once a package has been installed through the `Pkg.add` command, this command does not have to be repeated in further Julia sessions on the same machine.

Installing a Julia package is, however, not enough to use it. Before a package's features are used in a Julia session, it has to be loaded through the `using desired_package` command (where `desired_package` is the name of the package you wish to activate). This command has to be repeated whenever a Julia session is restarted.

We thus activate our three desired packages:

```@example ex1
using Catalyst
using DifferentialEquations
using Plots
```

For a more detailed introduction to Julia packages, please read [the Pkg documentation](https://docs.julialang.org/en/v1/stdlib/Pkg/).

## Simulating a basic Catalyst model
Now that we have some basic familiarity with Julia, and have installed and activated the required packages, we will create and simulate a basic chemical reaction network model through Catalyst.

Catalyst models are created through the `@reaction_network` *macro*. For more information on macros, please read [the Julia documentation on macros](https://docs.julialang.org/en/v1/manual/metaprogramming/#man-macros). This documentation is, however, rather advanced (and not required to use Catalyst). We instead recommend that you simply familiarise yourself with the Catalyst syntax, without studying in detail how macros work and what they are.

The `@reaction_network` command is followed by the `begin` keyword, which is followed by one line for each *reaction* of the model. Each reaction consists of a *reaction rate*, followed by the reaction itself. The reaction itself contains a set of *substrates* and a set of *products* (what is consumed and produced by the reaction, respectively). These are separated by a `-->` arrow. Finally, the model ends with the `end` keyword.

Here, we create a simple *birth-death* model, where a single species (*X*) is created at rate *b*, and degraded at rate *d*. The model is stored in the variable `rn`.

```@example ex1
rn = @reaction_network begin
    b, 0 --> X
    d, X --> 0
end
```

For more information on how to use the Catalyst model creator, please read [the corresponding documentation](https://docs.sciml.ai/Catalyst/stable/tutorials/dsl/).

Next, we wish to simulate our model. To do this, we need to provide some additional information to the simulator. This is
* The initial condition. That is the concentration or number of each species at the start of the simulation.
* The timespan. That is, the timeframe over which we wish to run the simulation.
* The parameter values. That is, the values of the model's parameters for this simulation.

The initial condition is given as a *Vector*. This is a type which collects several different values. To declare a vector, the values are specific within brackets, `[]`, and separated by `,`. Since we only have one species, the vector holds a single element. In this element, we set the value of *X* using the `:X => 1.0` syntax. Here, we first denote the name of the species, with a `:` pre-appended. Next follows `=>` and then the value of *X*. Since we wish to simulate the *concentration* of X over time, we will let the initial condition be decimal valued.
```@example ex1
u0 = [:X => 1.0]
```

The timespan sets the time point at which we start the simulation (typically `0.0` is used) and the final time point of the simulation. These are combined into a two-valued *Tuple*. Tuples are similar to vectors, but are enclosed by `()` and not `[]`. Again, we will let both time points be decimal valued.
```@example ex1
tspan = (0.0, 10.0)
```

Finally, the parameter values are, like the initial conditions, given as a vector. Since we have two parameters (*b* and *d*), the parameter vector has two values. We use a similar notation for setting the parameter values as the initial condition (first the parameter, then an arrow, then the value).
```@example ex1
params = [:b => 1.0, :d => 0.2]
```

Please read here for more information on [Vectors](https://docs.julialang.org/en/v1/manual/arrays/) and [Tuples](https://docs.julialang.org/en/v1/manual/types/#Tuple-Types).

Next, before we can simulate our model, we bundle all the required information together in a so-called `ODEProblem`. Note that the order in which the input (the model, the initial condition, the timespan, and the parameter values) is provided to the ODEProblem matters. E.g. the parameter values cannot be provided as the first argument, but have to be the last argument. Here, we save our `ODEProblem` in the `oprob` variable.


```@example ex1
oprob = ODEProblem(rn, u0, tspan, params)
```

We can now simulate our model. We do this by providing the `ODEProblem` to the `solve` function. We save the output to the `sol` variable.
```@example ex1
sol = solve(oprob)
```

Finally, we can plot the solution through the `plot` function.
```@example ex1
plot(sol)
```

Here, the plot shows the time evolution of the concentration of the species *X* from its initial condition.

For more information about the numerical simulation package, please see the [DifferentialEquation documentation](https://docs.sciml.ai/DiffEqDocs/stable/). For more information about the plotting package, please see the [Plots documentation](https://docs.juliaplots.org/stable/).

## Additional modelling example
To make this introduction more comprehensive, we here provide another example, using a more complicated model. In addition, instead of simulating our model as concentrations evolve over time, we will simulate the individual reaction events through the [Gillespie algorithm](https://en.wikipedia.org/wiki/Gillespie_algorithm). This is a way to add *noise* to our model.

Remember, unless we have restarted Julia, we do not need to activate our packages (through the `using` command) again.

This time, we will declare the so-called [SIR model for an infectious disease](https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology#The_SIR_model). Note that even if this model does not describe a set of chemical reactions, it can be modelled using the same dynamics. The model consists of 3 species:
* *S*, the amount of *susceptible* individuals.
* *I*, the amount of *infected* individuals.
* *R*, the amount of *recovered* (or *removed*) individuals.
It also has 2 reaction events:
* Infection, where a susceptible individual meets an infected individual, and also becomes infected.
* Recovery, where an infected individual recovers.
Each reaction is also associated with a specific rate (corresponding to a parameter).
* *b*, the infection rate.
* *k*, the recovery rate.

We declare the model using the `reaction_network` macro, and store it in the `sir_model` variable.
```@example ex1
sir_model = @reaction_network begin
    b, S + I --> 2I
    k, I --> R
end
```

Note that the first reaction contains two different substrates (separated by a `+` sign). While there is only a single product (*I*), two copies of *I* are produced. The *2* in front of the product *I* denotes this.

Next, we declare our initial condition, time span, and parameter values. Since we want to simulate the individual reaction events, that discretely change the state of our model, we want our initial conditions to be integer-valued. We will start with a mostly susceptible population, but where a single individual has been infected through some means.
```@example ex1
u0 = [:S => 50, :I => 1, :R => 0.0]
tspan = (0.0, 10.0)
params = [:b => 0.2, :k => 1.0]
```

Previously we have bundled this information into an `ODEProblem` (denoting a deterministic *ordinary differential equation*). Now we wish to simulate our model as a jump process (where each reaction event denotes a single jump in the state of the system). We do this by first creating a `DiscreteProblem`, and then using this as an input to a `JumpProblem`.
```@example ex1
dprob = DiscreteProblem(sir_model, u0, tspan, params)
jprob = JumpProblem(sir_model, dprob, Direct())
```

Again, the order in which the inputs are given to the `DiscreteProblem` and the `JumpProblem` is important. The last argument to the `JumpProblem` (`Direct()`) denotes which simulation method we wish to use. For now, we recommend the user simply use the `Direct()` option, and then consider alternative ones, see the [JumpProcesses.jl docs](https://docs.sciml.ai/JumpProcesses/stable/), when they are more familiar with modelling in Catalyst and Julia.

Finally, we can simulate our model using the `solve` function, and plot the solution using the `plot` function. Here, the `solve` function also has a second argument (`SSAStepper()`). This is a time stepping algorithm that calls the `Direct` solver to advance a simulation. Again, we recommend at this stage you simply use this option, and then explore further exactly what this means at a later stage.
```@example ex1
sol = solve(jprob, SSAStepper())
plot(sol)
```

**Exercise:** Try simulating the model several times. Note that the epidemic doesn't always take off, but sometimes dies out without spreading through the population. Try changing the infection rate (*b*), determining how this value affects the probability that the epidemic goes through the population.

## Feedback
If you are a new Julia user who has used this tutorial, and there was something you struggled with or would have liked to have explained better, please [raise an issue](https://github.com/SciML/Catalyst.jl/issues). That way, we can continue improving this tutorial.