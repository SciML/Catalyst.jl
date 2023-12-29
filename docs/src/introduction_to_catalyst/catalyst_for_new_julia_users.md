# [Introduction to Catalyst and Julia for New Julia users](@id catalyst_for_new_julia_users)
The Catalyst tool for the modelling of chemical reaction networks is based in the Julia programming language. While experience in Julia programming is advantageous for using Catalyst, it is not necessary for accessing most of its basic features. This tutorial serves as an introduction to Catalyst for those unfamiliar with Julia, while also introducing some basic Julia concepts. Anyone who plans on using Catalyst extensively is recommended to familiarise oneself more thoroughly with the Julia programming language. A collection of resources for learning Julia can be found [here](https://julialang.org/learning/), and a full documentation is available [here](https://docs.julialang.org/en/v1/). A more practical (but also extensive) guide to Julia programming can be found [here](https://modernjuliaworkflows.github.io/writing/).

Julia can be downloaded [here](https://julialang.org/downloads/). Generally, it is recommended to use the [*juliaup*](https://github.com/JuliaLang/juliaup) tool to install and update Julia. Furthermore, *Visual Studio Code* is a good IDE with [extensive Julia support](https://code.visualstudio.com/docs/languages/julia), and a good default choice.

*Users who are already familiar with Julia can skip to the [Introduction to Catalyst](@ref introduction_to_catalyst) tutorial.*

## Basic Julia usage
On the surface, Julia has many similarities to languages like MATLAB, Python, and R.

*Values* can be assigned to *variables* through `=` sign. Values (possibly stored in variables) can be used for most basic computations.
```@example ex1
length = 2.0
width = 4.0
area = length * width
```

*Functions* take one or more inputs (enclosed by `()`) and return some output. E.g. the `min` function returns the minimum of two values.
```@example ex1
min(1.0, 3.0)
```
Julia has a specific *help mode*, which can be [queried for information about any function](https://docs.julialang.org/en/v1/stdlib/REPL/#Help-mode) (including those defined by Catalyst).

Each Julia variable has a specific *type*, designating what type of value it contains. While not directly required to use Catalyst, this is useful to be aware of. To learn the type of a specific variable, use the `typeof` function. More information about types can be [found here](https://docs.julialang.org/en/v1/manual/types/).
```@example ex1
typeof(1.0)
```
Here, `Float64` denotes decimal-valued numbers. Integer-valued numbers instead have the `Int64` type.
```@example ex1
typeof(1)
```
There exists a large number of Julia types (with even more being defined by various packages). Additional examples include `String`s (defined by enclosing text within `" "`):
```@example ex1
"Hello world!"
```
and `Symbol`s (defined by pre-appending an expression with `:`):
```@example ex1
:Julia
```

Finally, we note that the first time some code is run in Julia, it has to be [*compiled*](https://en.wikipedia.org/wiki/Just-in-time_compilation). However, this is only required once per Julia session. Hence, the second time the same code is run, it runs much faster. E.g. try running this line of code first one time, and then one additional time. You will note that the second run is much faster.
```@example ex1
rand(100, 100)^3.5
nothing # hide
```
(This code creates a random 100x100 matrix, and takes it to the power of 3.5)

This is useful to know when you e.g. declare, simulate, or plot, a Catalyst model. The first time you run a command there might be a slight delay. However, subsequent runs will be much quicker. This holds even if you make minor adjustments before the second run (such as changing simulation initial conditions).

## [Installing and activating packages](@id catalyst_for_new_julia_users_packages_intro)
Due to its native package manager (Pkg), and a registry of almost all packages of relevancy, package management in Julia is unusually easy. Here, we will briefly describe how to install and activate Catalyst (and two additional packages relevant to this tutorial).

To import a Julia package into a session, you can use the `using PackageName` command (where `PackageName` is the name of the package you wish to import). However, before you can do so, it must first be installed on your computer. This is done through the `Pkg.add("PackageName")` command:
```julia
using Pkg
Pkg.add("Catalyst")
using Catalyst
```
Here, the Julia package manager package (`Pkg`) is by default installed on your computer when Julia is installed, and can be activated directly. Next, we also wish to install the `DifferentialEquations` and `Plots` packages (for numeric simulation of models, and plotting, respectively).
```julia
using Pkg
Pkg.activate(".")
```
Once a package has been installed through the `Pkg.add` command, this command does not have to be repeated if we restart our Julia session. We can now import all three packages into our current session with:
```@example ex2
using Catalyst
```
This will only make Catalyst available for the current Julia session. If you exit Julia, you will have to run `using Catalyst` again to use its features (however, `Pkg.add("Catalyst")` does not need to be rerun). Next, we wish to install the `DifferentialEquations` and `Plots` packages (for numeric simulation of models, and plotting, respectively):
```julia
Pkg.add("DifferentialEquations")
Pkg.add("Plots")
```
Once a package has been installed through the `Pkg.add` command, this command does not have to be repeated if we restart our Julia session. We can now import all three packages into our current session with:
```@example ex2
using Catalyst
using DifferentialEquations
using Plots
```
Here, if we restart Julia, these `using` commands *must be rerun*.

A more comprehensive (but still short) introduction to package management in Julia is provided at [the end of this documentation page](@ref catalyst_for_new_julia_users_packages). It contains some useful information and is hence highly recommended reading. For a more detailed introduction to Julia package management, please read [the Pkg documentation](https://docs.julialang.org/en/v1/stdlib/Pkg/).

## Simulating a basic Catalyst model
Now that we have some basic familiarity with Julia, and have installed and imported the required packages, we will create and simulate a basic chemical reaction network model using Catalyst.

Catalyst models are created through the `@reaction_network` *macro*. For more information on macros, please read [the Julia documentation on macros](https://docs.julialang.org/en/v1/manual/metaprogramming/#man-macros). This documentation is, however, rather advanced (and not required to use Catalyst). We instead recommend that you simply familiarise yourself with the Catalyst syntax, without studying in detail how macros work and what they are.

The `@reaction_network` command is followed by the `begin` keyword, which is followed by one line for each *reaction* of the model. Each reaction consists of a *reaction rate*, followed by the reaction itself. The reaction contains a set of *substrates* and a set of *products* (what is consumed and produced by the reaction, respectively). These are separated by a `-->` arrow. Finally, the model ends with the `end` keyword.

Here, we create a simple *birth-death* model, where a single species ($X$) is created at rate $b$, and degraded at rate $d$. The model is stored in the variable `rn`.
```@example ex2
rn = @reaction_network begin
 b, 0 --> X
 d, X --> 0
end
```
For more information on how to use the Catalyst model creator (also known as *the Catalyst DSL*), please read [the corresponding documentation](https://docs.sciml.ai/Catalyst/stable/catalyst_functionality/dsl_description/).

Next, we wish to simulate our model. To do this, we need to provide some additional information to the simulator. This is
* The initial condition. That is, the concentration (or copy numbers) of each species at the start of the simulation.
* The time span. That is, the time frame over which we wish to run the simulation.
* The parameter values. That is, the values of the model's parameters for this simulation.

The initial condition is given as a *Vector*. This is a type which collects several different values. To declare a vector, the values are specific within brackets, `[]`, and separated by `,`. Since we only have one species, the vector holds a single element. In this element, we set the value of $X$ using the `:X => 1.0` syntax. Here, we first denote the name of the species (with a `:` pre-appended, which creates a `Symbol`), next follows a `=>` and then the value of $X$. Since we wish to simulate the *concentration* of X over time, we will let the initial condition be decimal valued.
```@example ex2
u0 = [:X => 1.0]
```

The timespan sets the time point at which we start the simulation (typically `0.0` is used) and the final time point of the simulation. These are combined into a two-valued *tuple*. Tuples are similar to vectors, but are enclosed by `()` and not `[]`. Again, we will let both time points be decimal valued.
```@example ex2
tspan = (0.0, 10.0)
```

Finally, the parameter values are, like the initial conditions, given in a vector. Since we have two parameters ($b$ and $d$), the parameter vector has two values. We use a similar notation for setting the parameter values as the initial condition (first the parameter, then an arrow, then the value).
```@example ex2
params = [:b => 1.0, :d => 0.2]
```

Please read here for more information on [vectors](https://docs.julialang.org/en/v1/manual/arrays/) and [tuples](https://docs.julialang.org/en/v1/manual/types/#Tuple-Types).

Next, before we can simulate our model, we bundle all the required information together in a so-called `ODEProblem`. Note that the order in which the input (the model, the initial condition, the timespan, and the parameter values) is provided to the `ODEProblem` matters. E.g. the parameter values cannot be provided as the first argument, but have to be the fourth argument. Here, we save our `ODEProblem` in the `oprob` variable.
```@example ex2
oprob = ODEProblem(rn, u0, tspan, params)
```

We can now simulate our model. We do this by providing the `ODEProblem` to the `solve` function. We save the output to the `sol` variable.
```@example ex2
sol = solve(oprob)
```

Finally, we can plot the solution through the `plot` function.
```@example ex2
plot(sol)
```

Here, the plot shows the time evolution of the concentration of the species $X$ from its initial condition.

For more information about the numerical simulation package, please see the [DifferentialEquations documentation](https://docs.sciml.ai/DiffEqDocs/stable/). For more information about the plotting package, please see the [Plots documentation](https://docs.juliaplots.org/stable/).

## Additional modelling example
To make this introduction more comprehensive, we here provide another example, using a more complicated model. Instead of simulating our model as concentrations evolve over time, we will now simulate the individual reaction events through the [Gillespie algorithm](https://en.wikipedia.org/wiki/Gillespie_algorithm) (a common approach for adding *noise* to models).

Remember (unless we have restarted Julia) we do not need to activate our packages (through the `using` command) again.

This time, we will declare a so-called [SIR model for an infectious disease](https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology#The_SIR_model). Note that even if this model does not describe a set of chemical reactions, it can be modelled using the same framework. The model consists of 3 species:
* $S$, the amount of *susceptible* individuals.
* $I$, the amount of *infected* individuals.
* $R$, the amount of *recovered* (or *removed*) individuals.

It also has 2 reaction events:
* Infection, where a susceptible individual meets an infected individual and also becomes infected.
* Recovery, where an infected individual recovers from the infection.

Each reaction is also associated with a specific rate (corresponding to a parameter).
* *b*, the infection rate.
* *k*, the recovery rate.

We declare the model using the `@reaction_network` macro, and store it in the `sir_model` variable.
```@example ex2
sir_model = @reaction_network begin
 b, S + I --> 2I
 k, I --> R
end
```
Note that the first reaction contains two different substrates (separated by a `+` sign). While there is only a single product (*I*), two copies of *I* are produced. The *2* in front of the product *I* denotes this.

Next, we declare our initial condition, time span, and parameter values. Since we want to simulate the individual reaction events that discretely change the state of our model, we want our initial conditions to be integer-valued. We will start with a mostly susceptible population, but to which a single infected individual has been introduced.
```@example ex2
u0 = [:S => 50, :I => 1, :R => 0]
tspan = (0.0, 10.0)
params = [:b => 0.2, :k => 1.0]
nothing # hide
```

Previously we have bundled this information into an `ODEProblem` (denoting a deterministic *ordinary differential equation*). Now we wish to simulate our model as a jump process (where each reaction event corresponds to a single jump in the state of the system). We do this by first creating a `DiscreteProblem`, and then using this as an input to a `JumpProblem`.
```@example ex2
dprob = DiscreteProblem(sir_model, u0, tspan, params)
jprob = JumpProblem(sir_model, dprob, Direct())
```
Again, the order in which the inputs are given to the `DiscreteProblem` and the `JumpProblem` is important. The last argument to the `JumpProblem` (`Direct()`) denotes which simulation method we wish to use. For now, we recommend that users simply use the `Direct()` option, and then consider alternative ones (see the [JumpProcesses.jl docs](https://docs.sciml.ai/JumpProcesses/stable/)) when they are more familiar with modelling in Catalyst and Julia.

Finally, we can simulate our model using the `solve` function, and plot the solution using the `plot` function. Here, the `solve` function also has a second argument (`SSAStepper()`). This is a time-stepping algorithm that calls the `Direct` solver to advance a simulation. Again, we recommend at this stage you simply use this option, and then explore exactly what this means at a later stage.
```@example ex2
sol = solve(jprob, SSAStepper())
sol = solve(jprob, SSAStepper(); seed=1234) # hide
plot(sol)
```

**Exercise:** Try simulating the model several times. Note that the epidemic doesn't always take off, but sometimes dies out without spreading through the population. Try changing the infection rate (*b*), determining how this value affects the probability that the epidemic goes through the population.

## [Package management in Julia](@id catalyst_for_new_julia_users_packages)
We have previously introduced how to install and activate Julia packages. While this is enough to get started with Catalyst, for long-term users, there are some additional considerations for a smooth experience. These are described here.

### [Setting up a new Julia environment](@id catalyst_for_new_julia_users_packages_environments)
Whenever you run Julia, it will run in a specific *environment*. You can specify any folder on your computer as a Julia environment. Some modes of running Julia will automatically use the environment corresponding to the folder you start Julia in. Others (or if you start Julia in a folder without an environment), will use your *default* environment. In these cases you can, during your session, switch to another environment. While it is possible to not worry about environments (and always use the default one), this can lead to long-term problems as more packages are installed.

To activate your current folder as an environment, run the following commands:
```julia
using Pkg
Pkg.activate(".")
```
This will:
1. If your current folder (which can be displayed using the `pwd()` command) is not designated as a possible Julia environment, designate it as such.
2. Switch your current Julia session to use the current folder's environment.

!!! note
 If you check any folder which has been designated as a Julia environment, it contains a Project.toml and a Manifest.toml file. These store all information regarding the corresponding environment. For non-advanced users, it is recommended to never touch these files directly (and instead do so using various functions from the Pkg package, the important ones which are described in the next two subsections).

### [Installing and importing packages in Julia](@id catalyst_for_new_julia_users_packages_installing)
Package installation and import have been described [previously](@ref catalyst_for_new_julia_users_packages_intro). However, for the sake of this extended tutorial, let us repeat the description by demonstrating how to install the [Latexify.jl](https://github.com/korsbo/Latexify.jl) package (which enables e.g. displaying Catalyst models in Latex format). First, we import the Julia Package manager ([Pkg](https://github.com/JuliaLang/Pkg.jl)) (which is required to install Julia packages):
```@example ex3
using Pkg
```
Latexify is a registered package, so it can be installed directly using:
 ```julia
Pkg.add("Latexify")
```
Finally, to import Latexify into our current Julia session we use:
```@example ex3
using Latexify
```
Here, `using Latexify` must be rerun whenever you restart a Julia session. However, you only need to run `Pkg.add("Latexify")` once to install it on your computer (but possibly additional times to add it to new environments, see the next section).

### [Why environments are important](@id catalyst_for_new_julia_users_packages_environment_importance)
We have previously described how to set up new Julia environments, how to install Julia packages, and how to import them into a current session. Let us say that you were to restart Julia in a new folder and activate this as a separate environment. If you then try to import Latexify through `using Latexify` you will receive an error claiming that Latexify was not found. The reason is that the `Pkg.add("Latexify")` command actually carries out two separate tasks:
1. If Latexify is not already installed on your computer, install it.
2. Add Latexify as an available package to your current environment.

Here, while Catalyst has previously been installed on your computer, it has not been added to the new environment you created. To do so, simply run 
```julia
using Pkg
Pkg.add("Latexify")
```
after which Catalyst can be imported through `using Latexify`. You can get a list of all packages available in your current environment using:
```julia
Pkg.status()
```

So, why is this required, and why cannot we simply import any package installed on our computer? The reason is that most packages depend on other packages, and these dependencies may be restricted to only specific versions of these packages. This creates complicated dependency graphs that restrict what versions of what packages are compatible with each other. When you use `Pkg.add("PackageName")`, only a specific version of that package is actually added (the latest possible version as permitted by the dependency graph). Here, Julia environments both define what packages are available *and* their respective versions (these versions are also displayed by the `Pkg.status()` command). By doing this, Julia can guarantee that the packages (and their versions) specified in an environment are compatible with each other.

The reason why all this is important is that it is *highly recommended* to, for each project, define a separate environment. To these, only add the required packages. General-purpose environments with a large number of packages often, in the long term, produce package incompatibility issues. While these might not prevent you from installing all desired package, they often mean that you are unable to use the latest version of some packages.

!!! note
 A not-infrequent cause for reported errors with Catalyst (typically the inability to replicate code in tutorials) is package incompatibilities in large environments preventing the latest version of Catalyst from being installed. Hence, whenever an issue is encountered, it is useful to run `Pkg.status()` to check whenever the latest version of Catalyst is being used.

Some additional useful Pkg commands are:
- `Pk.rm("PackageName")` removes a package from the current environment.
- `Pkg.update("PackageName")`: updates the designated package.
- `Pkg.update()`: updates all packages.

!!! note
 A useful feature of Julia's environment system is that enables the exact definition of what packages and versions were used to execute a script. This supports e.g. reproducibility in academic research. Here, by providing the corresponding Project.toml and Manifest.toml files, you can enable someone to reproduce the exact program used to perform some set of analyses.


---
## Feedback
If you are a new Julia user who has used this tutorial, and there was something you struggled with or would have liked to have explained better, please [raise an issue](https://github.com/SciML/Catalyst.jl/issues). That way, we can continue improving this tutorial. The same goes for any part of the Catalyst documentation: It is written to help new users understand how to use the package, and if it is not doing so successfully we would like to know!

---
## References
[^1]: [Jeff Bezanson, Alan Edelman, Stefan Karpinski, Viral B. Shah, *Julia: A Fresh Approach to Numerical Computing*, SIAM Review (2017).](https://epubs.siam.org/doi/abs/10.1137/141000671)
[^2]: [Torkel E. Loman, Yingbo Ma, Vasily Ilin, Shashi Gowda, Niklas Korsbo, Nikhil Yewale, Chris Rackauckas, Samuel A. Isaacson, *Catalyst: Fast and flexible modeling of reaction networks*, PLOS Computational Biology (2023).](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011530)