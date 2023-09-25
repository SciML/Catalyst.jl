# [Parameter Fitting for ODEs using PEtab.jl](@id petab_parameter_fitting)
The [PEtab.jl package](https://github.com/sebapersson/PEtab.jl) implements the [PEtab format](https://petab.readthedocs.io/en/latest/) for parameter fitting within systems biology. It implements both methods for creating a cost function (determening a parameter set's ability to fit a model to data), and for minimizing said cost function. This format covers most cases of fitting deterministic (ODE) models to data. For other cases, an alternative approach should be used. This page describes how to use PEtab for by Catalyst defined `ReactionSystem`, with more extesnive documentation being avaiable at https://sebapersson.github.io/PEtab.jl/stable/.

## Introductory example
Let us consider a simple catalysis network, where an enzyme (`E`) turns a substrate (`S`) into a product (`P`):
```petab1
using Catalyst

rn = @reaction_network begin
    kB, S + E --> SE
    kD, SE --> S + E
    kP, SE --> P + E
end
```
From some known initial condition, and a parameter set (which we wish to fit) we generate synthetic data (on which we will demonstrate the fitting process).
```petab1
# Initial conditions an parameters.
u0 = [:S => 1.0, :E => 1.0, :SE => 1.0, :P => 0.0]
p_true = [:kB => 1.0, :kD => 0.1, :kP => 0.1]

# Simulate data.
using OrdinaryDiffEq, Random
sol = solve(ODEProblem(rn, u0, (0.0, 10.0), p_true); saveat=0.1)
data = 0.8 .+ 0.4*rand(10), sol[10:10:end, 3]

# Plots true and simualted data.
using Plots
plot(sol; idxs=3, label="True data")
plot(sol.t[10:10:end], data; label="Meassured data", seriestype=:scatter)
```

Generally, PEtab takes five different inputs to define an optimization problem (the minimum of which generates a fitted paraemter set):
1. **Model**: The model for which we want to fit parameters. This is specified as a by Catalyst defined `ReactionSystem`.
2. **Observables**: A list of possible observables that can be meassured. This is specified as a set of `PEtabObservable`s.
3. **Estimation parameters**: The parameters which you want to fit to data. This is specified as a set of `PEtabParameter`s. 
4. **Experimental conditions**: List of possible experimental conditions (each typically corresponding to a single experiment). This is specified as a set of `Dict`s (each setting the paraemter values used in a specific experiment).
5. **Measurements**: The meassurments to which the model is fitted. This is specified as a  `DataFrame`. 

In our example, we only have a single possible observable, the `P` species. We create a corresponding `PEtabObservable` which we then store in a dictionary with a corresponding reference tag:
```petab1
obs_P = PEtabObservable(P, 0.1)
observables = Dict("obs_P", obs_P)
nothing
```

For each of the three parameters we wish to fit (`kB`, `kD`, and `kP`). We create a single `PEtabParameter` for. Next, we store all in a single vector.
```petab1
par_kB = PEtabParameters(:kB , )
par_kD = PEtabParameters(:kD , )
par_kP = PEtabParameters(:kP , )
parameters = [par_kB, par_kD, par_kP]
nothing
```

The experimental conditions define parameters that varry between different experiments. Each condition is a dictionary with a list of parameter values, these are then strored in another dictionary, each given a tag. Here, since only a single experiment is performed, we do not need to define any conditions, and stores an empty dictionary:
```petab1
exp_1 = Dict()
experiments = Dict("exp_1" => exp_1)
nothing
```

Finally, we need to declare our meassurments. Each meassurements combines a value, a time point, an experiment, and an observable. These are then stored in a `DataFrame`. In our case, the same epxeirment and observable are used in all cases.
```petab1
measurements = DataFrame(simulation_id=["exp_1"], obs_id=["obs_P"], time=sol.t[10:10:end], measurement=data)
nothing
```

All this information can now be combined into a single `PEtabModel`. This one also takes the known initial conditions and parameters:
```petab1
petab_model = PEtabModel(rn, simulation_conditions, observables, measurements, parameters, state_map=u0, parameter_map=p)
nothing
```

Given all this information we can now fit parameters to data. We first create a `PEtabODEProblem` (combining information about the model, with how we want to simualte it).
```petab1
petab_problem = PEtabODEProblem(petab_model)
```

Parameters can be fitted using a selected Optimization method (PEtab provides a few alternatives):
```petab1
using Optim
petab_solution = calibrate_model(petab_problem, rand(4), IPNewton(), options=Optim.Options(iterations = 200))
```


## Multiple experiments


## Unknown initial conditions


## Additional arameter fitting options


## Citations
If you use this functionality in your reasearch, please cite the following papers to support the authors of the PEtab.jl package (currently there is no article associated with this paper) and the PEtab standard:

```
@misc{2023Petabljl,
  author       = {Ognissanti, Damiano and Arutjunjan, Rafael and Persson, Sebastian and Hasselgren, Viktor.
    
    Isaacson, S. A. and Ilin, V. and Rackauckas, C. V.},
  title        = {{2023Petabljl.jl}},
  howpublished = {\url{https://github.com/sebapersson/PEtab.jl/}},
  year         = {2023}
}
@article{SchmiesterSch2021,
  author    = {Schmiester, Leonard AND Schälte, Yannik AND Bergmann, Frank T. AND Camba, Tacio AND Dudkin, Erika AND Egert, Janine AND Fröhlich, Fabian AND Fuhrmann, Lara AND Hauber, Adrian L. AND Kemmer, Svenja AND Lakrisenko, Polina AND Loos, Carolin AND Merkt, Simon AND Müller, Wolfgang AND Pathirana, Dilan AND Raimúndez, Elba AND Refisch, Lukas AND Rosenblatt, Marcus AND Stapor, Paul L. AND Städter, Philipp AND Wang, Dantong AND Wieland, Franz-Georg AND Banga, Julio R. AND Timmer, Jens AND Villaverde, Alejandro F. AND Sahle, Sven AND Kreutz, Clemens AND Hasenauer, Jan AND Weindl, Daniel},
  journal   = {PLOS Computational Biology},
  title     = {PEtab—Interoperable specification of parameter estimation problems in systems biology},
  year      = {2021},
  month     = {01},
  number    = {1},
  pages     = {1-10},
  volume    = {17},
  doi       = {10.1371/journal.pcbi.1008646},
  publisher = {Public Library of Science},
  url       = {https://doi.org/10.1371/journal.pcbi.1008646},
}
```
## References
[^1]: [Schmiester, L et al. *PEtab—Interoperable specification of parameter estimation problems in systems biology*, PLOS Computational Biology (2021).](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008646)