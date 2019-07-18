# DiffEqBiological.jl

[![Join the chat at https://gitter.im/JuliaDiffEq/Lobby](https://badges.gitter.im/JuliaDiffEq/Lobby.svg)](https://gitter.im/JuliaDiffEq/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![Build Status](https://travis-ci.org/JuliaDiffEq/DiffEqBiological.jl.svg?branch=master)](https://travis-ci.org/JuliaDiffEq/DiffEqBiological.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/github/JuliaDiffEq/DiffEqBiological.jl?branch=master&svg=true)](https://ci.appveyor.com/project/ChrisRackauckas/diffeqbiological-jl/branch/master)
[![Coverage Status](https://coveralls.io/repos/github/JuliaDiffEq/DiffEqBiological.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaDiffEq/DiffEqBiological.jl?branch=master)
[![codecov.io](https://codecov.io/gh/JuliaDiffEq/DiffEqBiological.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaDiffEq/DiffEqBiological.jl)

DiffEqBiological.jl provides a domain specific language (DSL) for defining
chemical reaction networks in Julia. It interfaces with the broader
[DifferentialEquations.jl](http://juliadiffeq.org) infrastructure to enable the
easy generation and solution of corresponding mass action ODE models, Chemical
Langevin SDE models, and stochastic chemical kinetics jump models. These
generated models can also be used in higher level DifferentialEquations.jl
packages (e.g. for sensitivity analysis, parameter estimation, etc).

Here we give a brief introduction to using the DiffEqBiological package, with
a focus on how to define reaction networks, and a minimal example showing how to
create and solve ODE, SDE and jump models.

More detailed documentation is available from:
* Several DiffEqBiological tutorials are available as part of the
  [DiffEqTutorials Modeling
  Examples](https://github.com/JuliaDiffEq/DiffEqTutorials.jl). Both html and
  interactive IJulia notebook versions are provided there. These include
  * An introductory tutorial showing how to specify and solve both ODE and
  stochastic versions of the
  [repressilator](https://en.wikipedia.org/wiki/Repressilator).
  * A tutorial exploring the DiffEqBiological API for querying network
    properties, which also illustrates how to programmatically construct and
    solve a network model using the API.
* Full documentation of the DSL syntax, with information on the generated rate
  functions and models is available in the [DifferentialEquations.jl Chemical
  Reaction Models
  documentation](http://docs.juliadiffeq.org/latest/models/biological.html).
* API documentation showing how to retrieve network information from a
  generated `reaction_network` is available
  [here](http://docs.juliadiffeq.org/latest/apis/diffeqbio.html).

## The Reaction DSL

The `@reaction_network` DSL allows for the definition of reaction networks using
a simple format. Its input is a set of chemical reactions, from which it
generates a reaction network object which can be used as input to `ODEProblem`,
`SteadyStateProblem`, `SDEProblem` and `JumpProblem` constructors.

The basic syntax is
```julia
rn = @reaction_network rType begin
  2.0, X + Y --> XY               
  1.0, XY --> Z            
end
```
where each line corresponds to a chemical reaction. The (optional) input `rType`
designates the type of this instance (all instances will inherit from the
abstract type `AbstractReactionNetwork`).

The DSL has many features:
* It supports many different arrow types, corresponding to different directions
  of reactions and different rate laws:
  ```julia
  rn = @reaction_network begin
    1.0, X + Y --> XY               
    1.0, X + Y → XY      
    1.0, XY ← X + Y      
    2.0, X + Y ↔ XY               
  end
  ```
* It allows multiple reactions to be defined simultaneously on one line. The
  following two networks are equivalent:
  ```julia
  rn1 = @reaction_network begin
    (1.0,2.0), (S1,S2) → P             
  end
  rn2 = @reaction_network begin
    1.0, S1 → P     
    2.0, S2 → P
  end
  ```
* It allows the use of symbols to represent reaction rate parameters, with their
  numeric values specified during problem construction. i.e., the previous
  example could be given by
  ```julia
  rn2 = @reaction_network begin
    k1, S1 → P     
    k2, S2 → P
  end k1 k2 
  ```
  with `k1` and `k2` corresponding to the reaction rates.
* Rate law functions can be pre-defined and used within the DSL:
  ```julia
  @reaction_func myHill(x) = 2.0*x^3/(x^3+1.5^3)
  rn = @reaction_network begin
    myHill(X), ∅ → X
  end
  ```
* Pre-made rate laws for Hill and Michaelis-Menten reactions are provided:
  ```julia
  rn1 = @reaction_network begin
    hill(X,v,K,n), ∅ → X
    mm(X,v,K), ∅ → X
  end v K n
  ```
* Simple rate law functions of the species populations can be used within the DSL:
  ```julia
  rn = @reaction_network begin
    2.0*X^2, 0 → X + Y
    gamma(Y)/5, X → ∅
    pi*X/Y, Y → ∅
  end
  ```

For sufficiently large and structured network models it can often be easier to
specify some reactions through a programmatic API. For this reason the
`@min_reaction_network` and `@empty_reaction_network` macros, along with the
corresponding `addspecies!`, `addparam!` and `addreaction!` modifier functions,
are provided in the
[API](http://docs.juliadiffeq.org/latest/apis/diffeqbio.html#Functions-to-Add-Species,-Parameters-and-Reactions-to-a-Network-1).


## DiffEqBiological API for Querying Network Information

A variety of network information is calculated by the `reaction_network` macro,
and can then be retrieved using the [DiffEqBiological
API](http://docs.juliadiffeq.org/latest/apis/diffeqbio.html). This includes

* Orderings of species and reactions
  ```julia
    speciesmap(rn)
    paramsmap(rn)
  ```
* Reaction stoichiometries
  ```julia
    substratestoich(rn, rxidx)
    productstoich(rn, rxidx)
    netstoich(rn, rxidx)
  ```
* Expressions corresponding to the functions determining the deterministic and
  stochastic terms within resulting ODE, SDE or jump models
  ```julia
    ode_exprs = odeexprs(rn)
    jacobian_exprs = jacobianexprs(rn)
    noise_expr = noiseexprs(rn)
    rate_exprs,affect_exprs = jumpexprs(rn)  
  ```
  These can be used to generate LaTeX expressions corresponding to the system
  using packages such as [`Latexify`](https://github.com/korsbo/Latexify.jl).
* Dependency graphs
  ```julia
    rxtospecies_depgraph(rn)
    speciestorx_depgraph(rn)
    rxtorx_depgraph(rn)
  ```
and more.

## Simulating ODE, Steady-State, SDE and Jump Problems

Once a reaction network has been created it can be passed as input to either one
of the `ODEProblem`, `SteadyStateProblem`, `SDEProblem` or `JumpProblem`
constructors.
```julia
  probODE = ODEProblem(rn, args...; kwargs...)      
  probSS = SteadyStateProblem(rn, args...; kwargs...)
  probSDE = SDEProblem(rn, args...; kwargs...)
  probJump = JumpProblem(prob, Direct(), rn)
```
The output problems may then be used as input to the solvers of
[DifferentialEquations.jl](http://juliadiffeq.org/). *Note*, the noise used by
the `SDEProblem` will correspond to the Chemical Langevin Equations. 

As an example, consider models for a simple birth-death process:
```julia
rs = @reaction_network begin
  c1, X --> 2X
  c2, X --> 0
  c3, 0 --> X
end c1 c2 c3
params = (1.0,2.0,50.)
tspan = (0.,4.)
u0 = [5.]

# solve ODEs
oprob = ODEProblem(rs, u0, tspan, params)
osol  = solve(oprob, Tsit5())

# solve for Steady-States
ssprob = SteadyStateProblem(rs, u0, params)
sssol  = solve(ssprob, SSRootfind())

# solve Chemical Langevin SDEs
sprob = SDEProblem(rs, u0, tspan, params)
ssol  = solve(sprob, EM(), dt=.01)

# solve JumpProblem using Gillespie's Direct Method
u0 = [5]
dprob = DiscreteProblem(rs, u0, tspan, params)
jprob = JumpProblem(dprob, Direct(), rs)
jsol = solve(jprob, SSAStepper())
```

## Importing Predefined Networks 
[ReactionNetworkImporters.jl](https://github.com/isaacsas/ReactionNetworkImporters.jl)
can load several different types of predefined networks into DiffEqBiological
`reaction_network`s. These include 
  * A subset of BioNetGen .net files that can be generated from a BioNetGen language file (.bngl). (.net files can be generated using the `generate_network` command within BioNetGen.) 
  * Reaction networks specified by dense or sparse matrices encoding the stoichiometry of substrates and products within each reaction.
  * Networks defined by the basic file format used by the [RSSA](https://www.cosbi.eu/research/prototypes/rssa) group at COSBI in their [model collection](https://www.cosbi.eu/prototypes/jLiexDeBIgFV4zxwnKiW97oc4BjTtIoRGajqdUz4.zip).

## Finding steady states
The steady states of a reaction network can be found using homotopy continuation (as implemented by [HomotopyContinuation.jl](https://github.com/isaacsas/ReactionNetworkImporters.jl)). This method is limited to polynomial systems, which includes reaction network not containing non-polynomial rates in the reaction rates (such as logarithms and non integer exponents).

The basic syntax is
```julia
rn = @reaction_network begin 
  (1.,2.), 0 ↔ X               
end
ss = steady_states(rn,params)
```
and with parameters
```julia
rn = @reaction_network begin 
  (p,d), 0 ↔ X               
end p d
params = [1., 2.]
ss = steady_states(rn,params)
```
the stability of a steady state (or a vector of several) can be determined by the `stability` function:
```julia
stability(ss,rn,params)
```

Here the `@reaction_network` creates a multivariate polynomial and stores in the `equilibrate_content` field in the reaction network structure. The `steady_state` method the inserts the corresponding parameter values and solves the polynomial system. The exception is if there exists a parameter as an exponent (typically `n` in a hill function). In this case the steady state polynomial can first be created in the `steady_state` method. If one plans to solve a polynomial a large number of times with the same value of `n`, then one can get a speed-up by first fixing that value using the `fix_parameters` function:
```julia
rn = @reaction_network begin 
  (hill(X,v,K,n),d), 0 ↔ X               
end v K n d
fix_parameters(rn,n=4)
for i = 1:10000
  params = [i, 2.5, 4, 0.1]    #The value of 'n' here doesn't really matter, however, the field must exist.
  ss = steady_states(rn,params)
```

Some networks may have an infinite set of steady states, and which one is interested in depends on the initial conditions. For these networks some additional information is required (typically some concentrations which sums to a fixed value). This information can be added through the `@add_constraint` macro:
```julia
rn = @reaction_network begin 
  (k1,k2), X ↔ Y              
end k1 k2
params = [2.,1.]
@add_constraint X+Y=2.
steady_states(rn,params)
```
The `@add_constraint` macro may contain parameters, as long as these are declared in the network.
```julia
rn = @reaction_network begin 
  (k1,k2), X ↔ Y              
end k1 k2 C_tot
params = [2.,1.,2.]
@add_constraint X+Y=C_tot
steady_states(rn,params)
```
The `@add_constraints` macro can be used to add several constraints at the same time.
```julia
rn = @reaction_network begin 
  (k1,k2), X ↔ Y        
  (k3,k4), V ↔ W              
end k1 k2 k3 k4
params = [2.,1.,1.,2.]
@add_constraints begin
  X + Y = 2.
  V + W = 4.
end
steady_states(rn,params)
```

## Making bifurcation diagram
For any system for which we can find steady states, we can also make bifurcation diagrams.
```julia
rn = @reaction_network begin 
  (p,d), 0 ↔ X               
end p d
params = [1.,2.] #The value of 'p' here doesn't really matter, however, the field must exist.
bif = bifurcations(rn, params, :p, (0.1,5.))
```
These can then be plotted.
```julia
plot(bif)
```
In the plot blue values correspond to stable steady states, red to unstable. Also, cyan correspond to stable steady states with imaginary eigen values and orange to unstable steady states with imaginary eigen values.

In addition to the normal bifurcation diagram (varying a single parameter over a continuous range) there are three more types available.

A bifurcation grid varies a single parameter over a set of discrete values
```julia
bif_grid = bifurcation_grid(rn, params, :p, 1.:5.)
```
A two dimensional bifurcation grid varies two different parameters over a grid of discrete values.
```julia
bif_grid_2d = bifurcation_grid_2d(rn, params, :p, 1.:5. :d, 2.:10.)
```
A bifurcation diagram grid first varies a single variable over a discrete grid of values. Then, for each such value, in varies a second variable over a continuous interval to create a bifurcation grid.
```julia
bif_grid_dia = bifurcation_grid_diagram(rn, params, :p, 1.:5. :d, (2.,10.))
```
All of these can be plotted.
```julia
plot(bif_grid)
plot(bif_grid_2d)
plot(bif_grid_dia)
```
```
