# Catalyst.jl

[![Join the chat at https://gitter.im/JuliaDiffEq/Lobby](https://badges.gitter.im/JuliaDiffEq/Lobby.svg)](https://gitter.im/JuliaDiffEq/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![Build Status](https://travis-ci.com/SciML/Catalyst.jl.svg?branch=master)](https://travis-ci.com/SciML/Catalyst.jl)
[![Coverage Status](https://coveralls.io/repos/github/SciML/Catalyst.jl/badge.svg?branch=master)](https://coveralls.io/github/SciML/Catalyst.jl?branch=master)
[![codecov.io](https://codecov.io/gh/SciML/Catalyst.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/SciML/Catalyst.jl)

<!--- [![Build status](https://ci.appveyor.com/api/projects/status/github/SciML/Catalyst.jl?branch=master&svg=true)](https://ci.appveyor.com/project/ChrisRackauckas/Catalyst-jl/branch/master) --->
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://docs.sciml.ai/stable/models/catalyst/)
[![API Stable](https://img.shields.io/badge/API-stable-blue.svg)](https://docs.sciml.ai/stable/apis/catalyst_api/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://docs.sciml.ai/latest/models/catalyst/)
[![API Dev](https://img.shields.io/badge/API-dev-blue.svg)](https://docs.sciml.ai/latest/apis/catalyst_api/)

**Note for pre-version 5 users**: *Version 5 is a breaking release, with the DSL
now generating `ModelingToolkit.ReactionSystem`s and DiffEqBiological being
renamed to Catalyst.  As such, the `@reaction_network` macro no longer allows
the generation of custom types. Please see the updated documentation to
understand changes to the API and functionality. In particular, the earlier
bifurcation functionality has not yet been updated to the new system. If you
rely on this functionality please do not update at this time, or consider using
[BifurcationKit.jl](https://github.com/rveltz/BifurcationKit.jl).*

Catalyst.jl provides a domain specific language (DSL) for defining
chemical reaction networks in Julia, generating
[ModelingToolkit](https://github.com/SciML/ModelingToolkit.jl)
`ReactionSystem`s. These can be converted to ModelingToolkit-based systems of
ODEs, SDEs, jump processes and more. This allows for the easy generation and
solution of mass action ODE models, Chemical Langevin SDE models, stochastic
chemical kinetics jump process models, and more. The generated models can then
be used with solvers throughout the broader [SciML](https://sciml.ai) ecosystem,
including higher level SciML packages (e.g. for sensitivity analysis, parameter
estimation, machine learning applications, etc).

Here is a simple example of generating and solving an SIR ODE model:
```julia
using Catalyst, DiffEqBase, OrdinaryDiffEq
rn = @reaction_network begin
    k1, S + I --> R
    k2, I --> R
end k1 k2
p     = [.1/1000, .01]           # [k1,k2]
tspan = (0.0,250.0) 
u0    = [999.0,1.0,0.0]          # [S,I,R] at t=0
op    = ODEProblem(rn, u0, t, p)
sol   = solve(op, Tsit5())       # use Tsit5 ODE solver
```

Here we give a brief introduction to using the Catalyst package, with a
focus on how to define reaction networks, basic properties of the generated
`ReactionSystem`s, and a minimal example showing how to create and solve ODE,
SDE and jump models.

More detailed documentation includes:
<!-- * Several Catalyst tutorials are available as part of the
  [DiffEqTutorials Modeling
  Examples](https://github.com/JuliaDiffEq/DiffEqTutorials.jl). Both html and
  interactive IJulia notebook versions are provided there. These include
  * An introductory tutorial showing how to specify and solve both ODE and
  stochastic versions of the
  [repressilator](https://en.wikipedia.org/wiki/Repressilator).
  * A tutorial exploring the Catalyst API for querying network
    properties, which also illustrates how to programmatically and incrementally construct and
    solve a network model using the API.
  * A tutorial showing how to use the wrapped [HomotopyContinuation.jl](https://www.juliahomotopycontinuation.org/)
    functionality to find steady-states and make bifurcation plots. -->
* A full introduction to the DSL, with information on the generated
  `ReactionSystem`s and their conversion to ODE/SDE/jump process models is
  available in the [DifferentialEquations.jl Catalyst
  documentation](http://docs.sciml.ai/latest/models/catalyst).
* API documentation showing how to retrieve network information from a
  generated `reaction_network` is also available
  [here](http://docs.sciml.ai/latest/apis/catalyst_api).

## The Reaction DSL

The `@reaction_network` DSL allows for the definition of reaction networks using
a simple format. Its input is a set of chemical reactions, from which it
generates a `ModelingToolkit.ReactionSystem`. The latter can be converted to a
ModelingToolkit `ODESystem`, `SDESystem`, `JumpSystem` and more, which can
be used as input to building problems for use in corresponding solvers.

The basic syntax is
```julia
rn = @reaction_network begin
  2.0, X + Y --> XY               
  1.0, XY --> Z            
end
```
where each line corresponds to a chemical reaction. 

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
specify some reactions through a programmatic
[API](http://docs.sciml.ai/dev/catalyst_api/index.html#Functions-to-extend-a-Network).
In this case one can directly construct ModelingToolkit `Reaction`s and
`ReactionSystem`s, or construct an empty reaction network using
```julia
rn = make_empty_network()
```
or
```julia 
rn = @reaction_network
```
This can then be filled in using the API `addspecies!`, `addparam!` and
 `addreaction!` modifier functions. `merge` and `merge!` are also supported for
 constructing composed networks.


## Catalyst API for Querying Network Information

A variety of network information is generated by the `reaction_network` macro,
and can then be retrieved using the [Catalyst
API](http://docs.sciml.ai/dev/apis/catalyst_api/), the
`ModelingToolkit.AbstractSystem` API, or fields within the generated
`ReactionSystem` and `Reaction`s. The Catalyst API includes

* Basic properties of the `ReactionSystem`
  ```julia
    species(rn)              # vector of ModelingToolkit.Variables
    params(rn)               # vector of ModelingToolkit.Variables
    reactions(rn)            # vector of all Reactions within the system
  ```
* Integer ids of species and parameters
  ```julia
    speciesmap(rn)
    paramsmap(rn)
  ```
* Reaction information. Given a `Reaction` within a `ReactionSystem`, `rx`, we
  can extract basic information such as substrates and stoichiometries:
  ```julia
    rx = reactions(rn)[1]       # first reaction in rn
    rx.rate                     # ModelingToolki Operation representing the rate expression
    rx.substrates               # ModelingToolkit variables for each substrate
    rx.substoich                # vector of stoichiometric coefs for substrates
    rx.products
    rx.prodstoich
    rx.netstoich                # vector of Pairs of the form [X=>1]
  ```
* Properties of reactions:
  ```julia
  ismassaction(rx)                
  dependents(rx)                # species the full rate law for rx depends on
  oderatelaw(rx)
  jumpratelaw(rx)
  ```
* Functions and macros for extending and combining networks
  ```julia
  @parameter k                                  # create a ModelingToolkit parameter
  @variable S                                   # create a ModelingToolkit variable
  addspecies!(rn, S)
  addparam!(rn, k)
  addreaction!(rn, Reaction(k, [S], nothing))   # add S --> 0 into rn
  @add_reactions rn begin                       # add S --> 0 into rn
    k2, S --> 2S
    k3, 2S --> S
  end k2 k3
  merge!(rn, rn2)                               # merge rn2 into rn1
  rn3 = merge(rn, rn2)
  ```
In addition, leveraging `Modelingtoolkit` we can
* Convert `ReactionSystem`s to ODE, SDE, jump process models and more.
* Take advantage of Modelingtoolkit features for calculating Jacobians,
  sparse Jacobians, multithreaded ODE derivatives, encoding jumps into optimal types, etc.
* Obtain Julia `Expr`s corresponding to the functions determining the deterministic and
  stochastic terms within resulting ODE, SDE or jump models.
* Use [`Latexify`](https://github.com/korsbo/Latexify.jl) to generate LaTeX
  expressions corresponding to a given generated mathematical model.
* Obtain dependency graphs showing the relationship between reactions, species and parameters.

and more.

## Simulating ODE, Steady-State, SDE and Jump Problems

Once a reaction network has been created it can be passed as input to either one
of the `ODEProblem`, `SteadyStateProblem`, `SDEProblem` or `JumpProblem`
constructors.
```julia
  probODE = ODEProblem(rn, args...; kwargs...)      
  probSS = SteadyStateProblem(rn, args...; kwargs...)
  probSDE = SDEProblem(rn, args...; kwargs...)
  probJump = JumpProblem(rn, prob, Direct())
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
p = (1.0,2.0,50.)
tspan = (0.,4.)
u0 = [5.]

# solve ODEs
oprob = ODEProblem(rs, u0, tspan, p)
osol  = solve(oprob, Tsit5())

# solve for Steady-States
ssprob = SteadyStateProblem(rs, u0, p)
sssol  = solve(ssprob, SSRootfind())

# solve Chemical Langevin SDEs
sprob = SDEProblem(rs, u0, tspan, p)
ssol  = solve(sprob, EM(), dt=.01)

# solve JumpProblem using Gillespie's Direct Method
u0 = [5]
dprob = DiscreteProblem(rs, u0, tspan, p)
jprob = JumpProblem(rs, dprob, Direct())
jsol = solve(jprob, SSAStepper())
```

Finer control over the generation of such models can be obtained by explicitly
converting the `ReactionSystem` to another system type and then constructing a
problem. For example
```julia
osys = convert(ODESystem, rs)     
oprob = ODEProblem(osys, Pair.(species(rs),u0), tspan, Pair.(params(rs),p))
osol  = solve(oprob, Tsit5())
```
would be equivalent to solving via the direct `ODEProblem` construction above.
Using an `ODESystem` intermediate allows the possibility to modify the system
with further terms that are difficult to encode as a chemical reaction. For
example, suppose we wish to add a forcing term, `10*sin(10*t)`, to the ODE for
`dX/dt`. We can do so as:
```julia
dXdteq = equations(osys)[1]           
t      = independent_variable(osys)()    
dXdteq = Equation(dXdteq.lhs, dXdteq.rhs + 10*sin(10*t))   
osys2  = ODESystem([dXdteq], t, states(osys), parameters(osys))
oprob  = ODEProblem(osys2, Pair.(species(rs),u0), tspan, Pair.(params(rs),p))
osol   = solve(oprob, Tsit5())
```

<!-- ## Importing Predefined Networks 
[ReactionNetworkImporters.jl](https://github.com/isaacsas/ReactionNetworkImporters.jl)
can load several different types of predefined networks directly into
`ReactionSystem`s. These include 
  * A subset of BioNetGen .net files that can be generated from a BioNetGen language file (.bngl). (.net files can be generated using the `generate_network` command within BioNetGen.) 
  * Reaction networks specified by dense or sparse matrices encoding the stoichiometry of substrates and products within each reaction.
  * Networks defined by the basic file format used by the [RSSA](https://www.cosbi.eu/research/prototypes/rssa) group at COSBI in their [model collection](https://www.cosbi.eu/prototypes/jLiexDeBIgFV4zxwnKiW97oc4BjTtIoRGajqdUz4.zip). -->

<!-- ## Finding steady states
The steady states of a reaction network can be found using homotopy continuation (as implemented by [HomotopyContinuation.jl](https://github.com/isaacsas/ReactionNetworkImporters.jl)). This method is limited to polynomial systems, which includes reaction network not containing non-polynomial rates in the reaction rates (such as logarithms and non integer exponents). *Note, both the steady-state and the bifurcation diagram functionality only fully support Julia 1.1 and greater.*

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
@add_constraint rn X+Y=2.
steady_states(rn,params)
```
The `@add_constraint` macro may contain parameters, as long as these are declared in the network.
```julia
rn = @reaction_network begin 
  (k1,k2), X ↔ Y              
end k1 k2 C_tot
params = [2.,1.,2.]
@add_constraint rn X+Y=C_tot
steady_states(rn,params)
```
The `@add_constraints` macro can be used to add several constraints at the same time.
```julia
rn = @reaction_network begin 
  (k1,k2), X ↔ Y        
  (k3,k4), V ↔ W              
end k1 k2 k3 k4
params = [2.,1.,1.,2.]
@add_constraints rn begin
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
bif_grid_2d = bifurcation_grid_2d(rn, params, :p, 1.:5., :d, 2.:10.)
```
A bifurcation diagram grid first varies a single variable over a discrete grid of values. Then, for each such value, in varies a second variable over a continuous interval to create a bifurcation grid.
```julia
bif_grid_dia = bifurcation_grid_diagram(rn, params, :p, 1.:5., :d, (2.,10.))
```
All of these can be plotted.
```julia
plot(bif_grid)
plot(bif_grid_2d)
plot(bif_grid_dia)
```
``` -->
