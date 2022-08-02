# Bifurcation Diagrams
Bifurcation diagrams can be produced from Catalyst generated models through the use of the [BifurcationKit.jl](https://github.com/rveltz/BifurcationKit.jl/) package. This tutorial gives a simple example of how to create such a bifurcation diagram. 

First, we declare our model. For our example we will use a bistable switch, but which also contains a Hopf bifurcation.
```@example ex1
using Catalyst
rn = @reaction_network begin
    (v0+v*(S*X)^n/((S*X)^n+(D*A)^n+K^n),d), ∅ ↔ X
    (τ*X,τ), ∅ ↔ A
end S D τ v0 v K n d
```
Next, we specify the system parameters for which we wish to plot the bifurcation diagram. We also set the parameter we wish to vary in our bifurcation diagram, as well as over which interval. FInally, we set which variable we wish to plot the steady state values of in the diagram.
```@example ex1
p = [:S=>1., :D=>9., :τ=>1000., :v0=>0.01, :v=>2., :K=>20., :n=>3, :d=>0.05]
bif_par = :S      
p_span = (0.1,20.)   
plot_var = :X   
```
When creating a bifurcation diagram, we typically start in some point in parameter-phase space. For parameter space, we will simply select the beginning of the interval over which we wish to computer the bifurcation diagram. We thus create a modified parameter set where this is teh case. For this parameter set, we make a guess of an initial fixed point. While a good estimate could be provided through e.g. a simulation, the guess do not need to be very exact.
```@example ex1
p_bstart = setindex!(copy(p),bif_par => p_span[1],findfirst(first.(p).==bif_par))  
u0 = [:X=>1.0, :A=>1.0]
```
Finally, we extract the ode function and jacobian of our model in a form that BifurcationKit can read.
```@example ex1
oprob = ODEProblem(rn,u0,(0.0,0.0),p_bstart;jac=true)
F = (u,p) -> oprob.f(u,p,0)      
J = (u,p) -> oprob.f.jac(u,p,0)
```

Now, we fetch the required packages to create the bifurcation diagram. We also  bundle the information we have compiled so far into a "`BifurcationProblem`".
```@example ex1
using BifurcationKit, Plots, LinearAlgebra, Setfield
bprob = BifurcationProblem(F, oprob.u0, oprob.p, (@lens _[findfirst(first.(p).==bif_par)]); recordFromSolution = (x, p) -> x[findfirst(first.(u0).==plot_var)], J=J)
```
Next, we need to specify the input options for the pseudo-arclength continuation method which produces the diagram..
```@example ex1
bopts = ContinuationPar(dsmax = 0.05,        # Maximum arclength value of the pseudo-arc length continuation method.
                        dsmin = 1e-4,        # Minimum arclength value of the pseudo-arc length continuation method.
                        ds=0.001,            # Initial arclength value of the pseudo-arc length continuation method (should be positive).
                        maxSteps = 100000,   # The maximum number of steps.
                        pMin = p_span[1],    # Minimum p-vale (if hit, the method stops).
                        pMax = p_span[2],    # Maximum p-vale (if hit, the method stops).
                        detectBifurcation=3) # Value in {0,1,2,3} determining to what extent bifurcation points are detected (0 means nothing is done, 3 both them and there localisation are detected).
```
With all this done, we can compute the bifurcation diagram:
```@example ex1
bf = bifurcationdiagram(bprob, PALC(),2, (args...) -> bopts)
```
Finally, we can plot it:
```@example ex1
plot(bf)
```
![bifurcation_diagram1](../assets/bifurcation_diagram.svg)

Here, Hopf bifurcation is marked with a red dot and fold bifurcation with blue dots. The region with a thinner linewidth corresponds to unstable steady states. 

This tutorial demonstrates how to make a simple bifurcation diagram where all branches are connected. However, BifurcationKit.jl is a very powerful package capable of a lot more. For more details, please see that package's documentation: [BifurcationKit.jl](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/). 

