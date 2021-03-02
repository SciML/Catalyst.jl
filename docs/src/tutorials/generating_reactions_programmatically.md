## Population balance equations of the Smoluchowski coagulation model
This tutorial shows how to programmatically construct a `ReactionSystem` corresponding to the chemistry underlying the [Smoluchowski coagulation model](https://en.wikipedia.org/wiki/Smoluchowski_coagulation_equation) using [ModelingToolkit](https://mtk.sciml.ai/stable/)/[Catalyst](https://catalyst.sciml.ai/dev/). A jump process version of the model is then constructed from the `ReactionSystem`, and compared to the model's analytical solution obtained by the [method of Scott](https://journals.ametsoc.org/view/journals/atsc/25/1/1520-0469_1968_025_0054_asocdc_2_0_co_2.xml) (see also mentioned in reference [3](https://doi.org/10.1006/jcph.2002.7017).

The Smoluchowski coagulation equation describes a system of reactions in which monomers may collide to form dimers, monomers and dimers may collide to form trimers, and so on. This models a variety of chemical/physical processes, including polymerization and flocculation.

We begin by importing some necessary packages.
```julia
using ModelingToolkit, LinearAlgebra
using DiffEqBase, DiffEqJump, OrdinaryDiffEq
using LoopVectorization, Plots
using BenchmarkTools
using SpecialFunctions
```
Suppose the maximum cluster size is `N`. Lets initialize the system with some initial concentration `Nₒ`, initial number of monomers `uₒ` in the system. Since its a bimolecular chain of Reaction system(`nr` number of reactions), the bulk volume `V` of the system in which these binary collisions occur is important in the calculation of rate laws.
  
```julia
## Parameter
N = 10                       # maximum clusters size
Vₒ = (4π/3)*(10e-06*100)^3   # volume of a monomers in cm³
Nₒ = 1e-06/Vₒ                # initial conc. = (No. of init. monomers) / Volume of the bulk system
uₒ = 10000                   # No. of monomers initially
V = uₒ/Nₒ                    # Bulk volume of system in cm³

integ(x) = Int(floor(x));
n        = integ(N/2);
nr       = N%2 == 0 ? (n*(n + 1) - n) : (n*(n + 1)); # No. of forward reactions
```
Check the figure on [Smoluchowski coagulation equation](https://en.wikipedia.org/wiki/Smoluchowski_coagulation_equation) page, the `pair` of reactants that collide can be easily generated for `N` cluster size particles in the system. We also initialise the volumes of these colliding clusters as `volᵢ` and `volⱼ` for the reactants
  
```julia
## pairs of reactants
pair = [];
for i = 2:N
    push!(pair,[1:integ(i/2)  i .- (1:integ(i/2))])
end
pair = vcat(pair...);
vᵢ = @view pair[:,1];  # Reactant 1 index
vⱼ = @view pair[:,2];  # Reactant 2 index
volᵢ = Vₒ*vᵢ;    # cm⁻³
volⱼ = Vₒ*vⱼ;    # cm⁻³
sum_vᵢvⱼ = @. vᵢ + vⱼ;  # Product index
```
  - **4.)**  Specifying rate(kernel) at which reactants collide to form product. For simplicity we have used additive kernel, multiplicative kernel and constant kernel. The constants(`B`,`b` and `C`) used are adopted from the Scotts paper [2](https://journals.ametsoc.org/view/journals/atsc/25/1/1520-0469_1968_025_0054_asocdc_2_0_co_2.xml)
```julia
i = parse(Int, input("Enter 1 for additive kernel,
                2 for Multiplicative, 3 for constant"))
if i==1
    B = 1.53e03;    # s⁻¹
    kₛ = @. B*(volᵢ + volⱼ)/V;    # dividing by volume as its a bi-molecular reaction chain
elseif i==2
    b = 3.8e11;     #  cm⁻³ s⁻¹
    kₛ = @. b*(volᵢ*volⱼ)/V;
else
    C = 1.84e-04;   # cm³ s⁻¹
    kₛ = @. C/V;
end
```
  - **5.)**  Lets write-off the rates in `pₘₐₚ` as Pairs and initial condition with only monomers present initially in `u₀map` that we  will use in creating JumpSystems with massaction.
```julia
## Writing-off the parameter in Pairs in Sequence
@variables k[1:nr];   pₘₐₚ = Pair.(k, kₛ);
@parameters t;        @variables X[collect(1:N)](t);
if i == 1
    tspan = (0. ,2000.)   # time-span
elseif i == 2
    tspan = (0. ,3000.)
else
    tspan = (0. ,350.)
end
u₀ = zeros(Int64, N);   u₀[1] = uₒ;   # initial condition of monomers
u₀map = Pair.(X, u₀); # population of other polymers in zeros
```
  - **6.)** Here we generate the reactions programmatically. Push the Reactions(into empty reaction_network) as shown in figure at [here](https://en.wikipedia.org/wiki/Smoluchowski_coagulation_equation). When `vᵢ[n] == vⱼ[n]` ,we use rate as `2*k[n]` ,as coagulation kernel is related to deterministic rate-law in form, (to followed from a paper by [Laurenzi et.al](https://www.sciencedirect.com/science/article/pii/S0021999102970178))
          coagulation kernel = (1 + δᵢⱼ)*deterministic rate
                
```julia
rx = [];              # empty-reaction vector
reactant_stoich = Array{Array{Pair{Int64,Int64},1},1}(undef,nr);
net_stoich = Array{Array{Pair{Int64,Int64},1},1}(undef,nr);
@time for n = 1:nr
    if (vᵢ[n] == vⱼ[n])    # checking the reactants
        push!(rx, Reaction(2*k[n], [ X[vᵢ[n]] ] ,[ X[sum_vᵢvⱼ[n]] ] ,[2],[1]));

        reactant_stoich[n] = [vᵢ[n] => 2];
        net_stoich[n] = [vᵢ[n] => -2, sum_vᵢvⱼ[n] => 1];
    else
        push!(rx, Reaction(k[n], [ X[vᵢ[n]] , X[vⱼ[n]] ] ,[ X[sum_vᵢvⱼ[n]] ],
                                [1, 1],[1]));

        reactant_stoich[n] = [vᵢ[n] => 1 , vⱼ[n] => 1];
        net_stoich[n] = [vᵢ[n] => -1 , vⱼ[n] => -1 , sum_vᵢvⱼ[n] => 1];
    end
end
rs = ReactionSystem(rx, t, X, k);
```
  - **7.)**  Convert the reactionSystem into a JumpSystems and solve it using standard Jump solvers such as Gillespie process. For details, take a look at [DifferentialEquations](https://diffeq.sciml.ai/stable/) documentation 
```julia
## solving the system
jumpsys = convert(JumpSystem, rs; combinatoric_ratelaws = true);
dprob = DiscreteProblem(jumpsys, u₀map, tspan, pₘₐₚ; parallel = true);
alg = Direct();
stepper = SSAStepper();
mass_act_jump = MassActionJump(kₛ ,reactant_stoich, net_stoich);
jprob = @btime JumpProblem(dprob, alg ,mass_act_jump ,save_positions=(false,false));
jsol = @btime solve(jprob, stepper, saveat = 1.);
```
  - **8.)**  Lets check the results for only first three polymers/cluster sizes. The result is compared with analytical solution obtained for this system with additive, multiplicative and constant kernels(rate at which reactants collide)
```julia
## Results for first three polymers...i.e. monomers, dimers and trimers
v_res = [1;2;3]

## comparsion with analytical solution
if i == 1
    ϕ = @. 1 - exp(-B*Nₒ*Vₒ*jsol.t);        # normalised "time"
    sol = zeros(length(v_res), length(ϕ))
    for j in v_res
        sol[j,:] = @. Nₒ*(1 - ϕ)*(((j*ϕ)^(j-1))/gamma(j+1))*exp(-j*ϕ);
    end
elseif i == 2
    ϕ = @. (b*Nₒ*Vₒ*Vₒ*jsol.t);             # normalised "time"
    sol = zeros(length(v_res), length(ϕ))
    for j in v_res
        sol[j,:] = @. Nₒ*(((j*ϕ)^(j-1))/(j*gamma(j+1)))*exp(-j*ϕ);
    end
else
    ϕ = @. (C*Nₒ*jsol.t);                   # normalised "time"
    sol = zeros(length(v_res), length(ϕ))
    for j in v_res
        sol[j,:] = @. 4Nₒ*((ϕ^(j-1))/((ϕ + 2)^(j+1)));
    end
end
# plotting normalised concentration vs analytical solution
scatter(ϕ, jsol[1,:]/uₒ, lw = 2, xlabel = "Time (sec)",label = string("X",1))
plot!(ϕ, sol[1,:]/Nₒ, lw = 2, line = (:dot, 4) ,
        label = string("Analytical sol", "--X",1))

scatter!(ϕ, jsol[2,:]/uₒ, lw = 2, xlabel = "Time (sec)",label = string("X",2))
plot!(ϕ, sol[2,:]/Nₒ, lw = 2, line = (:dot, 4) ,
        label = string("Analytical sol", "--X",2))

scatter!(ϕ, jsol[3,:]/uₒ, lw = 2, xlabel = "Time (sec)",label = string("X",3))
plot!(ϕ, sol[3,:]/Nₒ, lw = 2, line = (:dot, 4) ,ylabel = "N(i,t)/N(1,0)"
        ,label = string("Analytical sol", "--X",3))

```
- **9.)** Example plot for **additive kernel** 
               ![](https://github.com/yewalenikhil65/Catalyst.jl/blob/master/docs/src/assets/additive_kernel.png)

## Sources
1. https://en.wikipedia.org/wiki/Smoluchowski_coagulation_equation
2. Scott, W. T. (1968). Analytic Studies of Cloud Droplet Coalescence I, Journal of Atmospheric Sciences, 25(1), 54-65. Retrieved Feb 18, 2021, from https://journals.ametsoc.org/view/journals/atsc/25/1/1520-0469_1968_025_0054_asocdc_2_0_co_2.xml
3. Ian J. Laurenzi, John D. Bartels, Scott L. Diamond, A General Algorithm for Exact Simulation of Multicomponent Aggregation Processes, Journal of Computational Physics, Volume 177, Issue 2, 2002, Pages 418-449, ISSN 0021-9991, https://doi.org/10.1006/jcph.2002.7017.
