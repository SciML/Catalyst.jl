# Network Representations in Catalyst

In this tutorial we introduce several of the Catalyst API functions for network
analysis. A complete summary of the exported functions is given in the API
section
[`Network-Analysis-and-Representations`](https://catalyst.sciml.ai/dev/api/catalyst_api/#Network-Analysis-and-Representations).
We illustrate these functions on the repressilator model from previous tutorials.

## Network representation of the Repressilator `ReactionSystem`
We first load Catalyst and construct our model of the repressilator
```@example s1
using Catalyst
repressilator = @reaction_network Repressilator begin
       hillr(P₃,α,K,n), ∅ --> m₁
       hillr(P₁,α,K,n), ∅ --> m₂
       hillr(P₂,α,K,n), ∅ --> m₃
       (δ,γ), m₁ <--> ∅
       (δ,γ), m₂ <--> ∅
       (δ,γ), m₃ <--> ∅
       β, m₁ --> m₁ + P₁
       β, m₂ --> m₂ + P₂
       β, m₃ --> m₃ + P₃
       μ, P₁ --> ∅
       μ, P₂ --> ∅
       μ, P₃ --> ∅
end α K n δ γ β μ
reactions(repressilator)
```
In the [Using Catalyst](https://catalyst.sciml.ai/dev/tutorials/using_catalyst/)
tutorial we showed how the above network could be visualized as a
species-reaction graph. There species are represented by the nodes of the graph
and edges show the reactions in which a given species is a substrate or product.
```julia
g = Graph(repressilator)
```
![Repressilator solution](../assets/repressilator.svg)

We also showed in the [Using
Catalyst](https://catalyst.sciml.ai/dev/tutorials/using_catalyst/) tutorial that
the reaction rate equation ODE model for the repressilator is
```math
\begin{aligned}
\frac{dm_1(t)}{dt} =& \frac{\alpha K^{n}}{K^{n} + \left( {P_3}\left( t \right) \right)^{n}} - \delta {m_1}\left( t \right) + \gamma \\
\frac{dm_2(t)}{dt} =& \frac{\alpha K^{n}}{K^{n} + \left( {P_1}\left( t \right) \right)^{n}} - \delta {m_2}\left( t \right) + \gamma \\
\frac{dm_3(t)}{dt} =& \frac{\alpha K^{n}}{K^{n} + \left( {P_2}\left( t \right) \right)^{n}} - \delta {m_3}\left( t \right) + \gamma \\
\frac{dP_1(t)}{dt} =& \beta {m_1}\left( t \right) - \mu {P_1}\left( t \right) \\
\frac{dP_2(t)}{dt} =& \beta {m_2}\left( t \right) - \mu {P_2}\left( t \right) \\
\frac{dP_3(t)}{dt} =& \beta {m_3}\left( t \right) - \mu {P_3}\left( t \right)
\end{aligned}
```

## Matrix-Vector Reaction Rate Equation Representation
In general, reaction rate equation (RRE) ODE models for chemical reaction networks can
be represented as a first order system of ODEs in a compact matrix-vector notation. Suppose
we have a reaction network with ``K`` reactions and ``M`` species, labelled by the state vector
```math
\mathbf{x}(t) = \begin{pmatrix} x_1(t) \\ \vdots \\ x_M(t)) \end{pmatrix}.
```
For the repressilator, ``\mathbf{x}(t)`` is just
```@example s1
x = species(repressilator)
```
The RRE ODEs satisfied by $\mathbf{x}(t)$ are then
```math
\frac{d\mathbf{x}}{dt} = N \mathbf{\nu}(\mathbf{x}(t),t),
```
where ``N`` is a constant ``M`` by ``K`` matrix with ``N_{m k}`` the net
stoichiometric coefficient of species ``m`` in reaction ``k``.
``\mathbf{\nu}(\mathbf{x}(t),t)`` is the rate law vector, with
``\nu_k(\mathbf{x}(t),t)`` the rate law for the ``k``th reaction. For example,
for the first reaction of the repressilator above, the rate law is
```math
\nu_1(\mathbf{x}(t),t) = \frac{\alpha K^{n}}{K^{n} + \left( P_3(t) \right)^{n}}
```
We can calculate each of these in Catalyst via
```@example s1
N = netstoichmat(repressilator)
```
and by using the [`oderatelaw`](@ref) function
```@example s1
rxs = reactions(repressilator)
ν = oderatelaw.(rxs)
```
Note, as [`oderatelaw`](@ref) takes just one reaction as input we use
broadcasting to apply it to each element of `rxs`.

Let's check this really gives the same ODEs as Catalyst. Here is what Catalyst
generates by converting to an `ODESystem`
```@example s1
osys = convert(ODESystem, repressilator)

# for display purposes we just pull out the right side of the equations
odes = [eq.rhs for eq in equations(osys)]
```
whereas our matrix-vector representation gives
```@example s1
odes2 = N * ν
```
Let's check these are equal symbolically
```@example s1
isequal(odes, odes2)
```

## Reaction Complex Representation
We now introduce a further decomposition of the RRE ODEs, which has been used to
facilitate analysis of a variety of reaction network properties. Consider a simple
reaction system like
```@example s1
rn = @reaction_network begin
 k*A, 2*A + 3*B --> A + 2*C + D
 b, C + D --> 2*A + 3*B
end k b
show(stdout, MIME"text/plain"(), rn) # hide
```
We can think of the first reaction as converting the *reaction complex*,
``2A:2B`` to the complex ``A:2C:D`` with rate ``2*A``. Suppose we order our
species the same way as Catalyst does, i.e.
```math
\begin{pmatrix}
x_1(t)\\
x_2(t)\\
x_3(t)\\
x_4(t)
\end{pmatrix} =
\begin{pmatrix}
A(t)\\
B(t)\\
C(t)\\
D(t)
\end{pmatrix},
```
which should be the same as
```@example s1
species(rn)
```
We can describe a given reaction complex by the stoichiometric coefficients of
each species within the complex. For the reactions in `rn` these vectors would
be
```math
\begin{align*}
2A:3B = \begin{pmatrix}
2\\
3\\
0\\
0
\end{pmatrix}, &&
A:2C:D = \begin{pmatrix}
1\\
0\\
2\\
1
\end{pmatrix},
 &&
C:D = \begin{pmatrix}
0\\
0\\
1\\
1
\end{pmatrix}
\end{align*}
```
Catalyst can calculate these representations as the columns of the complex
stoichiometry matrix,
```@example s1
Z = complexstoichmat(rn)
```
If we have ``C`` complexes, ``Z`` is a ``M`` by ``C`` matrix with ``Z_{m c}``
giving the stoichiometric coefficient of species ``m`` within complex ``c``.

We can use this representation to provide another representation of the RRE
ODEs. The net stoichiometry matrix can be factored as ``N = Z B``, where ``B``
is called the incidence matrix of the reaction network,
```@example s1
B = incidencemat(rn)
```
Here ``B`` is a ``C`` by ``K`` matrix with ``B_{c k} = 1`` if complex `c`
appears as a product of reaction `k`, and ``B_{c k} = -1`` if complex `c` is a
substrate of reaction `k`.

Using our decomposition of ``N``, the RRE ODEs become
```math
\frac{dx}{dt} = Z B \mathbf{\nu}(\mathbf{x}(t),t).
```
Let's verify that ``N = Z B``,
```@example s1
N = netstoichmat(rn)
N == Z*B
```

Reaction complexes give an alternative way to visualize a reaction network
graph. Catalyst's [`complexgraph`](@ref) command will calculate the complexes of
a network and then show how they are related. For example, for the repressilator
we find
```julia
complexgraph(repressilator)
```
Here ∅ represents the empty complex, black arrows show reactions converting
substrate complexes into product complexes where the rate is just a number or
parameter, and red arrows indicate conversion of substrate complexes into
product complexes where the rate is an expression involving chemical species.

![Repressilator complex](../assets/repressilator_complexgraph.svg)

## Aspects of Reaction Network Structure
Lets assume another example to highlight some more aspects of the reaction networks, using some more API functions
```@example s1
rn = @reaction_network begin
     (k1,k2), A + B <--> C
     k3, C --> D+E
     (k4,k5), D+E <--> F
     (k6,k7), 2A <--> B+G
     k8, B+G --> H
     k9, H --> 2A
end k1 k2 k3 k4 k5 k6 k7 k8 k9;

# compute incidence matrix as previously mentioned
B = reactioncomplexes(rn)[2];
```
```julia
# and visualise the graph using Graphviz based API function
complexgraph(rn)
```
![network_1](../assets/complex_rn.svg)

Lets compute a Graph from incidence matrix
```@example s1
incidencegraph = incidencematgraph(B);
```
#### Linkage classes and subnetworks from the reaction network
Figure above shows that the `ReactionSystem` `rn` is composed of two distinct “pieces,” one containing the complexes `A+B`, `C`,` D+E`, and
`F`, the other containing the complexes `2A`, `B + G`, and` H`. These two pieces can be viewed as non-linked set of complexes. These sets {`A+B`, `C`,` D+E`, `F`} and {`2A`, `B + G`,`H`} are the "linkage classes"of the reaction network.THe function `linkageclasses` returns vector of the indices of reaction complexes participating in each set of linkage-classes(Note: indices of reaction complexes found from function `reactioncomplexes` as previously explained)
```@example s1
lcs = linkageclasses(incidencegraph)
```
This however, does not tell us what are the individual `ReactionSystem`'s that form these linkage classes. This can be deduced from function `subnetworks` as follows. `subnetwork` returns a vector of `ReactionSystems`forming linkage classes.
```@example s1
subnets = subnetworks(rn, lcs)
# check the reactions in each subnetworks
reactions.(subnets)
```
```julia
 # or visualise them as
 complexgraph(subnets[1])
```
![subnetwork_1](../assets/complex_subnets1.svg)
and,
```julia
 complexgraph(subnets[2])
```
![subnetwork_2](../assets/complex_subnets2.svg)

#### Deficiency of the network
The rank of reaction network is defined as the subspace spanned by the net-stoichiometry of reaction-network. In other words, the number of uniquely represented "reactions vectors"(or the columns of net-stoichiometric matrix) is the rank of the reaction network, refer Feinberg [(1)](https://link.springer.com/book/10.1007/978-3-030-03858-8?noAccess=true).
This can be calculated as follows
```@example s1
using LinearAlgebra
s = rank(netstoichmat(rn))
```
Feinberg [(1)](https://link.springer.com/book/10.1007/978-3-030-03858-8?noAccess=true) shows that number of these uniquely represented "reaction vectors" cannot exceed the `no. of complexes - no. of linkage classes`. This puts an upper bound on the rank of reaction network, and allows us to define the deficiency of the reaction network
`δ = no. of complexes - no. of linkage classes - s`
This gives us a measure of how independent the reaction vectors are, provided the network’s linkage classes.
```@example s1
# crude way...note B i.e. incidence matrix has size of no. of complexes x no. of reactions
δ = size(B,1) - length(lcs) - rank(netstoichmat(rn))
```
or using Catalyst API function `deficiency`.
```@example s1
δ = deficiency(netstoichmat(rn), incidencegraph, lcs)
```
We may also define deficiencies for individual subnetworks in the linkage classes as follows,
```@example s1
linkage_δ = linkagedeficiencies(subnets, lcs)
```
It follows linear algebra that ,`∑ (linkage_δ) <= δ`

Quoting Feinberg [(1)](https://link.springer.com/book/10.1007/978-3-030-03858-8?noAccess=true),

> Deficiency zero networks are ones for which the reaction vectors are as independent as the partition of complexes into linkage classes will allow. And, any subnetwork of a deficiency zero network is also a deficiency zero network.

#### Reversibility of the network
Simply put, reaction network is defined reversible if, the "arrows" of the reactions are symmetric such that every reaction is accompanied by its backward reaction. Catalyst API provides function `isreversible` to determine whether the reaction network is reversible or not, based on the Graph constructed from incidence matrix for interaction of reaction complexes.
Example
```@example s1
rn = @reaction_network begin
  (k1,k2),A <--> B
  (k3,k4),A + C <--> D
  (k5,k6),D <--> B+E
  (k7,k8),B+E <--> A+C
end k1 k2 k3 k4 k5 k6 k7 k8;
# find the graph from incidence matrix for reaction complexes
incidencegraph = incidencematgraph(reactioncomplexes(rn)[2]);
isreversible(incidencegraph)
```
Consider another example,
```@example s1
rn = @reaction_network begin
  (k1,k2),A <--> B
  k3, A + C --> D
  k4, D --> B+E
  k5 ,B+E --> A+C
end k1 k2 k3 k4 k5;
incidencegraph = incidencematgraph(reactioncomplexes(rn)[2]);
isreversible(incidencegraph)
```
```julia
complexgraph(rn)
```
![reversibility](../assets/complex_reversibility.svg)

It is evident from the Figure above that the network is not "reversible", but has some sense of reversibility, such that each reaction-complex "ultimately reacts"(or indirectly reacts) to every other complex. This is known as "weakly reversible" system. One can test the network for "weak" reversibility by using function `isweaklyreversible` as follows
```@example s1
# need subnetworks from the reaction network first
B=reactioncomplexes(rn)[2]
lcs = linkageclasses(incidencematgraph(B));
subnets = subnetworks(rn, lcs)
isweaklyreversible(subnets)
```
Needless to say, every "reversible" network is also "weakly reversible", but the vice versa may or may not be true.

## Caching of Network Properties in `ReactionSystems`

## Sources
1) [Feinberg, M. *Foundations of Chemical Reaction Network Theory*, Applied Mathematical Sciences 202, Springer (2019).](https://link.springer.com/book/10.1007/978-3-030-03858-8?noAccess=true)
