#! format: off

### Fetch required packages and reaction networks ###
using Catalyst, Latexify

############################
### CURRENTLY NOT ACITVE ###
### REQUIRES REWRITING   ###
############################

### Tips for generating latex tests:
### Latexify has an unexported macro:
###
### Latexify.@generate_test
###
### which generates a test using a given latexify function.
### For example:
###
### Latexify.@generate_test latexify([1, 2, 3], [4, 5, 6]; env=:mdtable)
###
### This puts a ready-made test in your clipboard which you can paste into the
### test file.
###
### Just be sure to remove all such macros before you commit a change since it
### will cause issues with Travis.

r = @reaction_network begin
    @parameters v1 K1 n1 v2 K2 n2 v3 K3 n3 v4 K4 n4 v5 K5 n5 k1 k2 k3 k4 k5 k6 d1 d2 d3 d4 d5
    hillr(X2,v1,K1,n1)*hill(X4,v1,K1,n1), ∅ → X1
    hill(X5,v2,K2,n2), ∅ → X2
    hill(X3,v3,K3,n3), ∅ → X3
    hillr(X1,v4,K4,n4), ∅ → X4
    hill(X2,v5,K5,n5), ∅ → X5
    (k1,k2), X2 ⟷ X1 + 2X4
    (k3,k4), X4 ⟷ X3
    (k5,k6), 3X5 + X1 ⟷ X2
    (d1,d2,d3,d4,d5), (X1,X2,X3,X4,X5)  ⟶ ∅
end

# Latexify.@generate_test latexify(r)
@test latexify(r) == replace(
raw"\begin{align*}
\require{mhchem}
\ce{ \varnothing &->[$\frac{v1 X4^{n1}}{K1^{n1} + X4^{n1}} \frac{v1 K1^{n1}}{K1^{n1} + X2^{n1}}$] X1}\\
\ce{ \varnothing &->[$\frac{v2 X5^{n2}}{K2^{n2} + X5^{n2}}$] X2}\\
\ce{ \varnothing &->[$\frac{v3 X3^{n3}}{K3^{n3} + X3^{n3}}$] X3}\\
\ce{ \varnothing &->[$\frac{v4 K4^{n4}}{K4^{n4} + X1^{n4}}$] X4}\\
\ce{ \varnothing &->[$\frac{v5 X2^{n5}}{K5^{n5} + X2^{n5}}$] X5}\\
\ce{ X2 &<=>[$k1$][$k2$] X1 + 2 X4}\\
\ce{ X4 &<=>[$k3$][$k4$] X3}\\
\ce{ 3 X5 + X1 &<=>[$k5$][$k6$] X2}\\
\ce{ X1 &->[$d1$] \varnothing}\\
\ce{ X2 &->[$d2$] \varnothing}\\
\ce{ X3 &->[$d3$] \varnothing}\\
\ce{ X4 &->[$d4$] \varnothing}\\
\ce{ X5 &->[$d5$] \varnothing}
\end{align*}
", "\r\n"=>"\n")

# Latexify.@generate_test latexify(r, mathjax=false)
@test latexify(r, mathjax = false) == replace(
raw"\begin{align*}
\ce{ \varnothing &->[$\frac{v1 X4^{n1}}{K1^{n1} + X4^{n1}} \frac{v1 K1^{n1}}{K1^{n1} + X2^{n1}}$] X1}\\
\ce{ \varnothing &->[$\frac{v2 X5^{n2}}{K2^{n2} + X5^{n2}}$] X2}\\
\ce{ \varnothing &->[$\frac{v3 X3^{n3}}{K3^{n3} + X3^{n3}}$] X3}\\
\ce{ \varnothing &->[$\frac{v4 K4^{n4}}{K4^{n4} + X1^{n4}}$] X4}\\
\ce{ \varnothing &->[$\frac{v5 X2^{n5}}{K5^{n5} + X2^{n5}}$] X5}\\
\ce{ X2 &<=>[$k1$][$k2$] X1 + 2 X4}\\
\ce{ X4 &<=>[$k3$][$k4$] X3}\\
\ce{ 3 X5 + X1 &<=>[$k5$][$k6$] X2}\\
\ce{ X1 &->[$d1$] \varnothing}\\
\ce{ X2 &->[$d2$] \varnothing}\\
\ce{ X3 &->[$d3$] \varnothing}\\
\ce{ X4 &->[$d4$] \varnothing}\\
\ce{ X5 &->[$d5$] \varnothing}
\end{align*}
", "\r\n"=>"\n")


r = @reaction_network begin
    @parameters p_a k n d_a p_b d_b r_a r_b
    (hill(B, p_a, k, n), d_a), 0 ↔ A
    (p_b, d_b), 0 ↔ B
    (r_a, r_b), 3B ↔ A
end

# Latexify.@generate_test latexify(r)
@test latexify(r) == replace(
raw"\begin{align*}
\require{mhchem}
\ce{ \varnothing &<=>[$\frac{p_{a} B^{n}}{k^{n} + B^{n}}$][$d_{a}$] A}\\
\ce{ \varnothing &<=>[$p_{b}$][$d_{b}$] B}\\
\ce{ 3 B &<=>[$r_{a}$][$r_{b}$] A}
\end{align*}
", "\r\n"=>"\n")

# Latexify.@generate_test latexify(r, mathjax=false)
@test latexify(r, mathjax = false) == replace(
raw"\begin{align*}
\ce{ \varnothing &<=>[$\frac{p_{a} B^{n}}{k^{n} + B^{n}}$][$d_{a}$] A}\\
\ce{ \varnothing &<=>[$p_{b}$][$d_{b}$] B}\\
\ce{ 3 B &<=>[$r_{a}$][$r_{b}$] A}
\end{align*}
", "\r\n"=>"\n")

# test empty system
empty_rn = ReactionSystem(Reaction[]; name=:EmptySys)
# Latexify.@generate_test latexify(empty_rn)
@test latexify(empty_rn) == replace(
raw"ReactionSystem EmptySys has no reactions.", "\r\n"=>"\n")

# test for https://github.com/SciML/Catalyst.jl/issues/473
rn = @reaction_network begin
    @parameters k
    k*Y, Y --> ∅
end

# Latexify.@generate_test latexify(rn)
@test latexify(rn) == replace(
raw"\begin{align*}
\require{mhchem}
\ce{ Y &->[$Y k$] \varnothing}
\end{align*}
", "\r\n"=>"\n")

# Tests the `form` option
for rn in reaction_networks_standard
    @test latexify(rn)==latexify(rn; form=:reactions)
end

rn = @reaction_network begin 
    (p,d), 0 <--> X
    (kB,kD), 2X <--> X2
end
@test latexify(rn; form=:ode) == raw"$\begin{align}
\frac{\mathrm{d} X\left( t \right)}{\mathrm{d}t} =& p - \left( X\left( t \right) \right)^{2} kB - d X\left( t \right) + 2 kD \mathrm{X2}\left( t \right) \\
\frac{\mathrm{d} \mathrm{X2}\left( t \right)}{\mathrm{d}t} =& \frac{1}{2} \left( X\left( t \right) \right)^{2} kB - kD \mathrm{X2}\left( t \right)
\end{align}
$"