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
    hillR(X2,v1,K1,n1)*hill(X4,v1,K1,n1), ∅ → X1
    hill(X5,v2,K2,n2), ∅ → X2
    hill(X3,v3,K3,n3), ∅ → X3
    hillR(X1,v4,K4,n4), ∅ → X4
    hill(X2,v5,K5,n5), ∅ → X5
    (k1,k2), X2 ⟷ X1 + 2X4
    (k3,k4), X4 ⟷ X3
    (k5,k6), 3X5 + X1 ⟷ X2
    (d1,d2,d3,d4,d5), (X1,X2,X3,X4,X5)  ⟶ ∅
end v1 K1 n1 v2 K2 n2 v3 K3 n3 v4 K4 n4 v5 K5 n5 k1 k2 k3 k4 k5 k6 d1 d2 d3 d4 d5

# Latexify.@generate_test latexify(r)
@test latexify(r) == replace(
raw"\begin{align}
\require{mhchem}
\ce{ \varnothing &->[K1^{n1} v1^{2} \left( \mathrm{X4}\left( t \right) \right)^{n1} \mathrm{inv}\left( K1^{n1} + \left( \mathrm{X2}\left( t \right) \right)^{n1} \right) \mathrm{inv}\left( K1^{n1} + \left( \mathrm{X4}\left( t \right) \right)^{n1} \right)] X1}\\
\ce{ \varnothing &->[v2 \left( \mathrm{X5}\left( t \right) \right)^{n2} \mathrm{inv}\left( K2^{n2} + \left( \mathrm{X5}\left( t \right) \right)^{n2} \right)] X2}\\
\ce{ \varnothing &->[v3 \left( \mathrm{X3}\left( t \right) \right)^{n3} \mathrm{inv}\left( K3^{n3} + \left( \mathrm{X3}\left( t \right) \right)^{n3} \right)] X3}\\
\ce{ \varnothing &->[v4 K4^{n4} \mathrm{inv}\left( K4^{n4} + \left( \mathrm{X1}\left( t \right) \right)^{n4} \right)] X4}\\
\ce{ \varnothing &->[v5 \left( \mathrm{X2}\left( t \right) \right)^{n5} \mathrm{inv}\left( K5^{n5} + \left( \mathrm{X2}\left( t \right) \right)^{n5} \right)] X5}\\
\ce{ X2 &<=>[k1][k2] X1 + 2 X4}\\
\ce{ X4 &<=>[k3][k4] X3}\\
\ce{ 3 X5 + X1 &<=>[k5][k6] X2}\\
\ce{ X1 &->[d1] \varnothing}\\
\ce{ X2 &->[d2] \varnothing}\\
\ce{ X3 &->[d3] \varnothing}\\
\ce{ X4 &->[d4] \varnothing}\\
\ce{ X5 &->[d5] \varnothing}
\end{align}
", "\r\n"=>"\n")

# Latexify.@generate_test latexify(r, mathjax=false)
@test latexify(r, mathjax = false) == replace(
raw"\begin{align}
\ce{ \varnothing &->[$K1^{n1} v1^{2} \left( \mathrm{X4}\left( t \right) \right)^{n1} \mathrm{inv}\left( K1^{n1} + \left( \mathrm{X2}\left( t \right) \right)^{n1} \right) \mathrm{inv}\left( K1^{n1} + \left( \mathrm{X4}\left( t \right) \right)^{n1} \right)$] X1}\\
\ce{ \varnothing &->[$v2 \left( \mathrm{X5}\left( t \right) \right)^{n2} \mathrm{inv}\left( K2^{n2} + \left( \mathrm{X5}\left( t \right) \right)^{n2} \right)$] X2}\\
\ce{ \varnothing &->[$v3 \left( \mathrm{X3}\left( t \right) \right)^{n3} \mathrm{inv}\left( K3^{n3} + \left( \mathrm{X3}\left( t \right) \right)^{n3} \right)$] X3}\\
\ce{ \varnothing &->[$v4 K4^{n4} \mathrm{inv}\left( K4^{n4} + \left( \mathrm{X1}\left( t \right) \right)^{n4} \right)$] X4}\\
\ce{ \varnothing &->[$v5 \left( \mathrm{X2}\left( t \right) \right)^{n5} \mathrm{inv}\left( K5^{n5} + \left( \mathrm{X2}\left( t \right) \right)^{n5} \right)$] X5}\\
\ce{ X2 &<=>[$k1$][$k2$] X1 + 2 X4}\\
\ce{ X4 &<=>[$k3$][$k4$] X3}\\
\ce{ 3 X5 + X1 &<=>[$k5$][$k6$] X2}\\
\ce{ X1 &->[$d1$] \varnothing}\\
\ce{ X2 &->[$d2$] \varnothing}\\
\ce{ X3 &->[$d3$] \varnothing}\\
\ce{ X4 &->[$d4$] \varnothing}\\
\ce{ X5 &->[$d5$] \varnothing}
\end{align}
", "\r\n"=>"\n")


r = @reaction_network begin
    (hill(B, p_a, k, n), d_a), 0 ↔ A
    (p_b, d_b), 0 ↔ B
    (r_a, r_b), 3B ↔ A
end p_a k n d_a p_b d_b r_a r_b

# Latexify.@generate_test latexify(r)
@test latexify(r) == replace(
raw"\begin{align}
\require{mhchem}
\ce{ \varnothing &<=>[p_{a} \left( B\left( t \right) \right)^{n} \mathrm{inv}\left( k^{n} + \left( B\left( t \right) \right)^{n} \right)][d_{a}] A}\\
\ce{ \varnothing &<=>[p_{b}][d_{b}] B}\\
\ce{ 3 B &<=>[r_{a}][r_{b}] A}
\end{align}
", "\r\n"=>"\n")

# Latexify.@generate_test latexify(r, mathjax=false)
@test latexify(r, mathjax = false) == replace(
raw"\begin{align}
\ce{ \varnothing &<=>[$p_{a} \left( B\left( t \right) \right)^{n} \mathrm{inv}\left( k^{n} + \left( B\left( t \right) \right)^{n} \right)$][$d_{a}$] A}\\
\ce{ \varnothing &<=>[$p_{b}$][$d_{b}$] B}\\
\ce{ 3 B &<=>[$r_{a}$][$r_{b}$] A}
\end{align}
", "\r\n"=>"\n")

