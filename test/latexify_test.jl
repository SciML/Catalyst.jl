using Test
using Latexify
using DiffEqBiological

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


rn = @reaction_network begin
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


@test latexify(rn; noise=true) == 
raw"\begin{align}
\frac{dX1}{dt} =& \frac{v1 \cdot K1^{n1}}{K1^{n1} + X2^{n1}} \cdot \frac{v1 \cdot X4^{n1}}{K1^{n1} + X4^{n1}} + k1 \cdot X2 - \frac{k2}{2} \cdot X1 \cdot X4^{2} - \frac{k5}{6} \cdot X5^{3} \cdot X1 + k6 \cdot X2 - d1 \cdot X1 + \sqrt{\left\|\frac{v1 \cdot K1^{n1}}{K1^{n1} + X2^{n1}} \cdot \frac{v1 \cdot X4^{n1}}{K1^{n1} + X4^{n1}}\right\|} \cdot W_{1} + \sqrt{\left\|k1 \cdot X2\right\|} \cdot W_{6} - \sqrt{\left\|\frac{k2}{2} \cdot X1 \cdot X4^{2}\right\|} \cdot W_{7} - \sqrt{\left\|\frac{k5}{6} \cdot X5^{3} \cdot X1\right\|} \cdot W_{10} + \sqrt{\left\|k6 \cdot X2\right\|} \cdot W_{11} - \sqrt{\left\|d1 \cdot X1\right\|} \cdot W_{12} \\
\frac{dX2}{dt} =& \frac{v2 \cdot X5^{n2}}{K2^{n2} + X5^{n2}} - k1 \cdot X2 + \frac{k2}{2} \cdot X1 \cdot X4^{2} + \frac{k5}{6} \cdot X5^{3} \cdot X1 - k6 \cdot X2 - d2 \cdot X2 + \sqrt{\left\|\frac{v2 \cdot X5^{n2}}{K2^{n2} + X5^{n2}}\right\|} \cdot W_{2} - \sqrt{\left\|k1 \cdot X2\right\|} \cdot W_{6} + \sqrt{\left\|\frac{k2}{2} \cdot X1 \cdot X4^{2}\right\|} \cdot W_{7} + \sqrt{\left\|\frac{k5}{6} \cdot X5^{3} \cdot X1\right\|} \cdot W_{10} - \sqrt{\left\|k6 \cdot X2\right\|} \cdot W_{11} - \sqrt{\left\|d2 \cdot X2\right\|} \cdot W_{13} \\
\frac{dX3}{dt} =& \frac{v3 \cdot X3^{n3}}{K3^{n3} + X3^{n3}} + k3 \cdot X4 - k4 \cdot X3 - d3 \cdot X3 + \sqrt{\left\|\frac{v3 \cdot X3^{n3}}{K3^{n3} + X3^{n3}}\right\|} \cdot W_{3} + \sqrt{\left\|k3 \cdot X4\right\|} \cdot W_{8} - \sqrt{\left\|k4 \cdot X3\right\|} \cdot W_{9} - \sqrt{\left\|d3 \cdot X3\right\|} \cdot W_{14} \\
\frac{dX4}{dt} =& \frac{v4 \cdot K4^{n4}}{K4^{n4} + X1^{n4}} + 2 \cdot k1 \cdot X2 -2 \cdot \frac{k2}{2} \cdot X1 \cdot X4^{2} - k3 \cdot X4 + k4 \cdot X3 - d4 \cdot X4 + \sqrt{\left\|\frac{v4 \cdot K4^{n4}}{K4^{n4} + X1^{n4}}\right\|} \cdot W_{4} + 2 \cdot \sqrt{\left\|k1 \cdot X2\right\|} \cdot W_{6} -2 \cdot \sqrt{\left\|\frac{k2}{2} \cdot X1 \cdot X4^{2}\right\|} \cdot W_{7} - \sqrt{\left\|k3 \cdot X4\right\|} \cdot W_{8} + \sqrt{\left\|k4 \cdot X3\right\|} \cdot W_{9} - \sqrt{\left\|d4 \cdot X4\right\|} \cdot W_{15} \\
\frac{dX5}{dt} =& \frac{v5 \cdot X2^{n5}}{K5^{n5} + X2^{n5}} -3 \cdot \frac{k5}{6} \cdot X5^{3} \cdot X1 + 3 \cdot k6 \cdot X2 - d5 \cdot X5 + \sqrt{\left\|\frac{v5 \cdot X2^{n5}}{K5^{n5} + X2^{n5}}\right\|} \cdot W_{5} -3 \cdot \sqrt{\left\|\frac{k5}{6} \cdot X5^{3} \cdot X1\right\|} \cdot W_{10} + 3 \cdot \sqrt{\left\|k6 \cdot X2\right\|} \cdot W_{11} - \sqrt{\left\|d5 \cdot X5\right\|} \cdot W_{16}
\end{align}
"

r = @reaction_network begin
    (hill(B, p_a, k, n), d_a), 0 ↔ A
    (p_b, d_b), 0 ↔ B
    (r_a, r_b), 3B ↔ A
end p_a k n d_a p_b d_b r_a r_b


@test latexify(r; noise=true) == 
raw"\begin{align}
\frac{dA}{dt} =& \frac{p_{a} \cdot B^{n}}{k^{n} + B^{n}} - d_{a} \cdot A + \frac{r_{a}}{6} \cdot B^{3} - r_{b} \cdot A + \sqrt{\left\|\frac{p_{a} \cdot B^{n}}{k^{n} + B^{n}}\right\|} \cdot W_{1} - \sqrt{\left\|d_{a} \cdot A\right\|} \cdot W_{2} + \sqrt{\left\|\frac{r_{a}}{6} \cdot B^{3}\right\|} \cdot W_{5} - \sqrt{\left\|r_{b} \cdot A\right\|} \cdot W_{6} \\
\frac{dB}{dt} =& p_{b} - d_{b} \cdot B -3 \cdot \frac{r_{a}}{6} \cdot B^{3} + 3 \cdot r_{b} \cdot A + \sqrt{\left\|p_{b}\right\|} \cdot W_{3} - \sqrt{\left\|d_{b} \cdot B\right\|} \cdot W_{4} -3 \cdot \sqrt{\left\|\frac{r_{a}}{6} \cdot B^{3}\right\|} \cdot W_{5} + 3 \cdot \sqrt{\left\|r_{b} \cdot A\right\|} \cdot W_{6}
\end{align}
"

@test latexify(r; noise_only=true) == 
raw"\begin{align}
\frac{dA}{dt} =& \sqrt{\left\|\frac{p_{a} \cdot B^{n}}{k^{n} + B^{n}}\right\|} \cdot W_{1} - \sqrt{\left\|d_{a} \cdot A\right\|} \cdot W_{2} + \sqrt{\left\|\frac{r_{a}}{6} \cdot B^{3}\right\|} \cdot W_{5} - \sqrt{\left\|r_{b} \cdot A\right\|} \cdot W_{6} \\
\frac{dB}{dt} =& \sqrt{\left\|p_{b}\right\|} \cdot W_{3} - \sqrt{\left\|d_{b} \cdot B\right\|} \cdot W_{4} -3 \cdot \sqrt{\left\|\frac{r_{a}}{6} \cdot B^{3}\right\|} \cdot W_{5} + 3 \cdot \sqrt{\left\|r_{b} \cdot A\right\|} \cdot W_{6}
\end{align}
"

@test latexify(r; noise=true, noise_var=:Noise, bracket=true) == 
raw"\begin{align}
\frac{\mathrm{d}\left[A\right]}{dt} =& \frac{p_{a} \cdot \left[ B \right]^{n}}{k^{n} + \left[ B \right]^{n}} - d_{a} \cdot \left[ A \right] + \frac{r_{a}}{6} \cdot \left[ B \right]^{3} - r_{b} \cdot \left[ A \right] + \sqrt{\left\|\frac{p_{a} \cdot \left[ B \right]^{n}}{k^{n} + \left[ B \right]^{n}}\right\|} \cdot Noise_{1} - \sqrt{\left\|d_{a} \cdot \left[ A \right]\right\|} \cdot Noise_{2} + \sqrt{\left\|\frac{r_{a}}{6} \cdot \left[ B \right]^{3}\right\|} \cdot Noise_{5} - \sqrt{\left\|r_{b} \cdot \left[ A \right]\right\|} \cdot Noise_{6} \\
\frac{\mathrm{d}\left[B\right]}{dt} =& p_{b} - d_{b} \cdot \left[ B \right] -3 \cdot \frac{r_{a}}{6} \cdot \left[ B \right]^{3} + 3 \cdot r_{b} \cdot \left[ A \right] + \sqrt{\left\|p_{b}\right\|} \cdot Noise_{3} - \sqrt{\left\|d_{b} \cdot \left[ B \right]\right\|} \cdot Noise_{4} -3 \cdot \sqrt{\left\|\frac{r_{a}}{6} \cdot \left[ B \right]^{3}\right\|} \cdot Noise_{5} + 3 \cdot \sqrt{\left\|r_{b} \cdot \left[ A \right]\right\|} \cdot Noise_{6}
\end{align}
"

