using DiffEqBiological, Test, Latexify

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

@test latexify(r; noise=true, cdot=false) ==
raw"\begin{align}
\mathrm{dX1}\left( t \right) =& \left( \frac{v1 K1^{n1}}{K1^{n1} + X2^{n1}} \frac{v1 X4^{n1}}{K1^{n1} + X4^{n1}} + k1 X2 - \frac{k2}{2} X1 X4^{2} - \frac{k5}{6} X5^{3} X1 + k6 X2 - d1 X1 \right) dt + \sqrt{\left\|\frac{v1 K1^{n1}}{K1^{n1} + X2^{n1}} \frac{v1 X4^{n1}}{K1^{n1} + X4^{n1}}\right\|} dW_{1(t)} + \sqrt{\left\|k1 X2\right\|} dW_{6(t)} + \left(  - \sqrt{\left\|\frac{k2}{2} X1 X4^{2}\right\|} \right) dW_{7(t)} + \left(  - \sqrt{\left\|\frac{k5}{6} X5^{3} X1\right\|} \right) dW_{10(t)} + \sqrt{\left\|k6 X2\right\|} dW_{11(t)} + \left(  - \sqrt{\left\|d1 X1\right\|} \right) dW_{12(t)} \\
\mathrm{dX2}\left( t \right) =& \left( \frac{v2 X5^{n2}}{K2^{n2} + X5^{n2}} - k1 X2 + \frac{k2}{2} X1 X4^{2} + \frac{k5}{6} X5^{3} X1 - k6 X2 - d2 X2 \right) dt + \sqrt{\left\|\frac{v2 X5^{n2}}{K2^{n2} + X5^{n2}}\right\|} dW_{2(t)} + \left(  - \sqrt{\left\|k1 X2\right\|} \right) dW_{6(t)} + \sqrt{\left\|\frac{k2}{2} X1 X4^{2}\right\|} dW_{7(t)} + \sqrt{\left\|\frac{k5}{6} X5^{3} X1\right\|} dW_{10(t)} + \left(  - \sqrt{\left\|k6 X2\right\|} \right) dW_{11(t)} + \left(  - \sqrt{\left\|d2 X2\right\|} \right) dW_{13(t)} \\
\mathrm{dX3}\left( t \right) =& \left( \frac{v3 X3^{n3}}{K3^{n3} + X3^{n3}} + k3 X4 - k4 X3 - d3 X3 \right) dt + \sqrt{\left\|\frac{v3 X3^{n3}}{K3^{n3} + X3^{n3}}\right\|} dW_{3(t)} + \sqrt{\left\|k3 X4\right\|} dW_{8(t)} + \left(  - \sqrt{\left\|k4 X3\right\|} \right) dW_{9(t)} + \left(  - \sqrt{\left\|d3 X3\right\|} \right) dW_{14(t)} \\
\mathrm{dX4}\left( t \right) =& \left( \frac{v4 K4^{n4}}{K4^{n4} + X1^{n4}} + 2 k1 X2 -2 \frac{k2}{2} X1 X4^{2} - k3 X4 + k4 X3 - d4 X4 \right) dt + \sqrt{\left\|\frac{v4 K4^{n4}}{K4^{n4} + X1^{n4}}\right\|} dW_{4(t)} + 2 \sqrt{\left\|k1 X2\right\|} dW_{6(t)} -2 \sqrt{\left\|\frac{k2}{2} X1 X4^{2}\right\|} dW_{7(t)} + \left(  - \sqrt{\left\|k3 X4\right\|} \right) dW_{8(t)} + \sqrt{\left\|k4 X3\right\|} dW_{9(t)} + \left(  - \sqrt{\left\|d4 X4\right\|} \right) dW_{15(t)} \\
\mathrm{dX5}\left( t \right) =& \left( \frac{v5 X2^{n5}}{K5^{n5} + X2^{n5}} -3 \frac{k5}{6} X5^{3} X1 + 3 k6 X2 - d5 X5 \right) dt + \sqrt{\left\|\frac{v5 X2^{n5}}{K5^{n5} + X2^{n5}}\right\|} dW_{5(t)} -3 \sqrt{\left\|\frac{k5}{6} X5^{3} X1\right\|} dW_{10(t)} + 3 \sqrt{\left\|k6 X2\right\|} dW_{11(t)} + \left(  - \sqrt{\left\|d5 X5\right\|} \right) dW_{16(t)}
\end{align}
"


r = @reaction_network begin
    (hill(B, p_a, k, n), d_a), 0 ↔ A
    (p_b, d_b), 0 ↔ B
    (r_a, r_b), 3B ↔ A
end p_a k n d_a p_b d_b r_a r_b


@test latexify(r; noise=true, bracket=true) ==
raw"\begin{align}
\mathrm{\mathrm{d}\left[A\right]}\left( t \right) =& \left( \frac{p_{a} \cdot \left[ B \right]^{n}}{k^{n} + \left[ B \right]^{n}} - d_{a} \cdot \left[ A \right] + \frac{r_{a}}{6} \cdot \left[ B \right]^{3} - r_{b} \cdot \left[ A \right] \right) \cdot dt + \sqrt{\left\|\frac{p_{a} \cdot \left[ B \right]^{n}}{k^{n} + \left[ B \right]^{n}}\right\|} \cdot dW_{1(t)} + \left(  - \sqrt{\left\|d_{a} \cdot \left[ A \right]\right\|} \right) \cdot dW_{2(t)} + \sqrt{\left\|\frac{r_{a}}{6} \cdot \left[ B \right]^{3}\right\|} \cdot dW_{5(t)} + \left(  - \sqrt{\left\|r_{b} \cdot \left[ A \right]\right\|} \right) \cdot dW_{6(t)} \\
\mathrm{\mathrm{d}\left[B\right]}\left( t \right) =& \left( p_{b} - d_{b} \cdot \left[ B \right] -3 \cdot \frac{r_{a}}{6} \cdot \left[ B \right]^{3} + 3 \cdot r_{b} \cdot \left[ A \right] \right) \cdot dt + \sqrt{\left\|p_{b}\right\|} \cdot dW_{3(t)} + \left(  - \sqrt{\left\|d_{b} \cdot \left[ B \right]\right\|} \right) \cdot dW_{4(t)} -3 \cdot \sqrt{\left\|\frac{r_{a}}{6} \cdot \left[ B \right]^{3}\right\|} \cdot dW_{5(t)} + 3 \cdot \sqrt{\left\|r_{b} \cdot \left[ A \right]\right\|} \cdot dW_{6(t)}
\end{align}
"


@test latexify(r; noise_only=true, bracket=true) ==
raw"\begin{align}
\mathrm{\mathrm{d}\left[A\right]}\left( t \right) ∝& \sqrt{\left\|\frac{p_{a} \cdot \left[ B \right]^{n}}{k^{n} + \left[ B \right]^{n}}\right\|} \cdot dW_{1(t)} + \left(  - \sqrt{\left\|d_{a} \cdot \left[ A \right]\right\|} \right) \cdot dW_{2(t)} + \sqrt{\left\|\frac{r_{a}}{6} \cdot \left[ B \right]^{3}\right\|} \cdot dW_{5(t)} + \left(  - \sqrt{\left\|r_{b} \cdot \left[ A \right]\right\|} \right) \cdot dW_{6(t)} \\
\mathrm{\mathrm{d}\left[B\right]}\left( t \right) ∝& \sqrt{\left\|p_{b}\right\|} \cdot dW_{3(t)} + \left(  - \sqrt{\left\|d_{b} \cdot \left[ B \right]\right\|} \right) \cdot dW_{4(t)} -3 \cdot \sqrt{\left\|\frac{r_{a}}{6} \cdot \left[ B \right]^{3}\right\|} \cdot dW_{5(t)} + 3 \cdot \sqrt{\left\|r_{b} \cdot \left[ A \right]\right\|} \cdot dW_{6(t)}
\end{align}
"


@test latexify(r; noise_only=true, bracket=true, cdot=false) ==
raw"\begin{align}
\mathrm{\mathrm{d}\left[A\right]}\left( t \right) ∝& \sqrt{\left\|\frac{p_{a} \left[ B \right]^{n}}{k^{n} + \left[ B \right]^{n}}\right\|} dW_{1(t)} + \left(  - \sqrt{\left\|d_{a} \left[ A \right]\right\|} \right) dW_{2(t)} + \sqrt{\left\|\frac{r_{a}}{6} \left[ B \right]^{3}\right\|} dW_{5(t)} + \left(  - \sqrt{\left\|r_{b} \left[ A \right]\right\|} \right) dW_{6(t)} \\
\mathrm{\mathrm{d}\left[B\right]}\left( t \right) ∝& \sqrt{\left\|p_{b}\right\|} dW_{3(t)} + \left(  - \sqrt{\left\|d_{b} \left[ B \right]\right\|} \right) dW_{4(t)} -3 \sqrt{\left\|\frac{r_{a}}{6} \left[ B \right]^{3}\right\|} dW_{5(t)} + 3 \sqrt{\left\|r_{b} \left[ A \right]\right\|} dW_{6(t)}
\end{align}
"


@test latexify(r; noise=true, noise_var=:Noise, bracket=true) ==
raw"\begin{align}
\mathrm{\mathrm{d}\left[A\right]}\left( t \right) =& \left( \frac{p_{a} \cdot \left[ B \right]^{n}}{k^{n} + \left[ B \right]^{n}} - d_{a} \cdot \left[ A \right] + \frac{r_{a}}{6} \cdot \left[ B \right]^{3} - r_{b} \cdot \left[ A \right] \right) \cdot dt + \sqrt{\left\|\frac{p_{a} \cdot \left[ B \right]^{n}}{k^{n} + \left[ B \right]^{n}}\right\|} \cdot dNoise_{1(t)} + \left(  - \sqrt{\left\|d_{a} \cdot \left[ A \right]\right\|} \right) \cdot dNoise_{2(t)} + \sqrt{\left\|\frac{r_{a}}{6} \cdot \left[ B \right]^{3}\right\|} \cdot dNoise_{5(t)} + \left(  - \sqrt{\left\|r_{b} \cdot \left[ A \right]\right\|} \right) \cdot dNoise_{6(t)} \\
\mathrm{\mathrm{d}\left[B\right]}\left( t \right) =& \left( p_{b} - d_{b} \cdot \left[ B \right] -3 \cdot \frac{r_{a}}{6} \cdot \left[ B \right]^{3} + 3 \cdot r_{b} \cdot \left[ A \right] \right) \cdot dt + \sqrt{\left\|p_{b}\right\|} \cdot dNoise_{3(t)} + \left(  - \sqrt{\left\|d_{b} \cdot \left[ B \right]\right\|} \right) \cdot dNoise_{4(t)} -3 \cdot \sqrt{\left\|\frac{r_{a}}{6} \cdot \left[ B \right]^{3}\right\|} \cdot dNoise_{5(t)} + 3 \cdot \sqrt{\left\|r_{b} \cdot \left[ A \right]\right\|} \cdot dNoise_{6(t)}
\end{align}
"


@test latexify(r; bracket=true) ==
raw"\begin{align}
\frac{d\left[ A \right](t)}{dt} =& \frac{p_{a} \cdot \left[ B \right]^{n}}{k^{n} + \left[ B \right]^{n}} - d_{a} \cdot \left[ A \right] + \frac{r_{a}}{6} \cdot \left[ B \right]^{3} - r_{b} \cdot \left[ A \right] \\
\frac{d\left[ B \right](t)}{dt} =& p_{b} - d_{b} \cdot \left[ B \right] -3 \cdot \frac{r_{a}}{6} \cdot \left[ B \right]^{3} + 3 \cdot r_{b} \cdot \left[ A \right]
\end{align}
"


@test latexify(r; bracket=false) ==
raw"\begin{align}
\frac{dA(t)}{dt} =& \frac{p_{a} \cdot B^{n}}{k^{n} + B^{n}} - d_{a} \cdot A + \frac{r_{a}}{6} \cdot B^{3} - r_{b} \cdot A \\
\frac{dB(t)}{dt} =& p_{b} - d_{b} \cdot B -3 \cdot \frac{r_{a}}{6} \cdot B^{3} + 3 \cdot r_{b} \cdot A
\end{align}
"
