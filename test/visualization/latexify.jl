#! format: off

### Fetch Packages and Reaction Networks ###
using Catalyst, Latexify
include("../test_networks.jl")

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

### Basic Tests ###

let
    rn = @reaction_network begin
        hillr(X2,v1,K1,n1)*hill(X4,v1,K1,n1), ∅ → X1
        hill(X5,v2,K2,n2), ∅ → X2
        hill(X3,v3,K3,n3), ∅ → X3
        hillr(X1,v4,K4,n4), ∅ → X4
        hill(X2,v5,K5,n5), ∅ → X5
        hillar(X1,X6,v6,K6,n6), ∅ → X6
        (k1,k2), X2 ⟷ X1 + 2X4
        (k3,k4), X4 ⟷ X3
        (k5,k6), 3X5 + X1 ⟷ X2
        (d1,d2,d3,d4,d5,d6), (X1,X2,X3,X4,X5,X6)  ⟶ ∅
    end

    # Latexify.@generate_test latexify(rn)
    @test_broken latexify(rn; expand_functions = false) == replace(
    raw"\begin{align*}
    \varnothing &\xrightarrow{\frac{X4^{n1} v1^{2} K1^{n1}}{\left( K1^{n1} + X4^{n1} \right) \left( K1^{n1} + X2^{n1} \right)}} \mathrm{X1} \\
    \varnothing &\xrightarrow{\mathrm{hill}\left( X5, v2, K2, n2 \right)} \mathrm{X2} \\
    \varnothing &\xrightarrow{\mathrm{hill}\left( X3, v3, K3, n3 \right)} \mathrm{X3} \\
    \varnothing &\xrightarrow{\mathrm{hillr}\left( X1, v4, K4, n4 \right)} \mathrm{X4} \\
    \varnothing &\xrightarrow{\mathrm{hill}\left( X2, v5, K5, n5 \right)} \mathrm{X5} \\
    \varnothing &\xrightarrow{\mathrm{hillar}\left( X1, X6, v6, K6, n6 \right)} \mathrm{X6} \\
    \mathrm{X2} &\xrightleftharpoons[k2]{k1} \mathrm{X1} + 2 \mathrm{X4} \\
    \mathrm{X4} &\xrightleftharpoons[k4]{k3} \mathrm{X3} \\
    3 \mathrm{X5} + \mathrm{X1} &\xrightleftharpoons[k6]{k5} \mathrm{X2} \\
    \mathrm{X1} &\xrightarrow{d1} \varnothing \\
    \mathrm{X2} &\xrightarrow{d2} \varnothing \\
    \mathrm{X3} &\xrightarrow{d3} \varnothing \\
    \mathrm{X4} &\xrightarrow{d4} \varnothing \\
    \mathrm{X5} &\xrightarrow{d5} \varnothing \\
    \mathrm{X6} &\xrightarrow{d6} \varnothing  
    \end{align*}
    ", "\r\n"=>"\n")

    #Latexify.@generate_test latexify(rn; expand_functions=false)
    @test_broken latexify(rn; expand_functions = false) == replace(
    raw"\begin{align*}
    \varnothing &\xrightarrow{\frac{X4^{n1} v1^{2} K1^{n1}}{\left( K1^{n1} + X4^{n1} \right) \left( K1^{n1} + X2^{n1} \right)}} \mathrm{X1} \\
    \varnothing &\xrightarrow{\mathrm{mm}\left( X5, v2, K2 \right)} \mathrm{X2} \\
    \varnothing &\xrightarrow{\mathrm{mmr}\left( X3, v3, K3 \right)} \mathrm{X3} \\
    \varnothing &\xrightarrow{\mathrm{hillr}\left( X1, v4, K4, n4 \right)} \mathrm{X4} \\
    \varnothing &\xrightarrow{\mathrm{hill}\left( X2, v5, K5, n5 \right)} \mathrm{X5} \\
    \varnothing &\xrightarrow{\mathrm{hillar}\left( X1, X6, v6, K6, n6 \right)} \mathrm{X6} \\
    \mathrm{X2} &\xrightleftharpoons[k2]{k1} \mathrm{X1} + 2 \mathrm{X4} \\
    \mathrm{X4} &\xrightleftharpoons[k4]{k3} \mathrm{X3} \\
    3 \mathrm{X5} + \mathrm{X1} &\xrightleftharpoons[k6]{k5} \mathrm{X2} \\
    \mathrm{X1} &\xrightarrow{d1} \varnothing \\
    \mathrm{X2} &\xrightarrow{d2} \varnothing \\
    \mathrm{X3} &\xrightarrow{d3} \varnothing \\
    \mathrm{X4} &\xrightarrow{d4} \varnothing \\
    \mathrm{X5} &\xrightarrow{d5} \varnothing \\
    \mathrm{X6} &\xrightarrow{d6} \varnothing  
    \end{align*}
    ", "\r\n"=>"\n")

    # Latexify.@generate_test latexify(rn, mathjax=false)
    @test_broken latexify(rn, mathjax = false) == replace(
    raw"\begin{align*}
    \varnothing &\xrightarrow{\frac{X4^{n1} v1^{2} K1^{n1}}{\left( K1^{n1} + X4^{n1} \right) \left( K1^{n1} + X2^{n1} \right)}} \mathrm{X1} \\
    \varnothing &\xrightarrow{\frac{X5 v2}{K2 + X5}} \mathrm{X2} \\
    \varnothing &\xrightarrow{\frac{K3 v3}{K3 + X3}} \mathrm{X3} \\
    \varnothing &\xrightarrow{\frac{v4 K4^{n4}}{K4^{n4} + X1^{n4}}} \mathrm{X4} \\
    \varnothing &\xrightarrow{\frac{v5 X2^{n5}}{X2^{n5} + K5^{n5}}} \mathrm{X5} \\
    \varnothing &\xrightarrow{\frac{v6 X1^{n6}}{X6^{n6} + K6^{n6} + X1^{n6}}} \mathrm{X6} \\
    \mathrm{X2} &\xrightleftharpoons[k2]{k1} \mathrm{X1} + 2 \mathrm{X4} \\
    \mathrm{X4} &\xrightleftharpoons[k4]{k3} \mathrm{X3} \\
    3 \mathrm{X5} + \mathrm{X1} &\xrightleftharpoons[k6]{k5} \mathrm{X2} \\
    \mathrm{X1} &\xrightarrow{d1} \varnothing \\
    \mathrm{X2} &\xrightarrow{d2} \varnothing \\
    \mathrm{X3} &\xrightarrow{d3} \varnothing \\
    \mathrm{X4} &\xrightarrow{d4} \varnothing \\
    \mathrm{X5} &\xrightarrow{d5} \varnothing \\
    \mathrm{X6} &\xrightarrow{d6} \varnothing  
    \end{align*}
    ", "\r\n"=>"\n")
end

let
    rn = @reaction_network begin
        (hill(B, p_a, k, n), d_a), 0 ↔ A
        (p_b, d_b), 0 ↔ B
        (r_a, r_b), 3B ↔ A
    end

    # Latexify.@generate_test latexify(rn)
    @test_broken latexify(rn) == replace(
    raw"\begin{align*}
    \varnothing &\xrightleftharpoons[d_{a}]{\frac{p_{a} B^{n}}{k^{n} + B^{n}}} \mathrm{A} \\
    \varnothing &\xrightleftharpoons[d_{b}]{p_{b}} \mathrm{B} \\
    3 \mathrm{B} &\xrightleftharpoons[r_{b}]{r_{a}} \mathrm{A}  
    \end{align*}
    ", "\r\n"=>"\n")

    # Latexify.@generate_test latexify(rn, mathjax=false)
    @test_broken latexify(rn, mathjax = false) == replace(
    raw"\begin{align*}
    \varnothing &\xrightleftharpoons[d_{a}]{\frac{p_{a} B^{n}}{k^{n} + B^{n}}} \mathrm{A} \\
    \varnothing &\xrightleftharpoons[d_{b}]{p_{b}} \mathrm{B} \\
    3 \mathrm{B} &\xrightleftharpoons[r_{b}]{r_{a}} \mathrm{A}  
    \end{align*}
    ", "\r\n"=>"\n")
end

# Test empty system.
let
    empty_rn = ReactionSystem(Reaction[]; name=:EmptySys)

    # Latexify.@generate_test latexify(empty_rn)
    @test latexify(empty_rn) == replace(
    raw"ReactionSystem EmptySys has no reactions or equations.", "\r\n"=>"\n")
end

# Test for https://github.com/SciML/Catalyst.jl/issues/473.
let
    rn = @reaction_network begin
        k*Y, Y --> ∅
    end

    # Latexify.@generate_test latexify(rn)
    @test_broken latexify(rn) == replace(
    raw"\begin{align*}
    \varnothing &\xrightarrow{p} (m + n)\mathrm{X}  
    \end{align*}
    ", "\r\n"=>"\n")
end


# Test using various `env` options.
let
    rn = @reaction_network begin
        (p,d), 0 <--> X
    end
    chem_latex = latexify(rn; env = :arrows)
    @test chem_latex == latexify(rn; env = :chem)
    @test chem_latex == latexify(rn; env = :chemical)
    @test chem_latex == latexify(rn; env = :arrow)
    @test_throws Exception latexify(rn; env = :wrong_env)
end

# Tests that the `mathrm` option affects the output.
let
    rn = @reaction_network begin
        (k1,k2), 2X <--> X2
    end
    @test latexify(rn; mathrm = true) != latexify(rn; mathrm = false)
end

# Tests for system with parametric stoichiometry.
let
    rn = @reaction_network begin
        p, 0 --> (m + n)*X
    end
    
    @test_broken latexify(rn) == replace(
    raw"\begin{align*}
    \varnothing &\xrightarrow{p} (m + n)\mathrm{X}  
     \end{align*}
    ", "\r\n"=>"\n")
end

### Tests  `form` Option ###

# Check for large number of networks.
let
    for rn in reaction_networks_standard
        @test latexify(rn)==latexify(rn; form=:reactions)
        #@test_broken latexify(convert(ODESystem,rn)) == latexify(rn; form=:ode) # Slight difference due to some latexify weirdly. Both displays fine though
    end
end

# Check for specific network.
let
    rn = @reaction_network begin
        (p,d), 0 <--> X
        (kB,kD), 2X <--> X2
    end

    # Latexify.@generate_test latexify(rn; form=:ode)
    @test_broken latexify(rn; form = :ode) == replace(
    raw"$\begin{align}
    \frac{\mathrm{d} X\left( t \right)}{\mathrm{d}t} =& p - \left( X\left( t \right) \right)^{2} kB - d X\left( t \right) + 2 kD \mathrm{X2}\left( t \right) \\
    \frac{\mathrm{d} \mathrm{X2}\left( t \right)}{\mathrm{d}t} =& \frac{1}{2} \left( X\left( t \right) \right)^{2} kB - kD \mathrm{X2}\left( t \right)
    \end{align}
    $", "\r\n"=>"\n")

    # Currently latexify doesn't handle SDE systems properly, and they look identical to ode systems.
    # The "==" shoudl be a "!=", but due to latexify tests not working, for the broken test to work, I changed it.
    @test_broken latexify(rn; form=:sde) == replace(
    raw"$\begin{align}
    \frac{\mathrm{d} X\left( t \right)}{\mathrm{d}t} =& p - \left( X\left( t \right) \right)^{2} kB - d X\left( t \right) + 2 kD \mathrm{X2}\left( t \right) \\
    \frac{\mathrm{d} \mathrm{X2}\left( t \right)}{\mathrm{d}t} =& \frac{1}{2} \left( X\left( t \right) \right)^{2} kB - kD \mathrm{X2}\left( t \right)
    \end{align}
    $", "\r\n"=>"\n")

    # Tests that erroneous form gives error.
    @test_throws ErrorException latexify(rn; form=:xxx)
end


### Checks Reaction Network - Equations Combination ###

let
    base_network = @reaction_network begin
        k*r, X --> 0
    end
    @variables t r(t)
    @named decaying_rate = NonlinearSystem([r ~ -1], [r], [])
    extended = extend(decaying_rate, base_network)

    # Latexify.@generate_test latexify(extended)
    @test_broken latexify(extended) == replace(
    raw"\begin{align*}
    \mathrm{X} &\xrightarrow{k r} \varnothing  
    0 &= -1 - x\left( t \right)  
    \end{align*}
    ", "\r\n"=>"\n")

    # Latexify.@generate_test latexify(extended, mathjax=false)
    @test_broken latexify(extended, mathjax = false) == replace(
    raw"\begin{align*}
    \mathrm{X} &\xrightarrow{k r} \varnothing  
    0 &= -1 - x\left( t \right)  
    \end{align*}
    ", "\r\n"=>"\n")
end
