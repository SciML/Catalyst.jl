#! format: off

### Preparations ###

# Fetch packages.
using Catalyst, Latexify, Test

# Fetch test networks.
include("../test_networks.jl")

############################
### CURRENTLY NOT ACTIVE ###
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

# Generally, for all latexify tests, the lines after `@test latexify(rn) == replace(` must
# start without any tabs, hence the somewhat weird formatting.

### Basic Tests ###

# Tests functions on basic network (1).
# Tests on network with special functions (hill etc.).
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

# Removal of tab required for correctness.
# Latexify.@generate_test latexify(rn; expand_functions = false)
@test latexify(rn; expand_functions = false) == replace(
raw"\begin{align*}
\varnothing &\xrightarrow{\mathrm{hillr}\left( X2, v1, K1, n1 \right) \mathrm{hill}\left( X4, v1, K1, n1 \right)} \mathrm{X1} \\
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

# Latexify.@generate_test latexify(rn; expand_functions = true)
@test latexify(rn; expand_functions = true) == replace(
raw"\begin{align*}
\varnothing &\xrightarrow{\frac{X4^{n1} K1^{n1} v1^{2}}{\left( X2^{n1} + K1^{n1} \right) \left( K1^{n1} + X4^{n1} \right)}} \mathrm{X1} \\
\varnothing &\xrightarrow{\frac{X5^{n2} v2}{X5^{n2} + K2^{n2}}} \mathrm{X2} \\
\varnothing &\xrightarrow{\frac{X3^{n3} v3}{X3^{n3} + K3^{n3}}} \mathrm{X3} \\
\varnothing &\xrightarrow{\frac{K4^{n4} v4}{K4^{n4} + X1^{n4}}} \mathrm{X4} \\
\varnothing &\xrightarrow{\frac{X2^{n5} v5}{X2^{n5} + K5^{n5}}} \mathrm{X5} \\
\varnothing &\xrightarrow{\frac{X1^{n6} v6}{X1^{n6} + K6^{n6} + X6^{n6}}} \mathrm{X6} \\
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

# Latexify.@generate_test latexify(rn, mathjax = false)
@test latexify(rn, mathjax = false) == replace(
raw"\begin{align*}
\varnothing &\xrightarrow{\frac{X4^{n1} K1^{n1} v1^{2}}{\left( X2^{n1} + K1^{n1} \right) \left( K1^{n1} + X4^{n1} \right)}} \mathrm{X1} \\
\varnothing &\xrightarrow{\frac{X5^{n2} v2}{X5^{n2} + K2^{n2}}} \mathrm{X2} \\
\varnothing &\xrightarrow{\frac{X3^{n3} v3}{X3^{n3} + K3^{n3}}} \mathrm{X3} \\
\varnothing &\xrightarrow{\frac{K4^{n4} v4}{K4^{n4} + X1^{n4}}} \mathrm{X4} \\
\varnothing &\xrightarrow{\frac{X2^{n5} v5}{X2^{n5} + K5^{n5}}} \mathrm{X5} \\
\varnothing &\xrightarrow{\frac{X1^{n6} v6}{X1^{n6} + K6^{n6} + X6^{n6}}} \mathrm{X6} \\
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

# Tests basic functions on simple network (2).
let
    rn = @reaction_network begin
        (hill(B, p_a, k, n), d_a), 0 ↔ A
        (p_b, d_b), 0 ↔ B
        (r_a, r_b), 3B ↔ A
    end

# Latexify.@generate_test latexify(rn)
@test latexify(rn) == replace(
raw"\begin{align*}
\varnothing &\xrightleftharpoons[d_{a}]{\frac{B^{n} p_{a}}{B^{n} + k^{n}}} \mathrm{A} \\
\varnothing &\xrightleftharpoons[d_{b}]{p_{b}} \mathrm{B} \\
3 \mathrm{B} &\xrightleftharpoons[r_{b}]{r_{a}} \mathrm{A}  
 \end{align*}
", "\r\n"=>"\n")

# Latexify.@generate_test latexify(rn, mathjax = false)
@test latexify(rn, mathjax = false) == replace(
raw"\begin{align*}
\varnothing &\xrightleftharpoons[d_{a}]{\frac{B^{n} p_{a}}{B^{n} + k^{n}}} \mathrm{A} \\
\varnothing &\xrightleftharpoons[d_{b}]{p_{b}} \mathrm{B} \\
3 \mathrm{B} &\xrightleftharpoons[r_{b}]{r_{a}} \mathrm{A}  
 \end{align*}
", "\r\n"=>"\n")
end

# Tests for system with parametric stoichiometry.
let
    rn = @reaction_network begin
        p, 0 --> (m + n)*X
    end

# Latexify.@generate_test latexify(rn)
@test latexify(rn) == replace(
raw"\begin{align*}
\varnothing &\xrightarrow{p} (m + n)\mathrm{X}  
 \end{align*}
", "\r\n"=>"\n")
end

# Checks for systems with vector species/parameters.
# Technically tests would work, however, the display is non-ideal (https://github.com/SciML/Catalyst.jl/issues/932, https://github.com/JuliaSymbolics/Symbolics.jl/issues/1167).
let
    rn = @reaction_network begin
        @parameters k[1:2] x[1:2] [isconstantspecies=true]
        @species (X(t))[1:2] (K(t))[1:2]
        (k[1]*K[1],k[2]*K[2]), X[1] + x[1] <--> X[2] + x[2]
    end

    # Latexify.@generate_test latexify(rn)
    @test_broken false
end

### Tests the `form` Option ###

# Check option work for a large number of systems (and do not error).
let
    for rn in reaction_networks_standard
        @test latexify(rn)==latexify(rn; form=:reactions)
        #@test_broken latexify(ode_model(rn)) == latexify(rn; form=:ode) # Slight difference due to some latexify weirdly. Both displays fine though
    end
end

# Check for specific network.
let
    rn = @reaction_network begin
        (p,d), 0 <--> X
        (kB,kD), 2X <--> X2
    end

    # Note: The two trailing spaces on the last equation line are significant (from chemical_arrows).
    # Latexify.@generate_test latexify(rn; form=:ode)
    ode_expected = "\\begin{align*}\n\\frac{\\mathrm{d} X\\left( t \\right)}{\\mathrm{d}t} &= p - d X\\left( t \\right) + 2 kD \\mathrm{X2}\\left( t \\right) - \\left( X\\left( t \\right) \\right)^{2} kB \\\\\n\\frac{\\mathrm{d} \\mathrm{X2}\\left( t \\right)}{\\mathrm{d}t} &=  - kD \\mathrm{X2}\\left( t \\right) + \\frac{1}{2} \\left( X\\left( t \\right) \\right)^{2} kB  \n \\end{align*}\n"
    @test latexify(rn; form = :ode) == replace(ode_expected, "\r\n"=>"\n")

    # SDE form currently renders identically to ODE form (noise terms not yet shown).
    @test latexify(rn; form = :sde) == latexify(rn; form = :ode)

    # Tests that erroneous form gives error.
    @test_throws ErrorException latexify(rn; form=:xxx)
end

### Other Tests ###

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

# Test on an empty system.
let
    empty_rn = ReactionSystem(Reaction[]; name=:EmptySys)

    # Latexify.@generate_test latexify(empty_rn)
    @test latexify(empty_rn) == replace(
    raw"ReactionSystem EmptySys has no reactions or equations.", "\r\n"=>"\n")
end

# Test for https://github.com/SciML/Catalyst.jl/issues/473.
# (Error where there's an equation in the reaction rate)
let
    rn = @reaction_network begin
        k*Y, Y --> ∅
    end

 # Latexify.@generate_test latexify(rn)
@test latexify(rn) == replace(
raw"\begin{align*}
\mathrm{Y} &\xrightarrow{Y k} \varnothing  
 \end{align*}
", "\r\n"=>"\n")
end

# Test with reaction system containing a differential equation.
let
    # Declare model.
    rn = @reaction_network begin
        @equations D(V) ~ X - V
        (p/V,d/V), 0 <--> X
    end

# Latexify.@generate_test latexify(rn)
@test latexify(rn) == replace(
raw"\begin{align*}
\varnothing &\xrightleftharpoons[\frac{d}{V\left( t \right)}]{\frac{p}{V\left( t \right)}} \mathrm{X} \\
\frac{\mathrm{d} V\left( t \right)}{\mathrm{d}t} &=  - V\left( t \right) + X\left( t \right)  
 \end{align*}
", "\r\n"=>"\n")
end
