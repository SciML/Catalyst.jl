### Preparations ###

# Fetch packages.
using Catalyst, Latexify, Test

# These tests check the structure of the `:sde` and `:jump` forms: which equations there are,
# and which noise (dWⱼ) and Poisson process (Rⱼ) terms each has, with what sign and
# coefficient. They don't compare whole strings, whose details depend on how Latexify and
# Symbolics print expressions (which made the string comparisons in latexify.jl fragile).

# The rows of a latexified `align*` environment.
function latex_rows(ltx)
    body = split(split(string(ltx), "\\begin{align*}\n")[2], "\\end{align*}")[1]
    return strip.(split(body, "\\\\\n"))
end

# The right-hand side of a row, split into its top-level summands, each starting with its sign
# ("+ " for the first one when it has none).
function summands(row)
    rhs = strip(split(row, "&= "; limit = 2)[2])
    terms = String[]
    depth = 0
    start = 1
    for (i, c) in pairs(rhs)
        if c in ('{', '(')
            depth += 1
        elseif c in ('}', ')')
            depth -= 1
        elseif depth == 0 && i > 1 && c in ('+', '-') && rhs[prevind(rhs, i)] == ' '
            push!(terms, strip(rhs[start:prevind(rhs, i)]))
            start = i
        end
    end
    push!(terms, strip(rhs[start:end]))
    return [startswith(t, r"[+-] ") ? t : "+ " * t for t in terms if !isempty(t)]
end

# The summand of a row containing `marker` (e.g. "\mathrm{d}W_{2}"), or `nothing`.
function term_with(row, marker)
    idx = findall(t -> occursin(marker, t), summands(row))
    return isempty(idx) ? nothing : summands(row)[only(idx)]
end

dW(j) = "\\mathrm{d}W_{$j}"
R(j) = "R_{$j}\\left( \\int_{0}^{t}"

### Tests ###

# The chemical Langevin equation of a birth-death process with dimerization.
let
    rn = @reaction_network begin
        (p, d), 0 <--> X
        k, 2X --> X2
    end
    rows = latex_rows(latexify(rn; form = :sde))
    @test length(rows) == 2
    @test all(row -> startswith(row, "\\mathrm{d}"), rows)

    # X: drift, +√p dW₁, -√(dX) dW₂, -2√(kX²/2) dW₃.
    @test term_with(rows[1], "\\mathrm{d}t") !== nothing
    @test startswith(term_with(rows[1], dW(1)), "+ \\sqrt{")
    @test startswith(term_with(rows[1], dW(2)), "- \\sqrt{")
    @test startswith(term_with(rows[1], dW(3)), "- 2 \\sqrt{")
    @test length(summands(rows[1])) == 4

    # X2: drift, +√(kX²/2) dW₃.
    @test term_with(rows[2], "\\mathrm{d}t") !== nothing
    @test startswith(term_with(rows[2], dW(3)), "+ \\sqrt{")
    @test term_with(rows[2], dW(1)) === nothing
    @test length(summands(rows[2])) == 2

    # Not the ODE form, and the time arguments only with `show_time_arg = true`.
    @test latexify(rn; form = :sde) != latexify(rn; form = :ode)
    @test !occursin("\\left( t \\right)", string(latexify(rn; form = :sde)))
    rows_t = latex_rows(latexify(rn; form = :sde, show_time_arg = true))
    @test occursin("\\mathrm{d}W_{1}\\left( t \\right)", rows_t[1])
    @test startswith(rows_t[1], "\\mathrm{d}X\\left( t \\right)")
end

# The random time change representation of the same network.
let
    rn = @reaction_network begin
        (p, d), 0 <--> X
        k, 2X --> X2
    end
    rows = latex_rows(latexify(rn; form = :jump))
    @test length(rows) == 2
    @test startswith(rows[1], "X\\left( t \\right) &= ")

    # X(t) = X(0) + R₁(…) - R₂(…) - 2R₃(…).
    @test summands(rows[1])[1] == "+ X\\left( 0 \\right)"
    @test startswith(term_with(rows[1], R(1)), "+ R_{1}")
    @test startswith(term_with(rows[1], R(2)), "- R_{2}")
    @test startswith(term_with(rows[1], R(3)), "- 2 R_{3}")
    @test length(summands(rows[1])) == 4

    # X2(t) = X2(0) + R₃(…).
    @test startswith(term_with(rows[2], R(3)), "+ R_{3}")
    @test length(summands(rows[2])) == 2

    # The integrands are in s, with the species at s.
    @test occursin("X\\left( s \\right)", term_with(rows[1], R(2)))
    @test occursin("\\, \\mathrm{d}s \\right)", term_with(rows[1], R(2)))
    @test !occursin("X\\left( t \\right)", term_with(rows[1], R(2)))
end

# Symbolic stoichiometry, a time-dependent rate, noise scaling, and a parameter named `s`.
let
    rn = @reaction_network begin
        @parameters η
        k * t, m * A --> 2B
        b, B --> 0, [noise_scaling = η]
        s, 0 --> A
    end
    rows = latex_rows(latexify(rn; form = :sde))
    @test startswith(term_with(rows[1], dW(1)), "- m \\sqrt{")
    @test startswith(term_with(rows[2], dW(1)), "+ 2 \\sqrt{")
    @test occursin("\\eta", term_with(rows[2], dW(2)))

    # With a parameter named `s`, the integration variable is τ.
    rows = latex_rows(latexify(rn; form = :jump))
    @test startswith(term_with(rows[1], R(1)), "- m R_{1}")
    @test occursin("\\mathrm{d}\\tau", term_with(rows[1], R(1)))
    @test occursin("A\\left( \\tau \\right)", term_with(rows[1], R(1)))
end

# A coupled differential equation with a brownian: shown in the SDE form with its noise term;
# the jump form can't display it.
let
    rn = @reaction_network begin
        @parameters σ
        @brownians W
        @variables V(t)
        @equations D(V) ~ -V + σ * W
        (p, d), 0 <--> X
    end
    rows = latex_rows(latexify(rn; form = :sde))
    @test length(rows) == 2
    vrow = only(filter(row -> startswith(row, "\\mathrm{d}V"), rows))
    @test term_with(vrow, "\\mathrm{d}t") !== nothing
    @test occursin("\\sigma", term_with(vrow, "\\mathrm{d}W"))
    @test length(summands(vrow)) == 2
    @test_throws ErrorException latexify(rn; form = :jump)
end

# Options shared with the other forms, an empty system, and an unknown form.
let
    rn = @reaction_network begin
        (p, d), 0 <--> X
    end
    for form in (:sde, :jump)
        @test occursin("\\begin{align}", string(latexify(rn; form, starred = false)))
        delimited = string(latexify(rn; form, math_delimiters = true))
        @test startswith(delimited, "\$\$\n") && endswith(delimited, "\$\$")
    end
    empty_rn = @reaction_network begin end
    @test occursin("has no reactions or equations", string(latexify(empty_rn; form = :sde)))
    @test occursin("has no reactions or equations", string(latexify(empty_rn; form = :jump)))
    @test_throws ErrorException latexify(rn; form = :xxx)
end
