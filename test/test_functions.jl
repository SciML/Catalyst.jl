# Various functions that are useful for running the tests, and used across several test sets.

# Fetches the (required) Random package.
using Random, Test

### Initial Condition/Parameter Generators ###

# Generates a random initial condition (in the form of a map). Each value is a Float64.
function rnd_u0(sys, rng; factor = 1.0, min = 0.0)
    return [u => (min .+ factor .* rand(rng, size(u)...)) for u in unknowns(sys)]
end

# Generates a random initial condition (in the form of a map). Each value is a Int64.
function rnd_u0_Int64(sys, rng; n = 5, min = 0)
    return [u => (min .+ rand(rng, 1:n, size(u)...)) for u in unknowns(sys)]
end

# Generates a random parameter set (in the form of a map). Each value is a Float64.
function rnd_ps(sys, rng; factor = 1.0, min = 0.0)
    return [p => (min .+ factor .* rand(rng, size(p)...)) for p in parameters(sys)]
end

# Generates a random parameter set (in the form of a map). Each value is a Float64.
function rnd_ps_Int64(sys, rng; n = 5, min = 0)
    return [p => (min .+ rand(rng, 1:n, size(p)...)) for p in parameters(sys)]
end

# Used to convert a generated initial condition/parameter set to a vector that can be used for normal
# DiffEq functions (that are created for manual comparison). Requires order list of symbols.
function map_to_vec(map, syms)
    syms_dict = Dict([ModelingToolkit.getname(entry[1]) => entry[2] for entry in map])
    issetequal(keys(syms_dict), syms) || error("Map symbols ($(keys(syms_dict))) and symbol vector symbols ($(syms)) do not match.")
    return [syms_dict[sym] for sym in syms]
end

### System function evaluation ###

# Evaluates the the drift function of the ODE corresponding to a reaction network.
# Also checks that in place and out of place evaluations are identical.
function f_eval(rs::ReactionSystem, u, p, t; combinatoric_ratelaws = true)
    prob = ODEProblem(rs, u, 0.0, p; combinatoric_ratelaws)
    du = zeros(length(u))
    prob.f(du, prob.u0, prob.p, t)
    @test du == prob.f(prob.u0, prob.p, t)
    return du
end

# Evaluates the the Jacobian of the drift function of the ODE corresponding to a reaction network.
# Also checks that in place and out of place evaluations are identical.
function jac_eval(rs::ReactionSystem, u, p, t; combinatoric_ratelaws = true, sparse = false)
    prob = ODEProblem(rs, u, 0.0, p; jac = true, combinatoric_ratelaws, sparse)
    J = sparse ? deepcopy(prob.f.jac_prototype) : zeros(length(u), length(u))
    prob.f.jac(J, prob.u0, prob.p, t)
    @test J ≈ prob.f.jac(prob.u0, prob.p, t) atol = 1.0e-14 rtol = 1.0e-14
    return J
end

# Evaluates the the diffusion function of the SDE corresponding to a reaction network.
# Also checks that in place and out of place evaluations are identical.
function g_eval(rs::ReactionSystem, u, p, t; combinatoric_ratelaws = true)
    prob = SDEProblem(rs, u, 0.0, p; combinatoric_ratelaws)
    dW = zeros(length(u), numreactions(rs))
    prob.g(dW, prob.u0, prob.p, t)
    @test dW == prob.g(prob.u0, prob.p, t)
    return dW
end
