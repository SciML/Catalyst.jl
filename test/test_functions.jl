# Various functions that are useful for running the tests, and used across several test sets.


# Initial Condition/Parameter Generators ###

# Generates a random initial condition (in the form of a map). Each value is a Float64.
function rnd_u0(sys, rng; factor = 1.0, min = 0.0)
    return [u => min + factor * rand(rng) for u in unknowns(sys)]
end

# Generates a random initial condition (in the form of a map). Each value is a Int64.
function rnd_u0_Int64(sys, rng; n = 5, min = 0)
    return [u => min + rand(rng, 1:n) for u in unknowns(sys)]
end

# Generates a random parameter set (in the form of a map). Each value is a Float64.
function rnd_ps(sys, rng; factor = 1.0, min = 0.0)
    return [p => min + factor * rand(rng) for p in parameters(sys)]
end

# Generates a random parameter set (in the form of a map). Each value is a Float64.
function rnd_ps_Int64(sys, rng; n = 5, min = 0)
    return [p => min + rand(rng, 1:n) for p in parameters(sys)]
end

### System function evaluation ###

# Evaluates the the drift function of the ODE corresponding to a reaction network.
function f_eval(rs::ReactionSystem, u, p, τ)
    prob = ODEProblem(rs, u, (0.0, 0.0), p)
    return prob.f(prob.u0, prob.p, τ)
end

# Evaluates the the Jacobian of the drift function of the ODE corresponding to a reaction network.
function jac_eval(rs::ReactionSystem, u, p, t)
    prob = ODEProblem(rs, u, (0.0, 0.0), p; jac = true)
    return prob.f.jac(prob.u0, prob.p, t)
end

# Evaluates the the diffusion function of the SDE corresponding to a reaction network.
function g_eval(rs::ReactionSystem, u, p, t)
    prob = SDEProblem(rs, u, (0.0, 0.0), p)
    return prob.g(prob.u0, prob.p, t)
end
