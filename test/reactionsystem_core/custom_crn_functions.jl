### Prepares Tests ###

# Fetch packages.
using Catalyst, Test, LinearAlgebra
using ModelingToolkit: get_continuous_events, get_discrete_events
using Symbolics: derivative

# Sets stable rng number.
using StableRNGs
rng = StableRNG(12345)

# Sets the default `t` to use.
t = default_t()

# Fetch test functions.
include("../test_functions.jl")


### Basic Tests ###

# Compares network written with and without special functions.
let
    new_hill(x, v, k, n) = v * x^n / (k^n + x^n)
    new_poly(x, p1, p2) = 0.5 * p1 * x^2
    new_exp(x, p) = exp(-p * x)

    custom_function_network_1 = @reaction_network begin
        hill(X1, v1, K1, 2), X1 + Y1 --> Z1
        mm(X2, v2, K2), X2 + Y2 --> Z2
        p1 * X3^2 + p2, X3 + Y3 --> Z3
        exp(-p3 * Y4), X4 + Y4 --> Z4
        hillr(X5, v3, K3, 2), X5 + Y5 --> Z5
        mmr(X6, v4, K4), X6 + Y6 --> Z6
        hillar(X7, Y7, v5, K5, 2), X7 + Y7 --> Z7
    end

    custom_function_network_2 = @reaction_network begin
        new_hill(X1, v1, K1, 2), X1 + Y1 --> Z1
        v2 * X2 / (X2 + K2), X2 + Y2 --> Z2
        2 * new_poly(X3, p1, p2) + p2, X3 + Y3 --> Z3
        new_exp(Y4, p3), X4 + Y4 --> Z4
        v3 * (K3^2) / (K3^2 + X5^2), X5 + Y5 --> Z5
        v4 * K4 / (X6 + K4), X6 + Y6 --> Z6
        v5 * (X7^2) / (K5^2 + X7^2 + Y7^2), X7 + Y7 --> Z7
    end

    for factor in [1e-2, 1e-1, 1e0, 1e1, 1e2]
        u0 = rnd_u0(custom_function_network_1, rng; factor)
        ps = rnd_ps(custom_function_network_1, rng; factor)
        t = rand(rng)
        @test f_eval(custom_function_network_1, u0, ps, t) ≈ f_eval(custom_function_network_2, u0, ps, t)
        @test jac_eval(custom_function_network_1, u0, ps, t) ≈ jac_eval(custom_function_network_2, u0, ps, t)
        @test g_eval(custom_function_network_1, u0, ps, t) ≈ g_eval(custom_function_network_2, u0, ps, t)
    end
end

# Compares the symbolic derivatives with their manually computed forms.
let
    @variables X Y
    @parameters v K n

    @test isequal(derivative(Catalyst.mm(X, v, K), X), v * K / (K + X)^2)
    @test isequal(derivative(Catalyst.mm(X, v, K), v), X / (K + X))
    @test isequal(derivative(Catalyst.mm(X, v, K), K), -v * X / (K + X)^2)

    @test isequal(derivative(Catalyst.mmr(X, v, K), X), -v * K / (K + X)^2)
    @test isequal(derivative(Catalyst.mmr(X, v, K), v), K / (K + X))
    @test isequal(derivative(Catalyst.mmr(X, v, K), K), v * X / (K + X)^2)

    @test isequal(derivative(Catalyst.hill(X, v, K, n), X),
                  n * v * (K^n) * (X^(n - 1)) / (K^n + X^n)^2)
    @test isequal(derivative(Catalyst.hill(X, v, K, n), v), X^n / (K^n + X^n))
    @test isequal(derivative(Catalyst.hill(X, v, K, n), K),
                  -n * v * (K^(n - 1)) * (X^n) / (K^n + X^n)^2)
    @test isequal(derivative(Catalyst.hill(X, v, K, n), n),
                  v * (X^n) * (K^n) * (log(X) - log(K)) / (K^n + X^n)^2)

    @test isequal(derivative(Catalyst.hillr(X, v, K, n), X),
                  -n * v * (K^n) * (X^(n - 1)) / (K^n + X^n)^2)
    @test isequal(derivative(Catalyst.hillr(X, v, K, n), v), K^n / (K^n + X^n))
    @test isequal(derivative(Catalyst.hillr(X, v, K, n), K),
                  n * v * (K^(n - 1)) * (X^n) / (K^n + X^n)^2)
    @test isequal(derivative(Catalyst.hillr(X, v, K, n), n),
                  v * (X^n) * (K^n) * (log(K) - log(X)) / (K^n + X^n)^2)

    @test isequal(derivative(Catalyst.hillar(X, Y, v, K, n), X),
                  n * v * (K^n + Y^n) * (X^(n - 1)) / (K^n + X^n + Y^n)^2)
    @test isequal(derivative(Catalyst.hillar(X, Y, v, K, n), Y),
                  -n * v * (Y^(n - 1)) * (X^n) / (K^n + X^n + Y^n)^2)
    @test isequal(derivative(Catalyst.hillar(X, Y, v, K, n), v),
                  X^n / (K^n + X^n + Y^n))
    @test isequal(derivative(Catalyst.hillar(X, Y, v, K, n), K),
                  -n * v * (v^(n - 1)) * (X^n) / (K^n + X^n + Y^n)^2)
    @test isequal(derivative(Catalyst.hillar(X, Y, v, K, n), n),
                  v * (X^n) * ((K^n + Y^n) * log(X) - (K^n) * log(K) - (Y^n) * log(Y)) /
                  (K^n + X^n + Y^n)^2)
end

### Tests Function Expansion ###

# Tests `ReactionSystem`s.
let
    @species x(t) y(t)
    @parameters k v n
    rs1 = @reaction_network rs begin
        mm(x, v, k), 0 --> x
        mmr(x, v, k), 0 --> x
        hill(x, v, k, n), 0 --> x
        hillr(x, v, k, n), 0 --> x
        hillar(x, y, v, k, n), 0 --> x
    end
    rs2 = @reaction_network rs begin
        v * x / (x + k), 0 --> x
        v * k / (x + k), 0 --> x
        v * (x^n) / (x^n + k^n), 0 --> x
        v * (k^n) / (x^n + k^n), 0 --> x
        v * (x^n) / (x^n + y^n + k^n), 0 --> x
    end
    rs3 = Catalyst.expand_registered_functions(rs1)
    @test isequivalent(rs3, rs2; ignorenames = false, debug = true)
end

# Tests `Reaction`s.
let
    @species x(t) y(t)
    @parameters k v n

    r1 = @reaction mm(x, v, k), 0 --> x
    r2 = @reaction mmr(x, v, k), 0 --> x
    r3 = @reaction hill(x, v, k, n), 0 --> x
    r4 = @reaction hillr(x, v, k, n), 0 --> x
    r5 = @reaction hillar(x, y, v, k, n), 0 --> x + y

    @test isequal(Catalyst.expand_registered_functions(r1).rate, v * x / (x + k))
    @test isequal(Catalyst.expand_registered_functions(r2).rate, v * k / (x + k))
    @test isequal(Catalyst.expand_registered_functions(r3).rate, v * (x^n) / (x^n + k^n))
    @test isequal(Catalyst.expand_registered_functions(r4).rate, v * (k^n) / (x^n + k^n))
    @test isequal(Catalyst.expand_registered_functions(r5).rate, v * (x^n) / (x^n + y^n + k^n))
end

# Tests `Equation`s.
# Tests using non-default independent variable.
let
    @parameters T
    @variables X(T) Y(T)
    @parameters K V N

    eq1 = 0 ~ mm(X, V, K)
    eq2 = 0 ~ mmr(X, V, K)
    eq3 = 0 ~ hill(X, V, K, N)
    eq4 = 0 ~ hillr(X, V, K, N)
    eq5 = 0 ~ hillar(X, Y, V, K, N)

    @test isequal(Catalyst.expand_registered_functions(eq1), 0 ~ V * X / (X + K))
    @test isequal(Catalyst.expand_registered_functions(eq2), 0 ~ V * K / (X + K))
    @test isequal(Catalyst.expand_registered_functions(eq3), 0 ~ V * (X^N) / (X^N + K^N))
    @test isequal(Catalyst.expand_registered_functions(eq4), 0 ~ V * (K^N) / (X^N + K^N))
    @test isequal(Catalyst.expand_registered_functions(eq5), 0 ~ V * (X^N) / (X^N + Y^N + K^N))
end

# Ensures that original system is not modified.
let
    # Create model with a registered function.
    @species X(t)
    @variables V(t)
    @parameters v K
    eqs = [
        Reaction(mm(X,v,K), [], [X]),
        mm(V,v,K) ~ V + 1
    ]
    @named rs = ReactionSystem(eqs, t)

    # Check that `expand_registered_functions` does not mutate original model.
    rs_expanded_funcs = Catalyst.expand_registered_functions(rs)
    @test isequal(only(Catalyst.get_rxs(rs)).rate, Catalyst.mm(X,v,K))
    @test isequal(only(Catalyst.get_rxs(rs_expanded_funcs)).rate, v*X/(X + K))
    @test isequal(last(Catalyst.get_eqs(rs)).lhs, Catalyst.mm(V,v,K))
    @test isequal(last(Catalyst.get_eqs(rs_expanded_funcs)).lhs, v*V/(V + K))
end

# Tests on model with events.
let
    # Creates a model, saves it, and creates an expanded version.
    rs = @reaction_network begin
        @continuous_events begin
            [mm(X,v,K) ~ 1.0] => [X ~ X]
        end
        @discrete_events begin
            [1.0] => [X ~ mmr(X,v,K) + Y*(v + K)]
            1.0 => [X ~ X]
            (hill(X,v,K,n) > 1000.0) => [X ~ hillr(X,v,K,n) + 2]
        end
        v0 + hillar(X,Y,v,K,n), X --> Y
    end
    rs_saved = deepcopy(rs)
    rs_expanded = Catalyst.expand_registered_functions(rs)

    # Checks that the original model is unchanged (equality currently does not consider events).
    @test rs == rs_saved
    @test get_continuous_events(rs) == get_continuous_events(rs_saved)
    @test get_discrete_events(rs) == get_discrete_events(rs_saved)

    # Checks that the new system is expanded.
    @unpack v0, X, Y, v, K, n = rs
    continuous_events = [
        [v*X/(X + K) ~ 1.0] => [X ~ X]
    ]
    discrete_events = [
        [1.0] => [X ~ v*K/(X + K) + Y*(v + K)]
        1.0 => [X ~ X]
        (v * (X^n) / (X^n + K^n) > 1000.0) => [X ~ v * (K^n) / (X^n + K^n) + 2]
    ]
    continuous_events = ModelingToolkit.SymbolicContinuousCallback.(continuous_events)
    discrete_events = ModelingToolkit.SymbolicDiscreteCallback.(discrete_events)
    @test isequal(only(Catalyst.get_rxs(rs_expanded)).rate, v0 + v * (X^n) / (X^n + Y^n + K^n))
    @test isequal(get_continuous_events(rs_expanded), continuous_events)
    @test isequal(get_discrete_events(rs_expanded), discrete_events)
end

# test for hill function expansion
let
    rn = @reaction_network begin
        hill(X, v, K, n), X + Y --> Z
        mm(X, v, K), X + Y --> Z
        hillr(X, v, K, n), X + Y --> Z
        mmr(X, v, K), X + Y --> Z
    end
    osys = complete(make_rre_ode(rn; expand_catalyst_funs = false))
    t = default_t()
    D = default_time_deriv()
    @unpack X, v, K, n, Y, Z = rn
    osyseqs = equations(osys)
    eqs = [D(X) ~ -hill(X, v, K, n)*X*Y - mm(X,v,K)*X*Y - hillr(X,v,K,n)*X*Y - mmr(X,v,K)*X*Y,
           D(Y) ~ -hill(X, v, K, n)*X*Y - mm(X,v,K)*X*Y - hillr(X,v,K,n)*X*Y - mmr(X,v,K)*X*Y,
           D(Z) ~ hill(X, v, K, n)*X*Y + mm(X,v,K)*X*Y + hillr(X,v,K,n)*X*Y + mmr(X,v,K)*X*Y]
    reorder = [findfirst(eq -> isequal(eq.lhs, osyseq.lhs), eqs) for osyseq in osyseqs]
    for (osysidx,eqidx) in enumerate(reorder)
        @test iszero(simplify(eqs[eqidx].rhs - osyseqs[osysidx].rhs))
    end

    osys2 = complete(make_rre_ode(rn))
    hill2(x, v, k, n) = v * x^n / (k^n + x^n)
    mm2(X,v,K) = v*X / (X + K)
    mmr2(X,v,K) = v*K / (X + K)
    hillr2(X,v,K,n) = v * (K^n) / (X^n + K^n)
    eqs2 = [D(X) ~ -hill2(X, v, K, n)*X*Y - mm2(X,v,K)*X*Y - hillr2(X,v,K,n)*X*Y - mmr2(X,v,K)*X*Y,
        D(Y) ~ -hill2(X, v, K, n)*X*Y - mm2(X,v,K)*X*Y - hillr2(X,v,K,n)*X*Y - mmr2(X,v,K)*X*Y,
        D(Z) ~ hill2(X, v, K, n)*X*Y + mm2(X,v,K)*X*Y + hillr2(X,v,K,n)*X*Y + mmr2(X,v,K)*X*Y]
    osyseqs2 = equations(osys2)
    reorder = [findfirst(eq -> isequal(eq.lhs, osyseq.lhs), eqs2) for osyseq in osyseqs2]
    for (osysidx,eqidx) in enumerate(reorder)
        @test iszero(simplify(eqs2[eqidx].rhs - osyseqs2[osysidx].rhs))
    end

    nlsys = complete(make_rre_algeqs(rn; expand_catalyst_funs = false))
    nlsyseqs = equations(nlsys)
    eqs = [0 ~ -hill(X, v, K, n)*X*Y - mm(X,v,K)*X*Y - hillr(X,v,K,n)*X*Y - mmr(X,v,K)*X*Y,
           0 ~ -hill(X, v, K, n)*X*Y - mm(X,v,K)*X*Y - hillr(X,v,K,n)*X*Y - mmr(X,v,K)*X*Y,
           0 ~ hill(X, v, K, n)*X*Y + mm(X,v,K)*X*Y + hillr(X,v,K,n)*X*Y + mmr(X,v,K)*X*Y]
    for (i, eq) in enumerate(eqs)
        @test iszero(simplify(eq.rhs - nlsyseqs[i].rhs))
    end

    nlsys2 = complete(make_rre_algeqs(rn))
    nlsyseqs2 = equations(nlsys2)
    eqs2 = [0 ~ -hill2(X, v, K, n)*X*Y - mm2(X,v,K)*X*Y - hillr2(X,v,K,n)*X*Y - mmr2(X,v,K)*X*Y,
            0 ~ -hill2(X, v, K, n)*X*Y - mm2(X,v,K)*X*Y - hillr2(X,v,K,n)*X*Y - mmr2(X,v,K)*X*Y,
            0 ~ hill2(X, v, K, n)*X*Y + mm2(X,v,K)*X*Y + hillr2(X,v,K,n)*X*Y + mmr2(X,v,K)*X*Y]
    for (i, eq) in enumerate(eqs2)
        @test iszero(simplify(eq.rhs - nlsyseqs2[i].rhs))
    end

    sdesys = complete(make_cle_sde(rn; expand_catalyst_funs = false))
    sdesyseqs = equations(sdesys)
    eqs = [D(X) ~ -hill(X, v, K, n)*X*Y - mm(X,v,K)*X*Y - hillr(X,v,K,n)*X*Y - mmr(X,v,K)*X*Y,
           D(Y) ~ -hill(X, v, K, n)*X*Y - mm(X,v,K)*X*Y - hillr(X,v,K,n)*X*Y - mmr(X,v,K)*X*Y,
           D(Z) ~ hill(X, v, K, n)*X*Y + mm(X,v,K)*X*Y + hillr(X,v,K,n)*X*Y + mmr(X,v,K)*X*Y]
    reorder = [findfirst(eq -> isequal(eq.lhs, sdesyseq.lhs), eqs) for sdesyseq in sdesyseqs]
    for (sdesysidx,eqidx) in enumerate(reorder)
        @test iszero(simplify(eqs[eqidx].rhs - sdesyseqs[sdesysidx].rhs))
    end
    sdesysnoiseeqs = ModelingToolkit.get_noiseeqs(sdesys)
    neqvec = diagm(sqrt.(abs.([hill(X, v, K, n)*X*Y, mm(X,v,K)*X*Y, hillr(X,v,K,n)*X*Y, mmr(X,v,K)*X*Y])))
    neqmat = [-1 -1 -1 -1; -1 -1 -1 -1; 1 1 1 1]
    neqmat *= neqvec
    @test all(iszero, simplify.(sdesysnoiseeqs .- neqmat))

    sdesys = complete(make_cle_sde(rn))
    sdesyseqs = equations(sdesys)
    eqs = [D(X) ~ -hill2(X, v, K, n)*X*Y - mm2(X,v,K)*X*Y - hillr2(X,v,K,n)*X*Y - mmr2(X,v,K)*X*Y,
           D(Y) ~ -hill2(X, v, K, n)*X*Y - mm2(X,v,K)*X*Y - hillr2(X,v,K,n)*X*Y - mmr2(X,v,K)*X*Y,
           D(Z) ~ hill2(X, v, K, n)*X*Y + mm2(X,v,K)*X*Y + hillr2(X,v,K,n)*X*Y + mmr2(X,v,K)*X*Y]
    reorder = [findfirst(eq -> isequal(eq.lhs, sdesyseq.lhs), eqs) for sdesyseq in sdesyseqs]
    for (sdesysidx,eqidx) in enumerate(reorder)
        @test iszero(simplify(eqs[eqidx].rhs - sdesyseqs[sdesysidx].rhs))
    end
    sdesysnoiseeqs = ModelingToolkit.get_noiseeqs(sdesys)
    neqvec = diagm(sqrt.(abs.([hill2(X, v, K, n)*X*Y, mm2(X,v,K)*X*Y, hillr2(X,v,K,n)*X*Y, mmr2(X,v,K)*X*Y])))
    neqmat = [-1 -1 -1 -1; -1 -1 -1 -1; 1 1 1 1]
    neqmat *= neqvec
    @test all(iszero, simplify.(sdesysnoiseeqs .- neqmat))

    jsys = make_sck_jump(rn; expand_catalyst_funs = false)
    jsyseqs = equations(jsys).x[2]
    rates = getfield.(jsyseqs, :rate)
    affects = getfield.(jsyseqs, :affect!)
    reqs = [ Y*X*hill(X, v, K, n), Y*X*mm(X, v, K), hillr(X, v, K, n)*Y*X, Y*X*mmr(X, v, K)]
    affeqs = [Z ~ 1 + Z, Y ~ -1 + Y, X ~ -1 + X]
    @test all(iszero, simplify(rates .- reqs))
    @test all(aff -> isequal(aff, affeqs), affects)

    jsys = make_sck_jump(rn)
    jsyseqs = equations(jsys).x[2]
    rates = getfield.(jsyseqs, :rate)
    affects = getfield.(jsyseqs, :affect!)
    reqs = [ Y*X*hill2(X, v, K, n), Y*X*mm2(X, v, K), hillr2(X, v, K, n)*Y*X, Y*X*mmr2(X, v, K)]
    affeqs = [Z ~ 1 + Z, Y ~ -1 + Y, X ~ -1 + X]
    @test all(iszero, simplify(rates .- reqs))
    @test all(aff -> isequal(aff, affeqs), affects)
end
