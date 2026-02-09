#! format: off

### Prepares Tests ###

# Fetch packages.
using Catalyst, DynamicQuantities, DataInterpolations, Test
using Catalyst: catalyst_get_unit, SYM_UNITLESS, validate_units, assert_valid_units, UnitValidationError
using ModelingToolkitBase: get_iv

### Public API surface tests ###

let
    @test !Base.isexported(Catalyst, :validate_units)
    @test !Base.isexported(Catalyst, :assert_valid_units)
    @test !Base.isexported(Catalyst, :unit_validation_report)
    @test !Base.isexported(Catalyst, :UnitValidationError)
    @test !Base.isexported(Catalyst, :UnitValidationIssue)
    @test !Base.isexported(Catalyst, :UnitValidationReport)

    @static if VERSION >= v"1.11"
        @test Base.ispublic(Catalyst, :validate_units)
        @test Base.ispublic(Catalyst, :assert_valid_units)
        @test Base.ispublic(Catalyst, :unit_validation_report)
        @test Base.ispublic(Catalyst, :UnitValidationError)
        @test Base.ispublic(Catalyst, :UnitValidationIssue)
        @test Base.ispublic(Catalyst, :UnitValidationReport)
    end
end


### catalyst_get_unit Tests ###

# Tests basic unit access on individual variables and constants.
let
    @independent_variables t [unit=us"s"]
    @species A(t) [unit=us"M"] B(t) [unit=us"μM"]
    @parameters k1 [unit=us"s^(-1)"] k2 [unit=us"M^(-1)*s^(-1)"] γ

    # Species and parameter units.
    @test catalyst_get_unit(A) == us"M"
    @test catalyst_get_unit(B) == us"μM"
    @test catalyst_get_unit(k1) == us"s^(-1)"
    @test catalyst_get_unit(k2) == us"M^(-1)*s^(-1)"
    @test catalyst_get_unit(t) == us"s"

    # Unitless parameter (no [unit=...] metadata).
    @test catalyst_get_unit(γ) == SYM_UNITLESS

    # Numeric constants.
    @test catalyst_get_unit(1.0) == SYM_UNITLESS
    @test catalyst_get_unit(42) == SYM_UNITLESS
end

# Tests that concrete Dimensions (u"...") are rejected in favor of SymbolicDimensions (us"...").
let
    @independent_variables t [unit=us"s"]
    @species A(t) [unit=u"M"]
    @test_throws ArgumentError catalyst_get_unit(A)
end

# Tests catalyst_get_unit on expression types.
let
    @independent_variables t [unit=us"s"]
    @species X(t) [unit=us"M"]
    @parameters k [unit=us"s^(-1)"] v [unit=us"M/s"] K [unit=us"M"]

    # Multiplication.
    @test catalyst_get_unit(k * X) == us"M/s"

    # Division.
    @test catalyst_get_unit(X / K) == SYM_UNITLESS

    # Power.
    @test catalyst_get_unit(X^2) == us"M^2"

    # Compound expression.
    @test catalyst_get_unit(k * X^2) == us"M^2/s"

    # Addition (returns first term's unit).
    @test catalyst_get_unit(k + k) == us"s^(-1)"

    # Subtraction (represented as addition internally).
    @test catalyst_get_unit(v - k * X) == us"M/s"

    # Nested mixed arithmetic: (v + k*X) * K / X^2 = (M/s)*M/M^2 = 1/s.
    @test catalyst_get_unit((v + k * X) * K / X^2) == us"s^(-1)"

    # Differential.
    D = Differential(t)
    @test catalyst_get_unit(D(X)) == us"M/s"
end

# Tests catalyst_get_unit with registered functions.
let
    @independent_variables t [unit=us"s"]
    @species X(t) [unit=us"M"]
    @parameters v [unit=us"M/s"] K [unit=us"M"] n

    # mm(X, v, K) = v * X / (X + K), returns units of v.
    @test catalyst_get_unit(mm(X, v, K)) == us"M/s"

    # mmr(X, v, K) = v * K / (X + K), returns units of v.
    @test catalyst_get_unit(mmr(X, v, K)) == us"M/s"

    # hill(X, v, K, n) = v * X^n / (X^n + K^n), returns units of v.
    @test catalyst_get_unit(hill(X, v, K, n)) == us"M/s"

    # Registered vs traced form should match.
    mm_traced(X, v, K) = v * X / (X + K)
    @test catalyst_get_unit(mm(X, v, K)) == catalyst_get_unit(mm_traced(X, v, K))
end


### Basic Programmatic Tests (Molar Units) ###

# Tests that units work with programmatic model creation using M units.
let
    @independent_variables t [unit=us"s"]
    @species A(t) [unit=us"M"] B(t) [unit=us"M"] C(t) [unit=us"M"]
    @parameters k1 [unit=us"M/s"] k2 [unit=us"s^(-1)"] k3 [unit=us"M^(-1)*s^(-1)"]
    rxs = [Reaction(k1, nothing, [A]),
        Reaction(k2, [A], [B]),
        Reaction(k3, [A, B], [B], [1, 1], [2])]
    @test_nowarn ReactionSystem(rxs, t, [A, B, C], [k1, k2, k3]; name = :rs)
    @named rs = ReactionSystem(rxs, t, [A, B, C], [k1, k2, k3])
    rs = complete(rs)

    # Tests validate passes.
    @test validate_units(rs)

    # Tests that all reactions have the correct unit.
    for rx in reactions(rs)
        @test catalyst_get_unit(oderatelaw(rx)) == us"M/s"
        # We don't currently convert units, so they will be the same as for ODEs.
        @test catalyst_get_unit(jumpratelaw(rx)) == us"M/s"
    end

    # Tests that the system can be converted to models without warnings.
    @test_nowarn ode_model(rs)
    @test_nowarn sde_model(rs)
    @test_nowarn jump_model(rs)
    @test_nowarn ss_ode_model(rs)

    # Tests that creating Reactions with non-matching species units yields warnings.
    @species Bmol(t) [unit=us"mol"] D(t) [unit=us"kg"]
    bad_rx1 = Reaction(k1, [A], [D])
    bad_rx2 = Reaction(k1, [A], [Bmol, D])
    bad_rx3 = Reaction(k1, [A, D], [Bmol])
    @test (@test_logs (:warn, ) match_mode=:any validate_units(bad_rx1)) == false
    @test (@test_logs (:warn, ) match_mode=:any validate_units(bad_rx2)) == false
    @test (@test_logs (:warn, ) match_mode=:any validate_units(bad_rx3)) == false

    # Tests that Reactions created via the @reaction DSL with interpolation give warnings.
    bad_rx5 = @reaction $k1, $A --> $D
    bad_rx6 = @reaction $k1, $A --> $Bmol + $D
    bad_rx7 = @reaction $k1, $A + $D --> $Bmol
    @test (@test_logs (:warn, ) match_mode=:any validate_units(bad_rx5)) == false
    @test (@test_logs (:warn, ) match_mode=:any validate_units(bad_rx6)) == false
    @test (@test_logs (:warn, ) match_mode=:any validate_units(bad_rx7)) == false

    # Tests that creating systems with non-matching rate units yields warnings.
    @parameters k_bad [unit=us"kg"]
    bad_rxs = [
        Reaction(k_bad, nothing, [A]),
        Reaction(k_bad, [A], [B]),
        Reaction(k3, [A, B], [B], [1, 1], [2])
    ]
    @test_throws UnitValidationError ReactionSystem(bad_rxs, t, [A, B, C], [k1, k_bad, k3]; name = :bad_rs, unit_checks = true)
    @named bad_rs = ReactionSystem(bad_rxs, t, [A, B, C], [k1, k_bad, k3])
    @test (@test_logs (:warn, ) match_mode=:any validate_units(bad_rs)) == false
    @test_throws UnitValidationError assert_valid_units(bad_rs)
end


### DSL Tests (Molar Units) ###

# Tests that units work for DSL-based model creation with M units.
begin
    rs = @reaction_network begin
        @ivs t [unit=us"s"]
        @species begin
            A(t), [unit=us"M"]
            B(t), [unit=us"M"]
            C(t), [unit=us"M"]
        end
        @parameters begin
            k1, [unit=us"M/s"]
            k2, [unit=us"s^(-1)"]
            k3, [unit=us"M^(-1)*s^(-1)"]
        end
        k1, 0 --> A
        k2, A --> B
        k3, A + B --> 2B
    end

    # Checks that the ReactionSystem's content have the correct units.
    @test catalyst_get_unit(get_iv(rs)) == us"s"
    @test all(catalyst_get_unit.([rs.A, rs.B, rs.C]) .== us"M")
    @test catalyst_get_unit(rs.k1) == us"M/s"
    @test catalyst_get_unit(rs.k2) == us"s^(-1)"
    @test catalyst_get_unit(rs.k3) == us"M^(-1)*s^(-1)"
    for rx in reactions(rs)
        @test catalyst_get_unit(oderatelaw(rx)) == us"M/s"
    end

    # Checks that system declarations with erroneous units throw when unit_checks are enabled.
    @test_throws UnitValidationError @reaction_network begin
        @unit_checks true
        @ivs t [unit=us"1/s"] # Here, t's unit is wrong.
        @species begin
            A(t), [unit=us"M"]
            B(t), [unit=us"M"]
            C(t), [unit=us"M"]
        end
        @parameters begin
            k1, [unit=us"M/s"]
            k2, [unit=us"s^(-1)"]
            k3, [unit=us"M^(-1)*s^(-1)"]
        end
        k1, 0 --> A
        k2, A --> B
        k3, A + B --> 2B
    end
    @test_throws UnitValidationError @reaction_network begin
        @unit_checks true
        @ivs t [unit=us"s"]
        @species begin
            A(t), [unit=us"M"]
            B(t), [unit=us"M"]
            C(t), [unit=us"M"]
        end
        @parameters begin
            k1, [unit=us"M"] # Here, k1's unit is missing "/s".
            k2, [unit=us"s^(-1)"]
            k3, [unit=us"M^(-1)*s^(-1)"]
        end
        k1, 0 --> A
        k2, A --> B
        k3, A + B --> 2B
    end
    @test_throws UnitValidationError @reaction_network begin
        @unit_checks true
        @ivs t [unit=us"s"]
        @species begin
            A(t), [unit=us"M/s"] # Here, A's unit got an extra "/s".
            B(t), [unit=us"M"]
            C(t), [unit=us"M"]
        end
        @parameters begin
            k1, [unit=us"M/s"]
            k2, [unit=us"s^(-1)"]
            k3, [unit=us"M^(-1)*s^(-1)"]
        end
        k1, 0 --> A
        k2, A --> B
        k3, A + B --> 2B
    end

end

# Tests for the @unit_checks DSL option.
let
    # Checks that default (no @unit_checks) does NOT throw on bad units.
    @test_nowarn @reaction_network begin
        @ivs t [unit=us"s"]
        @species A(t) [unit=us"M"] B(t) [unit=us"M"]
        @parameters k_bad [unit=us"kg"]
        k_bad, A --> B
    end

    # Checks that @unit_checks true throws on bad units at construction time.
    @test_throws UnitValidationError @reaction_network begin
        @unit_checks true
        @ivs t [unit=us"s"]
        @species A(t) [unit=us"M"] B(t) [unit=us"M"]
        @parameters k_bad [unit=us"kg"]
        k_bad, A --> B
    end

    # Checks that @unit_checks false allows constructing a system with incompatible units,
    # and that the resulting system matches an equivalent programmatic construction.
    rs_dsl = @reaction_network begin
        @unit_checks false
        @ivs t [unit=us"s"]
        @species begin
            A(t), [unit=us"M"]
            B(t), [unit=us"M"]
        end
        @parameters begin
            k_bad, [unit=us"kg"]
        end
        k_bad, A --> B
    end

    @independent_variables t [unit=us"s"]
    @species A(t) [unit=us"M"] B(t) [unit=us"M"]
    @parameters k_bad [unit=us"kg"]
    rs_prog = complete(ReactionSystem(
        [Reaction(k_bad, [A], [B])], t; name=:rs
    ))
    @test Catalyst.isequivalent(rs_dsl, rs_prog)

    # Post-hoc validation still detects the bad units.
    @test validate_units(rs_dsl; warn=false) == false
    @test_throws UnitValidationError assert_valid_units(rs_dsl)

    # Checks that @unit_checks false with bad units still allows ODE problem creation.
    @unpack A, B, k_bad = rs_dsl
    u0 = [A => 1.0, B => 0.0]
    tspan = (0.0, 1.0)
    ps = [k_bad => 0.5]
    oprob = ODEProblem(rs_dsl, u0, tspan, ps)
    @test oprob isa ODEProblem

    # Checks that @unit_checks false works with @network_component (incomplete systems).
    rs_comp = @network_component begin
        @unit_checks false
        @ivs t [unit=us"s"]
        @species X(t) [unit=us"mol"]
        @parameters d [unit=us"m"]
        d, X --> 0
    end
    @test !Catalyst.iscomplete(rs_comp)
    @test length(species(rs_comp)) == 1
    @test length(reactions(rs_comp)) == 1

    # Checks @unit_checks false with a multi-reaction system where some reactions have
    # correct units and others don't, and verifies equivalence with programmatic construction.
    rs_mixed_dsl = @reaction_network begin
        @unit_checks false
        @ivs t [unit=us"s"]
        @species begin
            A(t), [unit=us"M"]
            B(t), [unit=us"M"]
            C(t), [unit=us"M"]
        end
        @parameters begin
            k_good, [unit=us"s^(-1)"]
            k_bad, [unit=us"kg"]
        end
        k_good, A --> B
        k_bad, B --> C
    end

    @independent_variables t [unit=us"s"]
    @species A(t) [unit=us"M"] B(t) [unit=us"M"] C(t) [unit=us"M"]
    @parameters k_good [unit=us"s^(-1)"] k_bad [unit=us"kg"]
    rs_mixed_prog = complete(ReactionSystem(
        [Reaction(k_good, [A], [B]), Reaction(k_bad, [B], [C])], t; name=:rs
    ))
    @test Catalyst.isequivalent(rs_mixed_dsl, rs_mixed_prog)
    @test validate_units(rs_mixed_dsl; warn=false) == false
end


### μM Unit Tests ###

# Tests with μM units (programmatic).
let
    @independent_variables t [unit=us"s"]
    @species A(t) [unit=us"μM"] B(t) [unit=us"μM"]
    @parameters k1 [unit=us"μM/s"] k2 [unit=us"s^(-1)"] k3 [unit=us"μM^(-1)*s^(-1)"]
    rxs = [Reaction(k1, nothing, [A]),
        Reaction(k2, [A], [B]),
        Reaction(k3, [A, B], [B], [1, 1], [2])]
    @named rs = ReactionSystem(rxs, t, [A, B], [k1, k2, k3])
    rs = complete(rs)
    @test validate_units(rs)
    @test_nowarn ode_model(rs)
    @test_nowarn sde_model(rs)
    @test_nowarn jump_model(rs)
    @test_nowarn ss_ode_model(rs)
end


### Compound Rate Tests ###

# Tests compound rates with M units (the key bug fix: non-SI compound rates
# previously failed due to FP precision loss in MTKBase's unit expansion).
let
    @independent_variables t [unit=us"s"]
    @species A(t) [unit=us"M"] B(t) [unit=us"M"]
    @parameters k [unit=us"M^(-2)*s^(-1)"]

    # k*A, 2A --> B: full rate = k*A * A^2 = k*A^3. Units: M^(-2)s^(-1) * M^3 = M/s ✓
    rx = Reaction(k * A, [A], [B], [2], [1])
    @named rs = ReactionSystem([rx], t, [A, B], [k])
    rs = complete(rs)
    @test validate_units(rs)
end

# Tests compound rates with μM units.
let
    @independent_variables t [unit=us"s"]
    @species A(t) [unit=us"μM"] B(t) [unit=us"μM"]
    @parameters k [unit=us"μM^(-2)*s^(-1)"]

    rx = Reaction(k * A, [A], [B], [2], [1])
    @named rs = ReactionSystem([rx], t, [A, B], [k])
    rs = complete(rs)
    @test validate_units(rs)
end


### only_use_rate Tests ###

# Tests that only_use_rate=true works correctly with units (programmatic).
let
    @independent_variables t [unit=us"s"]
    @species X(t) [unit=us"M"] Y(t) [unit=us"M"]
    @parameters v [unit=us"M/s"] K [unit=us"M"]

    # Custom rate law with only_use_rate=true: rate = mm(X, v, K).
    # Units of mm(...) = M/s directly (no substrate multiplication).
    rx = Reaction(mm(X, v, K), [X], [Y]; only_use_rate = true)
    @named rs = ReactionSystem([rx], t, [X, Y], [v, K])
    rs = complete(rs)
    @test validate_units(rs)
end

# Tests that only_use_rate=true works via DSL (=> arrow).
let
    rs = @reaction_network begin
        @ivs t [unit=us"s"]
        @species begin
            X(t), [unit=us"M"]
            Y(t), [unit=us"M"]
        end
        @parameters begin
            v, [unit=us"M/s"]
            K, [unit=us"M"]
        end
        mm(X, v, K), X => Y
    end
    @test validate_units(rs)
end


### Non-Trivial Rate Tests ###

# Tests non-trivial rates (no units).
let
    # Parametric rates.
    @test_nowarn @reaction_network begin
        k1, 2X1 --> Z1
        k2, n2*X2 --> m2*Z2
        k3, n3*X3 + m3*Y3 --> Z3
    end

    # Non-trivial rates with registered functions.
    @test_nowarn @reaction_network begin
        k1*X1, 2X1 --> Z1
        mm(X2, v2, K2), X2 --> Z2
        hill(X3, v3, K3, n3), X3 + Y3 --> Z3
    end
end

# Tests non-trivial rates with M units.
let
    rs = @reaction_network begin
        @ivs t [unit=us"s"]
        @species begin
            X1(t), [unit=us"M"]
            Z1(t), [unit=us"M"]
            X2(t), [unit=us"M"]
            Z2(t), [unit=us"M"]
            X3(t), [unit=us"M"]
            Y3(t), [unit=us"M"]
            Z3(t), [unit=us"M"]
        end
        @parameters begin
            k1, [unit=us"M^(-2)*s^(-1)"]
            v2, [unit=us"M^(-2)*s^(-1)"]
            K2, [unit=us"M"]
            v3, [unit=us"M^(-1)*s^(-1)"]
            K3, [unit=us"M"]
            n3
        end
        k1*X1, 2X1 --> Z1
        mm(X2, v2, K2), 3X2 --> Z2
        hill(X3, v3, K3, n3), X3 + Y3 --> Z3
    end
    @test validate_units(rs)
end


### Stochastic Chemical Kinetics Units ###

# Tests unitless species with time in seconds (stochastic CK convention):
# species as molecule counts (no unit metadata → dimensionless),
# time in seconds, all rate constants in s^(-1).
let
    @independent_variables t [unit=us"s"]
    @species S(t) I(t) R(t) # No unit metadata → unitless (molecule counts).
    @parameters k0 [unit=us"s^(-1)"] k1 [unit=us"s^(-1)"] k2 [unit=us"s^(-1)"]

    # All rates are s^(-1) since species are dimensionless.
    rxs = [
        Reaction(k0, nothing, [S]),             # 0th order
        Reaction(k1, [S], [I]),                 # 1st order
        Reaction(k2, [S, I], [R], [1, 1], [1])  # 2nd order
    ]
    @named rs = ReactionSystem(rxs, t, [S, I, R], [k0, k1, k2])
    rs = complete(rs)
    @test validate_units(rs)
    @test_nowarn ode_model(rs)
    @test_nowarn jump_model(rs)
end


### Model Conversion Tests ###

# Tests that model conversions work without warnings with M units.
let
    @independent_variables t [unit=us"s"]
    @species A(t) [unit=us"M"] B(t) [unit=us"M"]
    @parameters k1 [unit=us"M/s"] k2 [unit=us"s^(-1)"] k3 [unit=us"M^(-1)*s^(-1)"]
    rxs = [
        Reaction(k1, nothing, [A]),
        Reaction(k2, [A], [B]),
        Reaction(k3, [A, B], [B], [1, 1], [2])
    ]
    @named rs = ReactionSystem(rxs, t, [A, B], [k1, k2, k3])
    rs = complete(rs)
    @test_nowarn ode_model(rs)
    @test_nowarn sde_model(rs)
    @test_nowarn jump_model(rs)
    @test_nowarn ss_ode_model(rs)
end

# Tests that model conversions work without warnings with μM units.
let
    @independent_variables t [unit=us"s"]
    @species A(t) [unit=us"μM"] B(t) [unit=us"μM"]
    @parameters k1 [unit=us"μM/s"] k2 [unit=us"s^(-1)"]
    rxs = [
        Reaction(k1, nothing, [A]),
        Reaction(k2, [A], [B])
    ]
    @named rs = ReactionSystem(rxs, t, [A, B], [k1, k2])
    rs = complete(rs)
    @test_nowarn ode_model(rs)
    @test_nowarn sde_model(rs)
    @test_nowarn jump_model(rs)
    @test_nowarn ss_ode_model(rs)
end


### Functional Parameters with Units ###

# Tests the combination of functional parameters and units.
let
    ts = collect(0.0:0.01:1.0)
    spline = LinearInterpolation(1 ./ (1 .+ ts), ts)
    @parameters (pIn::typeof(spline))(..) [unit=us"1/(s*m^3)"]

    @independent_variables t [unit=us"s"]
    @species X(t) [unit=us"mol/m^3"]
    @parameters p_base [unit=us"mol"] d [unit=us"1/s"]
    rxs = [
        Reaction(p_base*pIn(t), [], [X])
        Reaction(d, [X], [])
    ]
    @named rs = ReactionSystem(rxs, t)

    # catalyst_get_unit works on both the bare callable parameter and the call form.
    @test catalyst_get_unit(pIn) == us"1/(s*m^3)"
    @test catalyst_get_unit(pIn(t)) == us"1/(s*m^3)"

    # Reaction rates have the correct units.
    rxs_out = reactions(rs)
    @test catalyst_get_unit(rxs_out[1].rate) == us"mol/(s*m^3)"
    @test catalyst_get_unit(rxs_out[2].rate) == us"1/s"
end


### Brownians with Units ###

# Tests brownian noise with M units. Brownians have effective units of
# time^(-1/2) (Wiener process derivative dW/dt), so the noise coefficient σ
# needs units M*s^(-1/2) for the equation D(V) ~ -k*V + σ*B to balance.
let
    @brownians B
    @independent_variables t [unit=us"s"]
    @variables V(t) [unit=us"M"]
    @species S(t) [unit=us"M"] P(t) [unit=us"M"]
    σ_unit = us"M" / us"s"^(1//2)
    @parameters k [unit=us"s^(-1)"] σ [unit=σ_unit]
    D = Differential(t)

    @named rn = ReactionSystem(
        [Reaction(k, [S], [P]), D(V) ~ -k * V + σ * B],
        t, [S, P, V], [k, σ], [B]
    )
    rn = complete(rn)
    @test validate_units(rn)
end


### Poissonians with Units ###

# Tests poissonian with M units. Poissonians inherit units from their rate
# parameter (λ has units M/s), so D(X) ~ dN balances as M/s = M/s.
let
    @independent_variables t [unit=us"s"]
    @parameters λ [unit=us"M/s"] k [unit=us"s^(-1)"]
    @species S(t) [unit=us"M"]
    @variables X(t) [unit=us"M"]
    @poissonians dN(λ)
    D = Differential(t)

    @named rn = ReactionSystem(
        [Reaction(k, [S], nothing), D(X) ~ dN],
        t, [S, X], [λ, k];
        poissonians = [dN]
    )
    rn = complete(rn)
    @test validate_units(rn)
end

# Tests poissonian with unitless species (stochastic CK). Rate λ has units
# s^(-1), so dN has units s^(-1), matching D(X) = unitless/s = s^(-1).
let
    @independent_variables t [unit=us"s"]
    @parameters λ [unit=us"s^(-1)"] k [unit=us"s^(-1)"]
    @species S(t)
    @variables X(t)
    @poissonians dN(λ)
    D = Differential(t)

    @named rn = ReactionSystem(
        [Reaction(k, [S], nothing), D(X) ~ dN],
        t, [S, X], [λ, k];
        poissonians = [dN]
    )
    rn = complete(rn)
    @test validate_units(rn)
end


### Array Variables with Units ###

# Tests unit inference for indexed array species and parameters.
let
    @independent_variables t [unit=us"s"]
    @species (X(t))[1:3] [unit=us"M"]
    @parameters (k)[1:3] [unit=us"s^(-1)"]

    # Individual array elements retain the parent's unit.
    @test catalyst_get_unit(X[1]) == us"M"
    @test catalyst_get_unit(X[3]) == us"M"
    @test catalyst_get_unit(k[1]) == us"s^(-1)"

    # Expressions with array elements.
    @test catalyst_get_unit(k[1] * X[2]) == us"M/s"
end

# Tests validation of a system with array species.
let
    @independent_variables t [unit=us"s"]
    @species (X(t))[1:2] [unit=us"M"]
    @parameters k [unit=us"s^(-1)"]
    rxs = [
        Reaction(k, [X[1]], [X[2]])
    ]
    @named rs = ReactionSystem(rxs, t)
    rs = complete(rs)
    @test validate_units(rs)
end


### Equation / Differential Unit Tests ###

# Tests catalyst_get_unit on differential expressions with various units.
let
    @independent_variables t [unit=us"s"]
    @species X(t) [unit=us"M"]
    @variables V(t) [unit=us"μM"]
    @parameters k [unit=us"s^(-1)"] v [unit=us"M/s"] K [unit=us"M"]
    D = Differential(t)

    # Basic first-order derivatives.
    @test catalyst_get_unit(D(X)) == us"M/s"
    @test catalyst_get_unit(D(V)) == us"μM/s"

    # Derivative of a compound expression: D(k*X) = s^(-1)*M/s = M/s^2.
    @test catalyst_get_unit(D(k * X)) == us"M/s^2"

    # Higher-order derivatives.
    @test catalyst_get_unit(D(D(X))) == us"M/s^2"
    @test catalyst_get_unit(D(D(V))) == us"μM/s^2"
end

# Tests basic equation validation with a single ODE equation (no reactions).
let
    @independent_variables t [unit=us"s"]
    @variables V(t) [unit=us"M"]
    @parameters k [unit=us"s^(-1)"]
    D = Differential(t)

    # D(V) ~ -k*V: M/s = s^(-1)*M = M/s ✓
    @named rs = ReactionSystem([D(V) ~ -k * V], t, [V], [k])
    rs = complete(rs)
    @test validate_units(rs)
end

# Tests equation validation with compound RHS expression (registered functions).
let
    @independent_variables t [unit=us"s"]
    @species X(t) [unit=us"M"]
    @variables V(t) [unit=us"M"]
    @parameters v [unit=us"M/s"] K [unit=us"M"] k [unit=us"s^(-1)"]
    D = Differential(t)

    # D(V) ~ mm(X, v, K) - k*V: M/s = M/s - s^(-1)*M = M/s ✓
    @named rs = ReactionSystem(
        [Reaction(k, [X], nothing), D(V) ~ mm(X, v, K) - k * V],
        t, [X, V], [v, K, k]
    )
    rs = complete(rs)
    @test validate_units(rs)
end

# Tests equation validation with μM units.
let
    @independent_variables t [unit=us"s"]
    @variables V(t) [unit=us"μM"]
    @parameters k [unit=us"s^(-1)"] p [unit=us"μM/s"]
    D = Differential(t)

    # D(V) ~ p - k*V: μM/s = μM/s - s^(-1)*μM = μM/s ✓
    @named rs = ReactionSystem([D(V) ~ p - k * V], t, [V], [k, p])
    rs = complete(rs)
    @test validate_units(rs)
end

# Tests system with multiple equations.
let
    @independent_variables t [unit=us"s"]
    @variables V(t) [unit=us"M"] W(t) [unit=us"M"]
    @parameters k1 [unit=us"s^(-1)"] k2 [unit=us"s^(-1)"] v [unit=us"M/s"]
    D = Differential(t)

    # Two coupled ODEs.
    @named rs = ReactionSystem(
        [D(V) ~ v - k1 * V, D(W) ~ k1 * V - k2 * W],
        t, [V, W], [k1, k2, v]
    )
    rs = complete(rs)
    @test validate_units(rs)
end

# Tests mixed reactions + equations system.
let
    @independent_variables t [unit=us"s"]
    @species S(t) [unit=us"M"] P(t) [unit=us"M"]
    @variables V(t) [unit=us"M"]
    @parameters k [unit=us"s^(-1)"] d [unit=us"s^(-1)"]
    D = Differential(t)

    @named rs = ReactionSystem(
        [Reaction(k, [S], [P]), D(V) ~ k * S - d * V],
        t, [S, P, V], [k, d]
    )
    rs = complete(rs)
    @test validate_units(rs)
end

# Tests equation with power and division expressions in RHS.
let
    @independent_variables t [unit=us"s"]
    @variables V(t) [unit=us"M"]
    @parameters k [unit=us"M^(-1)*s^(-1)"] K [unit=us"M"]
    D = Differential(t)

    # D(V) ~ k*V^2: M/s = M^(-1)s^(-1) * M^2 = M/s ✓
    @named rs = ReactionSystem([D(V) ~ k * V^2], t, [V], [k, K])
    rs = complete(rs)
    @test validate_units(rs)

    # Rational expression: D(V) ~ k*V^2*K/(K+V)
    # Units: M^(-1)s^(-1) * M^2 * M / M = M/s ✓
    @named rs2 = ReactionSystem([D(V) ~ k * V^2 * K / (K + V)], t, [V], [k, K])
    rs2 = complete(rs2)
    @test validate_units(rs2)
end

# Tests equation with unitless species (stochastic CK style).
let
    @independent_variables t [unit=us"s"]
    @variables X(t) Y(t) # No units → dimensionless.
    @parameters k1 [unit=us"s^(-1)"] k2 [unit=us"s^(-1)"]
    D = Differential(t)

    @named rs = ReactionSystem(
        [D(X) ~ -k1 * X, D(Y) ~ k1 * X - k2 * Y],
        t, [X, Y], [k1, k2]
    )
    rs = complete(rs)
    @test validate_units(rs)
end


### Bad Equation Unit Detection ###

# Tests that equations with mismatched LHS/RHS units are caught.
let
    @independent_variables t [unit=us"s"]
    @variables V(t) [unit=us"M"]
    @species S(t) [unit=us"M"]
    @parameters k [unit=us"s^(-1)"] p [unit=us"kg"]
    D = Differential(t)

    # D(V) has units M/s, but p has units kg → mismatch.
    @test_throws UnitValidationError ReactionSystem(
        [Reaction(k, [S], nothing), D(V) ~ p],
        t, [S, V], [k, p];
        name = :rn_bad_lhs_rhs, unit_checks = true
    )
    @named rn = ReactionSystem(
        [Reaction(k, [S], nothing), D(V) ~ p],
        t, [S, V], [k, p]
    )
    rn = complete(rn)
    @test (@test_logs (:warn,) match_mode=:any validate_units(rn)) == false
    @test_throws UnitValidationError assert_valid_units(rn)
end

# Tests that additive terms with mismatched units within an equation are caught.
let
    @independent_variables t [unit=us"s"]
    @variables V(t) [unit=us"M"]
    @species S(t) [unit=us"M"]
    @parameters k [unit=us"s^(-1)"] p [unit=us"kg/s"]
    D = Differential(t)

    # D(V) ~ -k*V + p: first term is M/s, second is kg/s → internal mismatch.
    @test_throws UnitValidationError ReactionSystem(
        [Reaction(k, [S], nothing), D(V) ~ -k * V + p],
        t, [S, V], [k, p];
        name = :rn_bad_add_terms, unit_checks = true
    )
    @named rn = ReactionSystem(
        [Reaction(k, [S], nothing), D(V) ~ -k * V + p],
        t, [S, V], [k, p]
    )
    rn = complete(rn)
    @test (@test_logs (:warn,) match_mode=:any validate_units(rn)) == false
    @test_throws UnitValidationError assert_valid_units(rn)
end

# Tests that equation-only systems are still unit checked.
let
    @independent_variables t [unit=us"s"]
    @variables V(t) [unit=us"M"]
    @parameters p [unit=us"kg"]
    D = Differential(t)

    @test_throws UnitValidationError ReactionSystem([D(V) ~ p], t, [V], [p]; name = :eq_only_bad, unit_checks = true)
    @named eq_only_bad = ReactionSystem([D(V) ~ p], t, [V], [p])
    @test (@test_logs (:warn,) match_mode=:any validate_units(eq_only_bad)) == false
end


### Conditional Expressions (ifelse) ###

# Tests unit inference through ifelse expressions.
let
    @independent_variables t [unit=us"s"]
    @species X(t) [unit=us"M"]
    @parameters v [unit=us"M/s"] K [unit=us"M"]

    # ifelse should return the unit of the true/false branches.
    cond_expr = ifelse(X > K, v, v / 2)
    @test catalyst_get_unit(cond_expr) == us"M/s"
end


### Comparison Expression Units ###

# Tests comparison expressions infer as unitless.
let
    @independent_variables t [unit=us"s"]
    @species X(t) [unit=us"M"]
    @parameters K [unit=us"M"] τ [unit=us"s"]

    @test catalyst_get_unit(X > K) == SYM_UNITLESS
    @test catalyst_get_unit(X < K) == SYM_UNITLESS
    @test catalyst_get_unit(X == K) == SYM_UNITLESS
    @test catalyst_get_unit(X > τ) == SYM_UNITLESS
end


### ifelse Validation ###

# Tests that ifelse condition and branch rules are unit-checked.
let
    @independent_variables t [unit=us"s"]
    @species X(t) [unit=us"M"]
    @variables V(t) [unit=us"M"]
    @parameters v [unit=us"M/s"] K [unit=us"M"] τ [unit=us"s"]
    D = Differential(t)

    # Comparison in condition must use compatible units.
    @test_throws UnitValidationError ReactionSystem(
        [D(V) ~ ifelse(X > τ, v, v / 2)],
        t, [X, V], [v, K, τ];
        name = :ifelse_bad_comparison, unit_checks = true
    )
    @named rs_bad_comparison = ReactionSystem(
        [D(V) ~ ifelse(X > τ, v, v / 2)],
        t, [X, V], [v, K, τ]
    )
    report_bad_comparison = Catalyst.unit_validation_report(rs_bad_comparison)
    @test !report_bad_comparison.valid
    @test any(issue -> issue.kind == :comparison_unit_mismatch, report_bad_comparison.issues)

    # Branches must have matching units.
    @test_throws UnitValidationError ReactionSystem(
        [D(V) ~ ifelse(X > K, v, X)],
        t, [X, V], [v, K, τ];
        name = :ifelse_bad_branches, unit_checks = true
    )
    @named rs_bad_branches = ReactionSystem(
        [D(V) ~ ifelse(X > K, v, X)],
        t, [X, V], [v, K, τ]
    )
    report_bad_branches = Catalyst.unit_validation_report(rs_bad_branches)
    @test !report_bad_branches.valid
    @test any(issue -> issue.kind == :conditional_branch_unit_mismatch, report_bad_branches.issues)
end


### Exponent Unit Validation ###

# Tests that exponents must be unitless and do not crash validation.
let
    @independent_variables t [unit=us"s"]
    @variables V(t) [unit=us"M"]
    D = Differential(t)

    @test_throws UnitValidationError ReactionSystem([D(V) ~ V^t], t, [V], []; name = :bad_exponent_units, unit_checks = true)
    @named rs_bad_exponent = ReactionSystem([D(V) ~ V^t], t, [V], [])
    report_bad_exponent = Catalyst.unit_validation_report(rs_bad_exponent)
    @test !report_bad_exponent.valid
    @test any(issue -> issue.kind == :exponent_unit_mismatch, report_bad_exponent.issues)
end


### Higher-Order Derivatives ###

# Tests unit inference for second-order derivatives.
let
    @independent_variables t [unit=us"s"]
    @species X(t) [unit=us"M"]
    D = Differential(t)

    # D(D(X)) should have units M/s^2.
    @test catalyst_get_unit(D(D(X))) == us"M/s^2"
end


### Symbolic Stoichiometry ###

# Tests that symbolic stoichiometry produces a validation issue rather than crashing.
let
    @independent_variables t [unit=us"s"]
    @species A(t) [unit=us"M"] B(t) [unit=us"M"]
    @parameters k [unit=us"s^(-1)"] n

    # Symbolic stoichiometry should produce a :symbolic_stoichiometry issue.
    @test_throws UnitValidationError ReactionSystem(
        [Reaction(k, [A], [B], [n], [1])], t, [A, B], [k, n];
        name = :sym_stoich, unit_checks = true
    )
    @named rs_sym = ReactionSystem([Reaction(k, [A], [B], [n], [1])], t, [A, B], [k, n])
    report = Catalyst.unit_validation_report(rs_sym)
    @test !report.valid
    @test any(issue -> issue.kind == :symbolic_stoichiometry, report.issues)

    # Non-symbolic (numeric) stoichiometry still works fine.
    @test_nowarn ReactionSystem(
        [Reaction(k, [A], [B])], t, [A, B], [k, n];
        name = :num_stoich, unit_checks = true
    )
end


### Symbolic Exponents ###

# Tests that symbolic exponents on unitful bases produce :symbolic_exponent issues
# rather than returning broken (non-Quantity) types from catalyst_get_unit.
let
    @independent_variables t [unit=us"s"]
    @species A(t) [unit=us"M"] B(t) [unit=us"M"]
    @parameters v [unit=us"M/s"] K [unit=us"M"] n k [unit=us"s^(-1)"]

    # Raw symbolic exponent in rate (full rate law, only_use_rate=true):
    # v * A^n / (K^n + A^n) should produce a :symbolic_exponent issue.
    rx_hill = Reaction(v * A^n / (K^n + A^n), [A], [B]; only_use_rate = true)
    @named rs_hill = ReactionSystem([rx_hill], t, [A, B], [v, K, n])
    report = Catalyst.unit_validation_report(rs_hill)
    @test !report.valid
    @test any(issue -> issue.kind == :symbolic_exponent, report.issues)

    # With unit_checks=true, the ReactionSystem constructor should throw.
    @test_throws UnitValidationError ReactionSystem(
        [rx_hill], t, [A, B], [v, K, n];
        name = :sym_exp, unit_checks = true
    )

    # Reaction-level validation also catches symbolic exponents in rate.
    rx_report = Catalyst.unit_validation_report(rx_hill)
    @test !rx_report.valid
    @test any(issue -> issue.kind == :symbolic_exponent, rx_report.issues)

    # Numeric exponents still work fine (full rate law, only_use_rate=true).
    rx_num = Reaction(v * A^2 / (K^2 + A^2), [A], [B]; only_use_rate = true)
    @named rs_num = ReactionSystem([rx_num], t, [A, B], [v, K, n])
    @test Catalyst.unit_validation_report(rs_num).valid

    # Registered hill function bypasses ispow and works correctly.
    # Note: uses only_use_rate=false (default) — the rate constant hill(A,v,K,n) has
    # units M^(-1)*s^(-1) here, so rate_const * A has correct units M/s.
    @parameters v2 [unit=us"s^(-1)"]
    rx_reg = Reaction(hill(A, v2, K, n), [A], [B])
    @named rs_reg = ReactionSystem([rx_reg], t, [A, B], [v2, K, n])
    @test Catalyst.unit_validation_report(rs_reg).valid

    # Symbolic exponent on unitless base is fine (no units to break).
    rx_unitless_base = Reaction(k * 2^n, [A], [B])
    @named rs_unitless = ReactionSystem([rx_unitless_base], t, [A, B], [k, n])
    @test Catalyst.unit_validation_report(rs_unitless).valid
end


### Reaction-Level unit_checks Tests ###

# Tests that Reaction constructor with unit_checks=true validates correctly.
let
    @independent_variables t [unit=us"s"]
    @species A(t) [unit=us"M"] B(t) [unit=us"M"] C(t) [unit=us"M"]
    @parameters k1 [unit=us"M/s"] k2 [unit=us"s^(-1)"] k3 [unit=us"M^(-1)*s^(-1)"]

    # Valid reactions pass with unit_checks=true.
    @test_nowarn Reaction(k1, nothing, [A]; unit_checks = true)
    @test_nowarn Reaction(k2, [A], [B]; unit_checks = true)
    @test_nowarn Reaction(k3, [A, B], [C]; unit_checks = true)

    # Valid reaction with non-trivial stoichiometries.
    @test_nowarn Reaction(k3, [A, B], [B], [1, 1], [2]; unit_checks = true)

    # Bad-unit Reaction does NOT throw with default (unit_checks=false).
    @parameters k_bad [unit=us"kg"]
    @test_nowarn Reaction(k_bad, [A], [B])

    # Bad-unit Reaction throws with unit_checks=true (species unit mismatch).
    @species D_kg(t) [unit=us"kg"]
    @test_throws UnitValidationError Reaction(k2, [A], [D_kg]; unit_checks = true)

    # Rate with internally inconsistent additive terms caught at Reaction construction.
    @parameters a [unit=us"M/s"] z [unit=us"kg/s"]
    @test_throws UnitValidationError Reaction(a + z, nothing, [A]; unit_checks = true)

    # Same inconsistent rate does NOT throw without unit_checks.
    @test_nowarn Reaction(a + z, nothing, [A])
end

# Tests Reaction unit_checks via DSL with @unit_checks true.
let
    # Valid DSL reaction network with unit_checks.
    @test_nowarn @reaction_network begin
        @unit_checks true
        @ivs t [unit=us"s"]
        @species A(t) [unit=us"M"] B(t) [unit=us"M"]
        @parameters k [unit=us"s^(-1)"]
        k, A --> B
    end

    # DSL with @unit_checks true catches bad rate units at Reaction construction.
    @test_throws UnitValidationError @reaction_network begin
        @unit_checks true
        @ivs t [unit=us"s"]
        @species A(t) [unit=us"M"] B(t) [unit=us"M"]
        @parameters k_bad [unit=us"kg"]
        k_bad, A --> B
    end
end
