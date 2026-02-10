#! format: off

### Prepares Tests ###

# Fetch packages.
using Catalyst, DynamicQuantities, Test
using ModelingToolkitBase: get_iv, get_unit, validate, ValidationError


### Basic Tests ###

# Checks that units work with programmatic model creation.
let
    # Creates a `ReactionSystem` programmatically, while designating units.
    @independent_variables t [unit=u"s"]
    @species A(t) [unit=u"mol/m^3"] B(t) [unit=u"mol/m^3"] C(t) [unit=u"mol/m^3"]
    @parameters k1 [unit=u"mol/(m^3*s)"] k2 [unit=u"s^(-1)"] k3 [unit=u"(m^3)/(s*mol)"]
    rxs = [Reaction(k1, nothing, [A]),
        Reaction(k2, [A], [B]),
        Reaction(k3, [A, B], [B], [1, 1], [2])]
    @test_nowarn ReactionSystem(rxs, t, [A, B, C], [k1, k2, k3]; name = :rs)
    @named rs = ReactionSystem(rxs, t, [A, B, C], [k1, k2, k3])
    rs = complete(rs)

    # Test that all reactions have the correct unit.
    for rx in reactions(rs)
        @test get_unit(oderatelaw(rx)) == u"mol/(s*m^3)"
        # we don't currently convert units, so they will be the same as for ODEs
        @test get_unit(jumpratelaw(rx)) == u"mol/(s*m^3)"
    end

    # Tests that the system can be converted to MTK systems without warnings.
    @test_nowarn ode_model(rs)
    @test_nowarn sde_model(rs)
    @test_nowarn jump_model(rs)
    @test_nowarn ss_ode_model(rs)

    # Tests that creating `Reaction`s with non-matching units yields warnings.
    @species B(t) [unit=u"mol"] D(t) [unit=u"kg"]
    bad_rx1 = Reaction(k1, [A], [D])
    bad_rx2 = Reaction(k1, [A], [B, D])
    bad_rx3 = Reaction(k1, [A, D], [B])
    bad_rx4 = Reaction(k1 + k2, [A], [B])
    @test (@test_logs (:warn, ) match_mode=:any validate(bad_rx1)) == false
    @test (@test_logs (:warn, ) match_mode=:any validate(bad_rx2)) == false
    @test (@test_logs (:warn, ) match_mode=:any validate(bad_rx3)) == false
    @test (@test_logs (:warn, ) match_mode=:any validate(bad_rx4)) == false

    # Tests that `Reactions` created via the `@reaction` DSL and interpolation gives warnings.
    bad_rx5 = @reaction $k1, $A --> $D
    bad_rx6 = @reaction $k1, $A --> $B + $D
    bad_rx7 = @reaction $k1, $A + $D --> $B
    bad_rx8 = @reaction $(k1+k2), $A --> $B
    @test (@test_logs (:warn, ) match_mode=:any validate(bad_rx5)) == false
    @test (@test_logs (:warn, ) match_mode=:any validate(bad_rx6)) == false
    @test (@test_logs (:warn, ) match_mode=:any validate(bad_rx7)) == false
    @test (@test_logs (:warn, ) match_mode=:any validate(bad_rx8)) == false

    # Tests that creating systems with non-matching units yields warnings.
    @parameters k2 [unit=u"kg"]
    bad_rxs = [
        Reaction(k2, nothing, [A]),
        Reaction(k2, [A], [B]),
        Reaction(k3, [A, B], [B], [1, 1], [2])
    ]
    @named bad_rs = ReactionSystem(bad_rxs, t, [A, B, C], [k1, k2, k3])
    @test (@test_logs (:warn, ) match_mode=:any validate(bad_rs)) == false
end

# Checks that units work for DSL-based model creation.
begin
    # Creates the model, while designating matching units.
    rs = @reaction_network begin
        @ivs t [unit=u"s"]
        @species begin
            A(t), [unit=u"mol/m^3"]
            B(t), [unit=u"mol/m^3"]
            C(t), [unit=u"mol/m^3"]
        end
        @parameters begin
            k1, [unit=u"mol/(s*m^3)"]
            k2, [unit=u"s^(-1)"]
            k3, [unit=u"(m^3)/(s*mol)"]
        end
        k1, 0 --> A
        k2, A --> B
        k3, A + B --> 2B
    end

    # Checks that the `ReactionSystem`'s content have the correct units.
    @test get_unit(get_iv(rs)) == u"s"
    @test all(get_unit.([rs.A, rs.B, rs.C]) .== [u"mol/m^3", u"mol/m^3", u"mol/m^3"])
    @test all(get_unit.([rs.k1, rs.k2, rs.k3]) .== [u"mol/(s*m^3)", u"s^(-1)", u"(m^3)/(s*mol)"])
    for rx in reactions(rs)
        @test get_unit(oderatelaw(rx)) == u"mol/(s*m^3)"
        # we don't currently convert units, so they will be the same as for ODEs
        @test get_unit(jumpratelaw(rx)) == u"mol/(s*m^3)"
    end

    # Checks that system declarations with erroneous units yields errors.
    @test_logs (:warn, ) match_mode=:any @reaction_network begin
        @ivs t [unit=u"1/s"] # Here, t's unit is wrong.
        @species begin
            A(t), [unit=u"mol/m^3"]
            B(t), [unit=u"mol/m^3"]
            C(t), [unit=u"mol/m^3"]
        end
        @parameters begin
            k1, [unit=u"mol/(s*m^3)"]
            k2, [unit=u"s^(-1)"]
            k3, [unit=u"(m^3)/(s*mol)"]
        end
        k1, 0 --> A
        k2, A --> B
        k3, A + B --> 2B
    end
    @test_logs (:warn, ) match_mode=:any @reaction_network begin
        @ivs t [unit=u"s"]
        @species begin
            A(t), [unit=u"mol/m^3"]
            B(t), [unit=u"mol/m^3"]
            C(t), [unit=u"mol/m^3"]
        end
        @parameters begin
            k1, [unit=u"mol/(m^3)"] # Here, k1's unit is missing "/s".
            k2, [unit=u"s^(-1)"]
            k3, [unit=u"(m^3)/(s*mol)"]
        end
        k1, 0 --> A
        k2, A --> B
        k3, A + B --> 2B
    end
    @test_logs (:warn, ) match_mode=:any @reaction_network begin
        @ivs t [unit=u"s"]
        @species begin
            A(t), [unit=u"mol/(s*m^3)"]  # Here, A's unit got an extra "/s".
            B(t), [unit=u"mol/m^3"]
            C(t), [unit=u"mol/m^3"]
        end
        @parameters begin
            k1, [unit=u"mol/(s*m^3)"]
            k2, [unit=u"s^(-1)"]
            k3, [unit=u"(m^3)/(s*mol)"]
        end
        k1, 0 --> A
        k2, A --> B
        k3, A + B --> 2B
    end
end

# Checks that units works for various non-trivial systems.
let
    # Parametric rates (no units).
    @test_nowarn @reaction_network begin
        k1, 2X1 --> Z1
        k2, n2*X2 --> m2*Z2
        k3, n3*X3 + m3*Y3--> Z3
    end

    # Parametric rates with units is not possible, since the unit depends on the parameter values.

    # Non-trivial rates (no units).
    @test_nowarn @reaction_network begin
        k1*X1, 2X1 --> Z1
        mm(X2, v2, K2), X2 --> Z2
        hill(X3, v3, K3, n3), X3 + Y3--> Z3
    end

    # Non-trivial rates (units).
    @test_nowarn @reaction_network begin
        @ivs t [unit=u"s"]
        @species begin
            X1(t), [unit=u"mol/m^3"]
            Z1(t), [unit=u"mol/m^3"]
            X2(t), [unit=u"mol/m^3"]
            Z2(t), [unit=u"mol/m^3"]
            X3(t), [unit=u"mol/m^3"]
            Y3(t), [unit=u"mol/m^3"]
            Z3(t), [unit=u"mol/m^3"]
        end
        @parameters begin
            k1, [unit=u"(m^6)/(s*mol^2)"]
            v2, [unit=u"(m^6)/(s*mol^2)"]
            K2, [unit=u"mol/m^3"]
            v3, [unit=u"(m^3)/(s*mol)"]
            K3, [unit=u"mol/m^3"]
            n3
        end
        k1*X1, 2X1 --> Z1
        mm(X2, v2, K2), 3X2 --> Z2
        hill(X3, v3, K3, n3), X3 + Y3--> Z3
    end
end


### Jumps Unit Validation Tests ###

# Tests that user-provided jumps with correct units pass validation.
let
    using ModelingToolkitBase: ConstantRateJump, VariableRateJump, MassActionJump, Pre

    t = default_t()
    @independent_variables t [unit=u"s"]
    @species S(t) [unit=u"mol/m^3"]
    @parameters k [unit=u"mol/(m^3*s)"]

    # ConstantRateJump with correct rate units (species/time).
    rate = k
    affect = [S ~ Pre(S) + 1]
    user_jump = ConstantRateJump(rate, affect)

    rxs = [Reaction(k, nothing, [S])]
    @named rs = ReactionSystem(rxs, t, [S], [k]; jumps = [user_jump])

    # Should validate without warnings.
    @test validate(rs) == true
end

# Tests that ConstantRateJump with incorrect rate units issues a warning.
let
    using ModelingToolkitBase: ConstantRateJump, Pre

    @independent_variables t [unit=u"s"]
    @species S(t) [unit=u"mol/m^3"]
    @parameters k_bad [unit=u"mol/m^3"] k_good [unit=u"mol/(m^3*s)"]

    # ConstantRateJump with wrong rate units (missing /s).
    rate = k_bad
    affect = [S ~ Pre(S) + 1]
    user_jump = ConstantRateJump(rate, affect)

    rxs = [Reaction(k_good, nothing, [S])]
    @named rs = ReactionSystem(rxs, t, [S], [k_bad, k_good]; jumps = [user_jump])

    # Should issue warning and return false.
    @test (@test_logs (:warn, ) match_mode=:any validate(rs)) == false
end

# Tests that VariableRateJump with correct units passes validation.
let
    using ModelingToolkitBase: VariableRateJump, Pre

    @independent_variables t [unit=u"s"]
    @species S(t) [unit=u"mol/m^3"]
    @parameters k [unit=u"mol/(m^3*s)"]

    # VariableRateJump with correct rate units.
    rate = k * S / S  # Still has mol/(m^3*s) units
    affect = [S ~ Pre(S) + 1]
    user_jump = VariableRateJump(rate, affect)

    rxs = [Reaction(k, nothing, [S])]
    @named rs = ReactionSystem(rxs, t, [S], [k]; jumps = [user_jump])

    # Should validate without warnings.
    @test validate(rs) == true
end

# Tests that VariableRateJump with incorrect rate units issues a warning.
let
    using ModelingToolkitBase: VariableRateJump, Pre

    @independent_variables t [unit=u"s"]
    @species S(t) [unit=u"mol/m^3"]
    @parameters k_bad [unit=u"kg"] k_good [unit=u"mol/(m^3*s)"]

    # VariableRateJump with wrong rate units.
    rate = k_bad
    affect = [S ~ Pre(S) + 1]
    user_jump = VariableRateJump(rate, affect)

    rxs = [Reaction(k_good, nothing, [S])]
    @named rs = ReactionSystem(rxs, t, [S], [k_bad, k_good]; jumps = [user_jump])

    # Should issue warning and return false.
    @test (@test_logs (:warn, ) match_mode=:any validate(rs)) == false
end


### Events Unit Validation Tests ###

# Tests that discrete events with consistent units pass validation.
let
    using ModelingToolkitBase: SymbolicDiscreteCallback

    @independent_variables t [unit=u"s"]
    @species S(t) [unit=u"mol/m^3"]
    @parameters k [unit=u"mol/(m^3*s)"]

    # Discrete event: reset S every 1 second.
    discrete_events = [1.0 => [S ~ 0.0u"mol/m^3"]]

    rxs = [Reaction(k, nothing, [S])]
    @named rs = ReactionSystem(rxs, t; discrete_events)

    # Should validate without warnings (unitless time trigger is fine).
    @test validate(rs) == true
end

# Tests that continuous events with consistent units pass validation.
let
    using ModelingToolkitBase: SymbolicContinuousCallback

    @independent_variables t [unit=u"s"]
    @species S(t) [unit=u"mol/m^3"]
    @parameters k [unit=u"mol/(m^3*s)"] S_threshold [unit=u"mol/m^3"]

    # Continuous event: when S reaches threshold, reset to 0.
    continuous_events = SymbolicContinuousCallback([S ~ S_threshold] => [S ~ 0.0u"mol/m^3"])

    rxs = [Reaction(k, nothing, [S])]
    @named rs = ReactionSystem(rxs, t; continuous_events)

    # Should validate without warnings.
    @test validate(rs) == true
end
