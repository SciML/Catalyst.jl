### Prepares Tests ###

# Fetch packages.
using Catalyst, DynamicQuantities, Test
using ModelingToolkit: get_iv, get_unit, validate, ValidationError


### Basic Tests ###

# Checks that units work with programmatic model creation.
# μM and M units currently do not work, so am using other units in the meantime.
let
    # Creates a `ReactionSystem` programmatically, while designating units.
    @variables t [unit=u"s"]
    @species A(t) [unit=u"m"] B(t) [unit=u"m"] C(t) [unit=u"m"]
    @parameters k1 [unit=u"m/s"] k2 [unit=u"s"^(-1)] k3 [unit=u"m*s"^(-1)]
    rxs = [Reaction(k1, nothing, [A]),
        Reaction(k2, [A], [B]),
        Reaction(k3, [A, B], [B], [1, 1], [2])]
    @test_nowarn @named rs = ReactionSystem(rxs, t, [A, B, C], [k1, k2, k3])
    rs = complete(rs)

    # Test that all reactions have the correct unit.
    for rx in reactions(rs)
        @test get_unit(oderatelaw(rx)) == u"m/s"
        # we don't currently convert units, so they will be the same as for ODEs
        @test get_unit(jumpratelaw(rx)) == u"m/s"
    end

    # Tests that the system can be converted to MTK systems without warnings.
    @test_nowarn convert(ODESystem, rs)
    @test_nowarn convert(SDESystem, rs)
    @test_nowarn convert(JumpSystem, rs)
    @test_nowarn convert(NonlinearSystem, rs)

    # Tests that creating `Reaction`s with non-matching units yields warnings.
    @species B(t) [unit=u"m"] D(t) [unit=u"A"]
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
    @parameters k2 [unit=u"A"]
    bad_rxs = [
        Reaction(k2, nothing, [A]),
        Reaction(k2, [A], [B]),
        Reaction(k3, [A, B], [B], [1, 1], [2])
    ]
    @named bad_rs = ReactionSystem(bad_rxs, t, [A, B, C], [k1, k2, k3])
    @test (@test_logs (:warn, ) match_mode=:any validate(bad_rs)) == false
end


# Checks that units work in the DSL.
begin
    # Creates the model, while designating matching units.
    rs = @reaction_network begin
        @ivs t [unit=u"s"]
        @species begin
            A(t), [unit=u"m"]
            B(t), [unit=u"m"]
            C(t), [unit=u"m"]
        end
        @parameters begin
            k1, [unit=u"m/s"]
            k2, [unit=u"s"^(-1)]
            k3, [unit=u"m*s"^(-1)]
        end
        k1, 0 --> A
        k2, A --> B
        k3, A + B --> 2B
    end

    # Checks that the `ReactionSystem`'s content have the correct units.
    @test get_unit(get_iv(rs)) == u"s"
    @test all(get_unit.([rs.A, rs.B, rs.C]) .== [u"m", u"m", u"m"])
    @test all(get_unit.([rs.k1, rs.k2, rs.k3]) .== [u"m/s", u"s"^(-1), u"m*s"^(-1)])
    for rx in reactions(rs)
        @test get_unit(oderatelaw(rx)) == u"m/s"
        # we don't currently convert units, so they will be the same as for ODEs
        @test get_unit(jumpratelaw(rx)) == u"m/s"
    end

    # Checks that system declarations with erroneous errors yields errors.
    @test_logs (:warn, ) match_mode=:any @reaction_network begin
        @ivs t [unit=u"1/s"]
        @species begin
            A(t), [unit=u"m"] 
            B(t), [unit=u"m"]
            C(t), [unit=u"m"]
        end
        @parameters begin
            k1, [unit=u"m/s"]
            k2, [unit=u"s"^(-1)]
            k3, [unit=u"m*s"^(-1)]
        end
        k1, 0 --> A
        k2, A --> B
        k3, A + B --> 2B
    end
    @test_logs (:warn, ) match_mode=:any @reaction_network begin
        @ivs t [unit=u"s"]
        @species begin
            A(t), [unit=u"m"]
            B(t), [unit=u"m"]
            C(t), [unit=u"m"]
        end
        @parameters begin
            k1, [unit=u"m"]
            k2, [unit=u"s"^(-1)]
            k3, [unit=u"m*s"^(-1)]
        end
        k1, 0 --> A
        k2, A --> B
        k3, A + B --> 2B
    end
    @test_logs (:warn, ) match_mode=:any @reaction_network begin
        @ivs t [unit=u"s"]
        @species begin
            A(t), [unit=u"m*s"]
            B(t), [unit=u"m"]
            C(t), [unit=u"m"]
        end
        @parameters begin
            k1, [unit=u"m/s"] 
            k2, [unit=u"s"^(-1)] 
            k3, [unit=u"m*s"^(-1)]
        end
        k1, 0 --> A
        k2, A --> B
        k3, A + B --> 2B
    end
end

