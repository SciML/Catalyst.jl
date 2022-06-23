using Catalyst, Unitful, Test
const MT = ModelingToolkit

@parameters α [unit = u"μM/s"] β [unit = u"s"^(-1)] γ [unit = u"μM*s"^(-1)]
@variables t [unit = u"s"] A(t) [unit = u"μM"] B(t) [unit = u"μM"] C(t) [unit = u"μM"]
rxs = [Reaction(α, nothing, [A]),
    Reaction(β, [A], [B]),
    Reaction(γ, [A, B], [B], [1, 1], [2])]
rs = ReactionSystem(rxs, t, [A, B, C], [α, β, γ], name = Symbol("unittester"))
@test_nowarn ReactionSystem(rxs, t, [A, B, C], [α, β, γ], name = Symbol("unittester"))

odeunit = u"μM/s"
#jumpunit = u"s^(-1)"
for rx in reactions(rs)
    @test MT.get_unit(oderatelaw(rx)) == odeunit

    # we don't currently convert units, so they will be the same as for ODEs
    @test MT.get_unit(jumpratelaw(rx)) == odeunit
end

@test_nowarn convert(ODESystem, rs)
@test_nowarn convert(SDESystem, rs)
@test_nowarn convert(JumpSystem, rs)

@parameters β [unit = u"M"]
rxs = [Reaction(α, nothing, [A]),
    Reaction(β, [A], [B]),
    Reaction(γ, [A, B], [B], [1, 1], [2])]
@test (@test_logs (:warn,) match_mode=:any ModelingToolkit.validate(ReactionSystem(rxs, t,
                                                                                   [
                                                                                       A,
                                                                                       B,
                                                                                       C,
                                                                                   ],
                                                                                   [
                                                                                       α,
                                                                                       β,
                                                                                       γ,
                                                                                   ],
                                                                                   name = Symbol("unittester")))) ==
      false

@parameters β [unit = u"s"^(-1)]
@variables B(t) [unit = u"M"]
rxs = [Reaction(α, nothing, [A]),
    Reaction(β, [A], [B]),
    Reaction(γ, [A, B], [B], [1, 1], [2])]
@test (@test_logs (:warn,) match_mode=:any ModelingToolkit.validate(ReactionSystem(rxs, t,
                                                                                   [
                                                                                       A,
                                                                                       B,
                                                                                       C,
                                                                                   ],
                                                                                   [
                                                                                       α,
                                                                                       β,
                                                                                       γ,
                                                                                   ],
                                                                                   name = Symbol("unittester")))) ==
      false

@variables B(t) [unit = u"μM"] D(t) [unit = u"M"]
badrx1 = Reaction(α, [A], [D])
@test (@test_logs (:warn,) match_mode=:any ModelingToolkit.validate(badrx1)) == false
badrx2 = Reaction(α, [A], [B, D])
@test (@test_logs (:warn,) match_mode=:any ModelingToolkit.validate(badrx2)) == false
badrx3 = Reaction(α, [A, D], [B])
@test (@test_logs (:warn,) match_mode=:any ModelingToolkit.validate(badrx3)) == false
badrx4 = Reaction(α + β, [A], [B])
@test (@test_logs (:warn,) match_mode=:any ModelingToolkit.validate(badrx4)) == false
