### Fetch required packages and reaction networks ###
using DiffEqBiological, Test, UnPack

### Tests construction of empty reaction networks ###
empty_network_1 = @reaction_network
@unpack eqs,iv,states,ps,name,systems = empty_network_1
@test length(eqs) == 0
@test iv.name == :t
@test length(states) == 0
@test length(ps) == 0

empty_network_2 = @reaction_network p1 p2 p3 p4 p5
@unpack eqs,iv,states,ps,name,systems = empty_network_2
@test length(eqs) == 0
@test iv.name == :t
@test length(states) == 0
@test length(ps) == 5
@test all(getproperty.(ps,:name) .== [:p1,:p2,:p3,:p4,:p5])


### All features not added yet, more tests will be added below as they are added ###
