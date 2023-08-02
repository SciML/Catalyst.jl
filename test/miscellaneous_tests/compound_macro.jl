using Catalyst, Test

# Test base funcationality in two cases.
let
    @variables t
    @species C(t) H(t) O(t)
    @compound C6H12O2(t) 6C 12H 2O

    @test iscompound(C6H12O2)
    @test isspecies(C6H12O2)
    @test !iscompound(C)
    @test !iscompound(H)
    @test !iscompound(O)

    @test isequal([C, H, O], components(C6H12O2))
    @test isequal([6, 12, 2], coefficients(C6H12O2))
    @test isequal([C => 6, H => 12, O => 2], component_coefficients(C6H12O2)) # Maybe chaneg "component_coefficients" name soem time.
    @test all(!iscompound(i) for i in components(C6H12O2))
    # Test threw exception
    # Expression: all(.!(iscompound.components(C6H12O2)))
    # type #iscompound has no field components
end

let
    @variables t
    @species O(t)
    @compound O2(t) 2O

    @test iscompound(O2)
    @test isspecies(O2)
    @test !iscompound(O)

    @test isequal([O], components(O2))
    @test isequal([2], coefficients(O2))
    @test isequal([O => 2], component_coefficients(O2)) # Maybe chaneg "component_coefficients" name soem time.
    @test all(!iscompound(i) for i in components(O2))
end

# Checks that compounds cannot be created from non-existing species.
let
    @variables t
    @species C(t) H(t)
    @test_throws Exception @compound C6H12O2(t) 6C 12H 2O    # @test_throws might behave weirdly here, if you think this fails as it should but the test comes out starneg, just report and we will figure it out.
end
let
    @variables t
    @test_throws Exception @compound O2(t) 2O    # @test_throws might behave weirdly here, if you think this fails as it should but the test comes out starneg, just report and we will figure it out.
end

# Checks that nested components works as expected.
let
    @variables t
    @species C(t) H(t) O(t)
    @compound OH(t) 1O 1H
    @compound C3H5OH3(t) 3C 5H 3OH

    @test !iscompound(O)
    @test !iscompound(H)
    @test iscompound(OH)
    @test iscompound(C3H5OH3)
    @test !(all(!iscompound(i) for i in components(C3H5OH3)))

    @test !iscompound(components(C3H5OH3)[1])
    @test !iscompound(components(C3H5OH3)[2])
    @test iscompound(components(C3H5OH3)[3])

    @test isequal([C, H, OH], components(C3H5OH3))
    @test isequal([O, H], components(components(C3H5OH3)[3]))
    @test isequal([3, 5, 3], coefficients(C3H5OH3))
end

# # Checks that interpolation works.
# let
#     @variables t
#     @species C(t) H(t) O(t)
#     s = C
#     @compound C6H12O2_1(t) 6s 12H 2O
#     @compound C6H12O2_2(t) 6C 12H 2O

#     @test isequal(components(C6H12O2_1), components(C6H12O2_2))
#     @test isequal(coefficients(C6H12O2_1), coefficients(C6H12O2_2))
# end

# let
#     @variables t
#     @species C(t) H(t)
#     @compound Cyclopentadiene(t) 5C 6H
#     C5H6 = Cyclopentadiene
#     @compound C10H12(t) 2C5H6

#     @test iscompound(C10H12)
#     @test iscompound(components(C10H12)[1])

#     @test isequal(components(C10H12)[1], C5H6)
#     @test isequal(components(C10H12)[1], Cyclopentadiene)
#     @test isequal(coefficients(C10H12)[1], 2)
# end

#Check that balancing works.
let 
    @variables t
    @parameters k
    @species H(t) O(t)
    @compound H2(t) 2H
    @compound O2(t) 2O
    @compound H2O(t) 2H 1O

    rx = Reaction(k,[H2,O2],[H2O])

    @test isequal(create_matrix(rx),[2 0 -2; 0 2 -1; 0 0 1;])
    @test isequal(get_stoich(rx),[2, 1, 2])
    
    balanced_rx = Reaction(k,[H2,O2],[H2O],[2,1],[2])
    @test isequal(balanced_rx, balance(rx))

end

let 
    @variables t
    @parameters k
    @species C(t) H(t) O(t) 
    @compound O2(t) 2O
    @compound CO2(t) 1C 2O
    @compound H2O(t) 2H 1O
    @compound C6H12O6(t) 6C 12H 6O

    rx = Reaction(k,[CO2,H2O],[C6H12O6,O2])

    @test isequal(create_matrix(rx),[ 1 0 -6 0; 2 1 -6 -2; 0 2 -12 0; 0 0 0 1;])
    @test isequal(get_stoich(rx),[6, 6, 1, 6])
    
    balanced_rx = Reaction(k,[CO2,H2O],[C6H12O6,O2],[6, 6], [1,6])
    @test isequal(balanced_rx, balance(rx))
end