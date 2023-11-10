using Catalyst, Test

# Test base funcationality in two cases.
let
    @variables t
    @species C(t) H(t) + O(t)
    @compound C6H12O2(t) = 6C + 12H + 2O

    @test iscompound(C6H12O2)
    @test isspecies(C6H12O2)
    @test !iscompound(C)
    @test !iscompound(H)
    @test !iscompound(O)

    @test isequal([C, H, O], components(C6H12O2))
    @test isequal([6, 12, 2], coefficients(C6H12O2))
    @test isequal([C => 6, H => 12, O => 2], Catalyst.component_coefficients(C6H12O2)) 
    @test all(!iscompound(i) for i in components(C6H12O2))
end

let
    @variables t
    @species O(t)
    @compound O2(t) = 2O

    @test iscompound(O2)
    @test isspecies(O2)
    @test !iscompound(O)

    @test isequal([O], components(O2))
    @test isequal([2], coefficients(O2))
    @test isequal([O => 2], Catalyst.component_coefficients(O2)) 
    @test all(!iscompound(i) for i in components(O2))
end

# Checks that compounds cannot be created from non-existing species.
let
    @variables t
    @species C(t) H(t)
    @test_throws Exception @compound C6H12O2(t) = 6C + 12H + 2O    
end
let
    @variables t
    @test_throws Exception @compound O2(t) = 2O    
end

# Checks that nested components works as expected.
let
    @variables t
    @species C(t) H(t) O(t)
    @compound OH(t) = 1O + 1H
    @compound C3H5OH3(t) = 3C + 5H + 3OH

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

# Checks that interpolation works.
let
    @variables t
    @species C(t) H(t) O(t)
    s = C
    @compound C6H12O2_1(t) = 6s + 12H + 2O
    @compound C6H12O2_2(t) = 6C + 12H + 2O

    @test iscompound(C6H12O2_1)
    @test iscompound(C6H12O2_2)

    @test isequal(components(C6H12O2_1), components(C6H12O2_2))
    @test isequal(coefficients(C6H12O2_1), coefficients(C6H12O2_2))
    @test isequal(component_coefficients(C6H12O2_1), component_coefficients(C6H12O2_2))
end

let
    @variables t
    @species C(t) H(t)
    @compound Cyclopentadiene(t) = 5C + 6H
    C5H6 = Cyclopentadiene
    @compound C10H12(t) = 2C5H6

    @test iscompound(C10H12)
    @test iscompound(components(C10H12)[1])

    @test isequal(components(C10H12)[1], C5H6)
    @test isequal(components(C10H12)[1], Cyclopentadiene)
    @test isequal(coefficients(C10H12)[1], 2)
end

let
    @variables t
    @species H(t)

    alpha = 2
    @compound H2_1(t) = alpha*H
    @compound H2_2(t) = 2H

    @test iscompound(H2_1)
    @test iscompound(H2_2)

    @test isequal(components(H2_1),components(H2_2))
    @test isequal(coefficients(H2_1),coefficients(H2_2))
end

let
    @variables t
    @parameters alpha = 2
    @species H(t)

    @compound H2_1(t) = alpha*H
    @compound H2_2(t) = 2H

    @test iscompound(H2_1)
    @test iscompound(H2_2)

    @test isequal(components(H2_1),components(H2_2))
    @test isequal(coefficients(H2_1), @parameters alpha = 2)
end

let 
    @variables t
    @species A(t)
    B = A
    @compound A2(t) = 2A
    @compound B2(t) = 2B

    @test iscompound(A2)
    @test iscompound(B2)

    @test isequal(components(A2),components(B2))
    @test isequal(coefficients(A2), coefficients(B2))
    @test isequal(component_coefficients(A2), component_coefficients(B2))
end

###  Check @compounds Macro ###

# Basic syntax.
let 
    @variables t
    @species C(t) H(t) O(t)
    @compound OH(t) = 1O + 1H
    @compound C3H5OH3(t) = 3C + 5H + 3OH

    @compounds begin
        OH_alt(t) = 1O + 1H
        C3H5OH3_alt(t) = 3C + 5H + 3OH
    end

    @test iscompound(OH_alt)
    @test iscompound(C3H5OH3_alt)

    @test isequal(components(OH),components(OH_alt))
    @test isequal(coefficients(OH), coefficients(OH_alt))
    @test isequal(component_coefficients(OH), component_coefficients(OH_alt))
    @test isequal(components(C3H5OH3),components(C3H5OH3_alt))
    @test isequal(coefficients(C3H5OH3), coefficients(C3H5OH3_alt))
    @test isequal(component_coefficients(C3H5OH3), component_coefficients(C3H5OH3_alt))
end

# Interpolation
let 
    @variables t
    @species s1(t) s2(t) s3(t)
    s2_alt = s2
    s3_alt = s3

    @compounds begin
        comp(t) = s1 + s2 + 4s3
        comp_alt(t) = s1 + s2_alt + 4s3_alt
    end

    @test iscompound(comp)
    @test iscompound(comp_alt)

    @test isequal(components(comp),components(comp_alt))
    @test isequal(coefficients(comp), coefficients(comp_alt))
    @test isequal(component_coefficients(comp), component_coefficients(comp_alt))
end