using Catalyst, Test

### Tests Main Macro Creation Forms ### 
let
    @variables t
    @species C(t) H(t) O(t) 
    @parameters p1 p2

    # Basic cases that should pass:
    @compound H2O_1 ~ 2H + O
    @compound H2O_2, [output=true] ~ 2H + O
    @compound (H2O_3 = 1.5) ~ 2H + O
    @compound (H2O_4 = 4, [output=true]) ~ 2H + O
    @compound (H2O_5 = p1, [output=true]) ~ 2H + p2*O
    @test iscompound(H2O_1)
    @test iscompound(H2O_2)
    @test iscompound(H2O_3)
    @test iscompound(H2O_4)
    @test iscompound(H2O_5)

    # Independent variable given.
    @test_throws LoadError @eval @compound H2O(t) ~ 2H + O
    @test_throws LoadError @eval @compound H2O(t), [output=true] ~ 2H + O
    @test_throws LoadError @eval @compound (H2O(t) = 1.5) ~ 2H + O
    @test_throws LoadError @eval @compound (H2O(t) = 4, [output=true]) ~ 2H + O
    @test_throws LoadError @eval @compound (H2O(t) = p1, [output=true]) ~ 2H + p2*O

    # Other errors.
    @test_throws LoadError @eval @compound H2O(t) = 2H + O
    @test_throws LoadError @eval @compound H2O(t), [output=true] = 2H + O
    @test_throws LoadError @eval @compound H2O(t) = 1.5 ~ 2H + O
    @test_throws LoadError @eval @compound H2O(t) = 4, [output=true] ~ 2H + O
    @test_throws LoadError @eval @compound H2O(t) = p1, [output=true] ~ 2H + p2*O

    # Compounds created in block notation.
    @compounds begin
        CO2_1 ~ 2H + O
    end
    @compounds begin
        CO2_2, [output=true] ~ 2H + O
        (CO2_3 = 1.5) ~ 2H + O
        (CO2_4 = 4, [output=true]) ~ 2H + O
        (CO2_5 = p1, [output=true]) ~ 2H + p2*O
    end
    @test iscompound(CO2_1)
    @test iscompound(CO2_2)
    @test iscompound(CO2_3)
    @test iscompound(CO2_4)
    @test iscompound(CO2_5)

    # Declares stuff in the DSL.
    rn = @reaction_network begin
        @species N(t) H(t) 
        @parameters p1 p2
        @compounds begin
            NH3_1 ~ N + 3H
            NH3_2, [output=true] ~ N + 3H
            (NH3_3 = 1.5) ~ N + 3H
            (NH3_4 = 4, [output=true]) ~ N + 3H
            (NH3_5 = p1, [output=true]) ~ N + p2*H
        end
    end
    @test iscompound(rn.NH3_1)
    @test iscompound(rn.NH3_2)
    @test iscompound(rn.NH3_3)
    @test iscompound(rn.NH3_4)
    @test iscompound(rn.NH3_5)
end

### Other Minor Tests ###

# Test base functionality in two cases.
let
    @variables t
    @species C(t) H(t) O(t)
    @compound C6H12O2 ~ 6C + 12H + 2O

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
    @compound O2 ~ 2O

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
    @test_throws Exception @compound C6H12O2 ~ 6C + 12H + 2O    
end
let
    @variables t
    @test_throws Exception @compound O2 ~ 2O    
end

# Checks that nested components works as expected.
let
    @variables t
    @species C(t) H(t) O(t)
    @compound OH ~ 1O + 1H
    @compound C3H5OH3 ~ 3C + 5H + 3OH

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
    @compound C6H12O2_1 ~ 6s + 12H + 2O
    @compound C6H12O2_2 ~ 6C + 12H + 2O

    @test iscompound(C6H12O2_1)
    @test iscompound(C6H12O2_2)

    @test isequal(components(C6H12O2_1), components(C6H12O2_2))
    @test isequal(coefficients(C6H12O2_1), coefficients(C6H12O2_2))
    @test isequal(component_coefficients(C6H12O2_1), component_coefficients(C6H12O2_2))
end

let
    @variables t
    @species C(t) H(t)
    @compound Cyclopentadiene ~ 5C + 6H
    C5H6 = Cyclopentadiene
    @compound C10H12 ~ 2C5H6

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
    h = H
    @compound H2_1 ~ 2*H
    @compound H2_2 ~ alpha*H
    @compound H2_3 ~ 2*h
    @compound H2_4 ~ alpha*H

    @test iscompound(H2_1)
    @test iscompound(H2_2)
    @test iscompound(H2_2)
    @test iscompound(H2_4)

    @test isequal(components(H2_1),components(H2_2))
    @test isequal(components(H2_2),components(H2_3))
    @test isequal(components(H2_3),components(H2_4))
    @test isequal(coefficients(H2_1),coefficients(H2_2))
    @test isequal(coefficients(H2_2),coefficients(H2_3))
    @test isequal(coefficients(H2_3),coefficients(H2_4))
end

let
    @variables t
    @parameters alpha = 2
    @species H(t)

    @compound H2_1 ~ alpha*H
    @compound H2_2 ~ 2H

    @test iscompound(H2_1)
    @test iscompound(H2_2)

    @test isequal(components(H2_1),components(H2_2))
    @test isequal(coefficients(H2_1), @parameters alpha = 2)
end

let 
    @variables t
    @species A(t)
    B = A
    @compound A2 ~ 2A
    @compound B2 ~ 2B

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
    @compound OH ~ 1O + 1H
    @compound C3H5OH3 ~ 3C + 5H + 3OH

    @compounds begin
        OH_alt ~ 1O + 1H
        C3H5OH3_alt ~ 3C + 5H + 3OH
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
        comp ~ s1 + s2 + 4s3
        comp_alt ~ s1 + s2_alt + 4s3_alt
    end

    @test iscompound(comp)
    @test iscompound(comp_alt)

    @test isequal(components(comp),components(comp_alt))
    @test isequal(coefficients(comp), coefficients(comp_alt))
    @test isequal(component_coefficients(comp), component_coefficients(comp_alt))
end

### Compounds in DSL ###

# Checks with a single compound.
# Checks using @unpack.
# Check where compounds and components does not occur in reactions.
let 
    rn = @reaction_network begin
        @species C(t) O(t)
        @compounds begin
            CO2 ~ C + 2O
        end
    end
    @unpack C, O, CO2 = rn
    
    @test length(species(rn)) == 3
    @test iscompound(CO2)
    @test isequal([C, O], components(CO2))
    @test isequal([1, 2], coefficients(CO2))
    @test isequal([C => 1, O => 2], component_coefficients(CO2)) 
end

# Test using multiple compounds.
# Test using rn. notation to fetch species.
let 
    rn = complete(@reaction_network begin
        @species C(t) O(t) H(t)
        @compounds begin
            CH4 ~ C + 4H
            O2 ~ 2O
            CO2 ~ C + 2O
            H2O ~ 2H + O
        end
        k, CH4 + O2 --> CO2 + H2O
    end)
    species(rn)
    
    @test length(species(rn)) == 7
    @test isequal([rn.C, rn.H], components(rn.CH4))
    @test isequal([1, 4], coefficients(rn.CH4))
    @test isequal([rn.C => 1, rn.H => 4], component_coefficients(rn.CH4)) 
    @test isequal([rn.O], components(rn.O2))
    @test isequal([2], coefficients(rn.O2))
    @test isequal([rn.O => 2], component_coefficients(rn.O2)) 
    @test isequal([rn.C, rn.O], components(rn.CO2))
    @test isequal([1, 2], coefficients(rn.CO2))
    @test isequal([rn.C => 1, rn.O => 2], component_coefficients(rn.CO2)) 
    @test isequal([rn.H, rn.O], components(rn.H2O))
    @test isequal([2, 1], coefficients(rn.H2O))
    @test isequal([rn.H => 2, rn.O => 1], component_coefficients(rn.H2O)) 
end

# Tests using compounds of compounds.
# Tests where species are part of reactions and not declared using "@species".
let
    rn = @reaction_network begin
        @compounds begin
            SO2 ~ S + 2O
            S2O4 ~ 2SO2
        end
        dS, S --> 0
        dO, O --> 0
    end
    species(rn)
    @unpack S, O, SO2, S2O4 = rn
    
    @test length(species(rn)) == 4
    
    @test isequal([S, O], components(SO2))
    @test isequal([1, 2], coefficients(SO2))
    @test isequal([S => 1, O => 2], component_coefficients(SO2)) 
    @test isequal([SO2], components(S2O4))
    @test isequal([2], coefficients(S2O4))
    @test isequal([SO2 => 2], component_coefficients(S2O4)) 
end