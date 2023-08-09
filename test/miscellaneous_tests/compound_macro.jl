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
    @test isequal([C => 6, H => 12, O => 2], Catalyst.component_coefficients(C6H12O2)) 
    @test all(!iscompound(i) for i in components(C6H12O2))
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
    @test isequal([O => 2], Catalyst.component_coefficients(O2)) 
    @test all(!iscompound(i) for i in components(O2))
end

# Checks that compounds cannot be created from non-existing species.
let
    @variables t
    @species C(t) H(t)
    @test_throws Exception @compound C6H12O2(t) 6C 12H 2O    
end
let
    @variables t
    @test_throws Exception @compound O2(t) 2O    
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

# Checks that interpolation works.
let
    @variables t
    @species C(t) H(t) O(t)
    s = C
    @compound C6H12O2_1(t) 6s 12H 2O
    @compound C6H12O2_2(t) 6C 12H 2O

    @test iscompound(C6H12O2_1)
    @test iscompound(C6H12O2_2)

    @test isequal(components(C6H12O2_1), components(C6H12O2_2))
    @test isequal(coefficients(C6H12O2_1), coefficients(C6H12O2_2))
end

let
    @variables t
    @species C(t) H(t)
    @compound Cyclopentadiene(t) 5C 6H
    C5H6 = Cyclopentadiene
    @compound C10H12(t) 2C5H6

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
    @compound H2_1(t) alpha*H
    @compound H2_2(t) 2H

    @test iscompound(H2_1)
    @test iscompound(H2_2)

    @test isequal(components(H2_1),components(H2_2))
    @test isequal(coefficients(H2_1),coefficients(H2_2))
end

let
    @variables t
    @parameters alpha = 2
    @species H(t)

    @compound H2_1(t) alpha*H
    @compound H2_2(t) 2H

    @test iscompound(H2_1)
    @test iscompound(H2_2)

    @test isequal(components(H2_1),components(H2_2))
    @test isequal(coefficients(H2_1), @parameters alpha = 2)
end

let 
    @variables t
    @species A(t)
    B = A
    @compound A2(t) 2A
    @compound B2(t) 2B

    @test iscompound(A2)
    @test iscompound(B2)

    @test isequal(components(A2),components(B2))
    @test isequal(coefficients(A2), coefficients(B2))
end

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
    @test isequal(balanced_rx, balance_reaction(rx))

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

    balanced_rx = Reaction(k,[CO2,H2O],[C6H12O6,O2],[6, 6], [1, 6])
    @test isequal(balanced_rx, balance_reaction(rx))
end

# @reaction k, H2O --> H2O 
let 
    @variables t
    @species H(t) O(t)
    @compound H2O(t) 2H O

    rx = Reaction(1.0, [H2O], [H2O], [2], [2])

    balanced_rx = Reaction(1.0, [H2O], [H2O], [1], [1])
    @test isequal(balanced_rx, balance_reaction(rx))
end

# @reaction k, 2H + 1O --> H2O   
let 
    @variables t
    @species H(t) O(t)
    @compound H2O(t) 2H O

    rx = Reaction(1.0, [H, O], [H2O], [23, 1], [7])

    balanced_rx = Reaction(1.0, [H,O], [H2O], [2, 1], [1])
    @test isequal(balanced_rx, balance_reaction(rx))
end

# @reaction k, 1CH4 + 2O2 --> 1CO2 + 2H2O
let 
    @variables t
    @species H(t) O(t) C(t)
    @compound CH4(t) C 4H
    @compound O2(t) 2O
    @compound CO2(t) C 2O
    @compound H2O(t) 2H O

    rx = Reaction(1.0, [CH4, O2], [CO2, H2O])

    balanced_rx = Reaction(1.0, [CH4, O2], [CO2, H2O], [1, 2], [1, 2])
    @test isequal(balanced_rx, balance_reaction(rx))
end

# @reaction k, N2 + 3H2 --> 2NH3
let
    @variables t
    @species H(t) N(t)
    @compound N2(t) 2N
    @compound H2(t) 2H
    @compound NH3(t) N 3H

    rx = Reaction(1.0, [N2, H2], [NH3])

    balanced_rx = Reaction(1.0, [N2, H2], [NH3], [1 ,3], [2])
    @test isequal(balanced_rx, balance_reaction(rx))
end

# @reaction k, C2H5OH + CH3COOH --> C4H8O2 + H2O 
let
    @variables t
    @species C(t) H(t) O(t)
    @compound C2H5OH(t) 2C 6H O
    @compound CH3COOH(t) 2C 4H 2O
    @compound C4H8O2(t) 4C 8H 2O
    @compound H2O(t) 2H O

    rx = Reaction(1.0, [C2H5OH, CH3COOH], [C4H8O2, H2O])

    balanced_rx = Reaction(1.0, [C2H5OH, CH3COOH], [C4H8O2, H2O], [1, 1], [1, 1])
    @test isequal(balanced_rx, balance_reaction(rx))
end

# @reaction k, 2Ca3PO42 --> 6CaO + 1P4O10  
let 
    @variables t
    @species Ca(t) P(t) O(t)   
    @compound Ca3PO42(t) 3Ca 2P 8O
    @compound CaO(t) Ca O
    @compound P4O10(t) 4P 10O

    rx = Reaction(1.0, [Ca3PO42], [CaO, P4O10])

    balanced_rx = Reaction(1.0, [Ca3PO42], [CaO, P4O10], [2], [6, 1])
    @test isequal(balanced_rx, balance_reaction(rx))
end

# @reaction k, 4Fe + 3O2 + 6H2O --> 4FeOH3
let
    @variables t
    @species Fe(t) O(t) H(t)
    @compound O2(t) 2O 
    @compound H2O(t) 2H O
    @compound FeOH3(t) Fe 3H 3O

    rx = Reaction(1.0, [Fe, O2, H2O], [FeOH3])

    balanced_rx = Reaction(1.0, [Fe, O2, H2O], [FeOH3], [4, 3, 6], [4])
    @test isequal(balanced_rx, balance_reaction(rx))
end

# @reaction k, 2NaOH + H2SO4 --> Na2SO4 + 2H2O
let
    @variables t
    @species Na(t) O(t) H(t) S(t)
    @compound SO4(t) S 4O
    @compound NaOH(t) Na O H
    @compound H2SO4(t) 2H 1S 4O
    @compound Na2SO4(t) 2Na 1S 4O
    @compound H2O(t) 2H O

    rx = Reaction(1.0, [NaOH,H2SO4], [Na2SO4,H2O])

    balanced_rx = Reaction(1.0, [NaOH,H2SO4], [Na2SO4,H2O], [2, 1], [1, 2])
    @test isequal(balanced_rx, balance_reaction(rx))
end

# @reaction k, 2NO2 --> 1N2O4 
let
    @variables t
    @species N(t) O(t)
    @compound NO2(t) N 2O
    @compound N2O4(t) 2N 4O

    rx = Reaction(1.0, [NO2], [N2O4])

    balanced_rx = Reaction(1.0, [NO2], [N2O4], [2], [1])
    @test isequal(balanced_rx, balance_reaction(rx))
end

# @reaction k, 1CaCO3 + 2HCl --> CaCl2  + H2O + CO2
let 
    @variables t
    @species C(t) H(t) O(t) Ca(t) Cl(t)
    @compound H2O(t) 2H 1O
    @compound CO2(t) 1C 2O
    @compound CaCO3(t) 1Ca 1C 3O
    @compound HCl(t) 1H 1Cl
    @compound CaCl2(t) 1Ca 2Cl

    rx = Reaction(1.0,[CaCO3,HCl],[CaCl2,CO2,H2O])
    balanced_rx = Reaction(1.0,[CaCO3,HCl],[CaCl2,CO2,H2O], [1, 2], [1, 1, 1])
    @test isequal(balanced_rx, balance_reaction(rx))
end
