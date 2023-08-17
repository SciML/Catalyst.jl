using Catalyst, Test

#Check that balancing works.
let 
    @variables t
    @parameters k
    @species H(t) O(t)
    @compound H2(t) 2H
    @compound O2(t) 2O
    @compound H2O(t) 2H 1O

    rx = Reaction(k,[H2,O2],[H2O])

    @test isequal(Catalyst.create_matrix(rx),[2 0 -2; 0 2 -1;])
    @test isequal(Catalyst.get_stoich(rx),[2, 1, 2])

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

    @test isequal(Catalyst.create_matrix(rx),[ 1 0 -6 0; 2 1 -6 -2; 0 2 -12 0;])
    @test isequal(Catalyst.get_stoich(rx),[6, 6, 1, 6])

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
