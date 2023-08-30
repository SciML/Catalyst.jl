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

# @reaction k, SiCl4 + 4H2O → H4SiO4 + 4HCl
let 
    @variables t
    @species Si(t) Cl(t) H(t) O(t)
    @compound SiCl4(t) 1Si 4Cl
    @compound H2O(t) 2H O
    @compound H4SiO4(t) 4H Si 4O
    @compound HCl(t) H Cl

    rx = Reaction(1.0,[SiCl4,H2O],[H4SiO4,HCl])
    balanced_rx = Reaction(1.0,[SiCl4,H2O],[H4SiO4,HCl], [1,4], [1,4])
    @test isequal(balanced_rx, balance_reaction(rx))
end

# @reaction k, 2Al + 6HCl → 2AlCl3 + 3H2
let 
    @variables t
    @species Al(t) Cl(t) H(t)
    @compound HCl(t) H Cl
    @compound AlCl3(t) Al 3Cl
    @compound H2(t) 2H

    rx = Reaction(1.0,[Al,HCl],[AlCl3,H2])
    balanced_rx = Reaction(1.0,[Al,HCl],[AlCl3,H2],[2,6], [2,3])
    @test isequal(balanced_rx, balance_reaction(rx))
end

# @reaction k, Na2CO3 + 2HCl → 2NaCl + H2O + CO2
let 
    @variables t
    @species Na(t) C(t) O(t) H(t) Cl(t)
    @compound Na2CO3(t) 2Na C 3O
    @compound HCl(t) H Cl
    @compound NaCl(t) Na Cl
    @compound H2O(t) 2H O
    @compound CO2(t) C 2O

    rx = Reaction(1.0,[Na2CO3,HCl],[NaCl,H2O,CO2])
    balanced_rx = Reaction(1.0,[Na2CO3,HCl],[NaCl,H2O,CO2], [1,2], [2,1,1])
    @test isequal(balanced_rx, balance_reaction(rx))
end

# @reaction k, 2C7H6O2 + 15O2 → 14CO2 + 6H2O
let 
    @variables t
    @species C(t) H(t) O(t)
    @compound C7H6O2(t) 7C 6H 2O
    @compound O2(t) 2O
    @compound CO2(t) C 2O
    @compound H2O(t) 2H O

    rx = Reaction(1.0,[C7H6O2,O2],[CO2,H2O])
    balanced_rx = Reaction(1.0,[C7H6O2,O2],[CO2,H2O], [2,15], [14,6])
    @test isequal(balanced_rx, balance_reaction(rx))
end

# @reaction k,Fe2(SO4)3 + 6KOH → 3K2SO4 + 2Fe(OH)3
let 
    @variables t
    @species Fe(t) S(t) O(t) H(t) K(t)
    @compound Fe2S3O12(t) 2Fe 3S 12O
    @compound KOH(t) K O H
    @compound K2SO4(t) 2K S 4O
    @compound FeO3H3(t) Fe 3O 3H

    rx = Reaction(1.0,[Fe2S3O12,KOH],[K2SO4,FeO3H3]) #5x4 matrix
    balanced_rx = Reaction(1.0,[Fe2S3O12,KOH],[K2SO4,FeO3H3], [1,6], [3,2])
    @test isequal(balanced_rx, balance_reaction(rx))
end

# @reaction k, 2Ca3(PO4)2 + 6SiO2 → P4O10 + 6CaSiO3
let 
    @variables t
    @species Ca(t) P(t) O(t) Si(t)
    @compound Ca3P2O8(t) 3Ca 2P 8O
    @compound SiO2(t) Si 2O
    @compound P4O10(t) 4P 10O
    @compound CaSiO3(t) Ca Si 3O

    rx = Reaction(1.0,[Ca3P2O8,SiO2],[P4O10,CaSiO3]) #5x4 matrix
    balanced_rx = Reaction(1.0,[Ca3P2O8,SiO2],[P4O10,CaSiO3], [2,6] , [1,6])
    @test isequal(balanced_rx, balance_reaction(rx))
end

# @reaction k, 4KClO3 → 3KClO4 + KCl
let 
    @variables t
    @species K(t) Cl(t) O(t)
    @compound KClO3(t) K Cl 3O
    @compound KClO4(t) K Cl  4O
    @compound KCl(t) K Cl

    rx = Reaction(1.0,[KClO3],[KClO4,KCl])
    balanced_rx = Reaction(1.0,[KClO3],[KClO4,KCl], [4], [3,1]) 
    @test isequal(balanced_rx, balance_reaction(rx))
end

# @reaction k, Al2(SO4)3 + 3Ca(OH)2 → 2Al(OH)3 + 3CaSO4
let 
    @variables t
    @species Al(t) S(t) O(t) Ca(t) O(t) (H)
    @compound Al2S3O12(t) 2Al 3S 12O
    @compound CaO2H2(t) Ca 2O 2H
    @compound AlO3H3(t) Al 3O 3H
    @compound CaSO4(t) Ca S 4O

    rx = Reaction(1.0,[Al2S3O12,CaO2H2],[AlO3H3,CaSO4]) 
    balanced_rx = Reaction(1.0,[Al2S3O12,CaO2H2],[AlO3H3,CaSO4], [1,3], [2,3])
    @test isequal(balanced_rx, balance_reaction(rx))
end

# @reaction k, H2SO4 + 8HI → H2S + 4I2 + 4H2O
let 
    @variables t
    @species H(t) S(t) O(t) I(t)
    @compound H2SO4(t) 2H S 4O
    @compound HI(t) H I
    @compound H2S(t) 2H S
    @compound I2(t) 2I
    @compound H2O(t) 2H O

    rx = Reaction(1.0,[H2SO4,HI],[H2S,I2,H2O]) 
    balanced_rx = Reaction(1.0,[H2SO4,HI],[H2S,I2,H2O], [1,8], [1,4,4]) 
    @test isequal(balanced_rx, balance_reaction(rx))
end

# @reaction k, C2H4 + 3O2 ↔ 2CO2 + 2H2O
let 
    @variables t
    @species C(t) H(t) O(t)
    @compound C2H4(t) 2C 4H
    @compound O2(t) 2O
    @compound CO2(t) C 2O
    @compound H2O(t) 2H O

    rx = Reaction(1.0,[C2H4,O2],[CO2,H2O])
    balanced_rx = Reaction(1.0,[C2H4,O2],[CO2,H2O],[1,3],[2,2])
    @test isequal(balanced_rx, balance_reaction(rx))
end

# Infinite solutions
let 
    @variables t
    @species C(t) H(t) O(t)
    @compound CO(t) C O
    @compound CO2(t) C 2O
    @compound H2(t) 2H
    @compound CH4(t) C 4H
    @compound H2O(t) 2H O

    rx = Reaction(1.0,[CO,CO2,H2],[CH4,H2O])
    @test typeof(balance_reaction(rx)) == Vector{Reaction}
    @test length(balance_reaction(rx)) == 2
end

# No way to balance
let 
    @variables t
    @species Fe(t) S(t) O(t) H(t) N(t)

    @compound FeS2(t) Fe 2S
    @compound HNO3(t) H N 3O
    @compound Fe2S3O12(t) 2Fe 3S 12O
    @compound NO(t) N O
    @compound H2SO4(t) 2H S 4O

    rx = Reaction(1.0,[FeS2,HNO3],[Fe2S3O12,NO,H2SO4])
    @test_throws ErrorException("Chemical equation cannot be balanced") balance_reaction(rx)
end
