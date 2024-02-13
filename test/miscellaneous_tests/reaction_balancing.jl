using Catalyst, Test
import Catalyst: t_nounits as t

#Check that balancing works.
let
    @parameters k
    @species H(t) O(t)
    @compounds begin
        H2 ~ 2H
        O2 ~ 2O
        H2O ~ 2H + 1O
    end

    rx = Reaction(k,[H2,O2],[H2O])

    @test isequal(Catalyst.create_matrix(rx),[2 0 -2; 0 2 -1;])
    @test isequal(Catalyst.get_balanced_stoich(rx),[[2, 1, 2]])

    balanced_rx = Reaction(k,[H2,O2],[H2O],[2,1],[2])
    brxs = balance_reaction(rx)
    @test length(brxs) == 1
    @test isequal(balanced_rx, first(brxs))

    rx_macro = @reaction k, $H2 + $O2 --> $H2O
    brxs_macro = balance_reaction(rx_macro)
    @test isequal(brxs, brxs_macro)
end

let
    @parameters k
    @species C(t) H(t) O(t)
    @compound O2 ~ 2O
    @compound CO2 ~ 1C + 2O
    @compound H2O ~ 2H + 1O
    @compound C6H12O6 ~ 6C + 12H + 6O

    rx = Reaction(k,[CO2,H2O],[C6H12O6,O2])

    @test isequal(Catalyst.create_matrix(rx),[ 1 0 -6 0; 2 1 -6 -2; 0 2 -12 0;])
    @test isequal(Catalyst.get_balanced_stoich(rx),[[6, 6, 1, 6]])

    balanced_rx = Reaction(k,[CO2,H2O],[C6H12O6,O2],[6, 6], [1, 6])
    brxs = balance_reaction(rx)
    @test length(brxs) == 1
    @test isequal(balanced_rx, first(brxs))

    rx_macro = @reaction k, $CO2 + $H2O --> $C6H12O6 + $O2
    brxs_macro = balance_reaction(rx_macro)
    @test isequal(brxs, brxs_macro)
end

# @reaction k, H2O --> H2O
let
    @species H(t) O(t)
    @compound H2O ~ 2H + O

    rx = Reaction(1.0, [H2O], [H2O], [2], [2])

    balanced_rx = Reaction(1.0, [H2O], [H2O], [1], [1])
    brxs = balance_reaction(rx)
    @test length(brxs) == 1
    @test isequal(balanced_rx, first(brxs))
end

# @reaction k, 2H + 1O --> H2O
let
    @species H(t) O(t)
    @compound H2O ~ 2H + O

    rx = Reaction(1.0, [H, O], [H2O], [23, 1], [7])

    balanced_rx = Reaction(1.0, [H,O], [H2O], [2, 1], [1])
    brxs = balance_reaction(rx)
    @test length(brxs) == 1
    @test isequal(balanced_rx, first(brxs))
end

# @reaction k, 1CH4 + 2O2 --> 1CO2 + 2H2O
let
    @species H(t) O(t) C(t)
    @compounds begin
        CH4 ~ C + 4H
        O2 ~ 2O
        CO2 ~ C + 2O
        H2O ~ 2H + O
    end

    rx = Reaction(1.0, [CH4, O2], [CO2, H2O])

    balanced_rx = Reaction(1.0, [CH4, O2], [CO2, H2O], [1, 2], [1, 2])
    brxs = balance_reaction(rx)
    @test length(brxs) == 1
    @test isequal(balanced_rx, first(brxs))
end

# @reaction k, N2 + 3H2 --> 2NH3
let
    @species H(t) N(t)
    @compound N2 ~ 2N
    @compound H2 ~ 2H
    @compound NH3 ~ N + 3H

    rx = Reaction(1.0, [N2, H2], [NH3])

    balanced_rx = Reaction(1.0, [N2, H2], [NH3], [1 ,3], [2])
    brxs = balance_reaction(rx)
    @test length(brxs) == 1
    @test isequal(balanced_rx, first(brxs))
end

# @reaction k, C2H5OH + CH3COOH --> C4H8O2 + H2O
let
    @species C(t) H(t) O(t)
    @compound C2H5OH ~ 2C + 6H + O
    @compound CH3COOH ~ 2C + 4H + 2O
    @compound C4H8O2 ~ 4C + 8H + 2O
    @compound H2O ~ 2H + O

    rx = Reaction(1.0, [C2H5OH, CH3COOH], [C4H8O2, H2O])

    balanced_rx = Reaction(1.0, [C2H5OH, CH3COOH], [C4H8O2, H2O], [1, 1], [1, 1])
    brxs = balance_reaction(rx)
    @test length(brxs) == 1
    @test isequal(balanced_rx, first(brxs))
end

# @reaction k, 2Ca3PO42 --> 6CaO + 1P4O10
let
    @species Ca(t) P(t) O(t)
    @compound Ca3PO42 ~ 3Ca + 2P + 8O
    @compound CaO ~ Ca + O
    @compound P4O10 ~ 4P + 10O

    rx = Reaction(1.0, [Ca3PO42], [CaO, P4O10])

    balanced_rx = Reaction(1.0, [Ca3PO42], [CaO, P4O10], [2], [6, 1])
    brxs = balance_reaction(rx)
    @test length(brxs) == 1
    @test isequal(balanced_rx, first(brxs))
end

# @reaction k, 4Fe + 3O2 + 6H2O --> 4FeOH3
let
    @species Fe(t) O(t) H(t)
    @compound O2 ~ 2O
    @compound H2O ~ 2H + O
    @compound FeOH3 ~ Fe + 3H + 3O

    rx = Reaction(1.0, [Fe, O2, H2O], [FeOH3])

    balanced_rx = Reaction(1.0, [Fe, O2, H2O], [FeOH3], [4, 3, 6], [4])
    brxs = balance_reaction(rx)
    @test length(brxs) == 1
    @test isequal(balanced_rx, first(brxs))
end

# @reaction k, 2NaOH + H2SO4 --> Na2SO4 + 2H2O
let
    @species Na(t) O(t) H(t) S(t)
    @compound SO4 ~ S + 4O
    @compound NaOH ~ Na + O + H
    @compound H2SO4 ~ 2H + 1S + 4O
    @compound Na2SO4 ~ 2Na + 1S + 4O
    @compound H2O ~ 2H + O

    rx = Reaction(1.0, [NaOH,H2SO4], [Na2SO4,H2O])

    balanced_rx = Reaction(1.0, [NaOH,H2SO4], [Na2SO4,H2O], [2, 1], [1, 2])
    brxs = balance_reaction(rx)
    @test length(brxs) == 1
    @test isequal(balanced_rx, first(brxs))
end

# @reaction k, 2NO2 --> 1N2O4
let
    @species N(t) O(t)
    @compound NO2 ~ N + 2O
    @compound N2O4 ~ 2N + 4O

    rx = Reaction(1.0, [NO2], [N2O4])

    balanced_rx = Reaction(1.0, [NO2], [N2O4], [2], [1])
    brxs = balance_reaction(rx)
    @test length(brxs) == 1
    @test isequal(balanced_rx, first(brxs))
end

# @reaction k, 1CaCO3 + 2HCl --> CaCl2  + H2O + CO2
let
    @species C(t) H(t) O(t) Ca(t) Cl(t)
    @compound H2O ~ 2H + 1O
    @compound CO2 ~ 1C + 2O
    @compound CaCO3 ~ 1Ca + 1C + 3O
    @compound HCl ~ 1H + 1Cl
    @compound CaCl2 ~ 1Ca + 2Cl

    rx = Reaction(1.0,[CaCO3,HCl],[CaCl2,CO2,H2O])
    balanced_rx = Reaction(1.0,[CaCO3,HCl],[CaCl2,CO2,H2O], [1, 2], [1, 1, 1])
    brxs = balance_reaction(rx)
    @test length(brxs) == 1
    @test isequal(balanced_rx, first(brxs))
end

# @reaction k, SiCl4 + 4H2O → H4SiO4 + 4HCl
let
    @species Si(t) Cl(t) H(t) O(t)
    @compound SiCl4 ~ 1Si + 4Cl
    @compound H2O ~ 2H + O
    @compound H4SiO4 ~ 4H + Si + 4O
    @compound HCl ~ H + Cl

    rx = Reaction(1.0,[SiCl4,H2O],[H4SiO4,HCl])
    balanced_rx = Reaction(1.0,[SiCl4,H2O],[H4SiO4,HCl], [1,4], [1,4])
    brxs = balance_reaction(rx)
    @test length(brxs) == 1
    @test isequal(balanced_rx, first(brxs))
end

# @reaction k, 2Al + 6HCl → 2AlCl3 + 3H2
let
    @species Al(t) Cl(t) H(t)
    @compound HCl ~ H + Cl
    @compound AlCl3 ~ Al + 3Cl
    @compound H2 ~ 2H

    rx = Reaction(1.0,[Al,HCl],[AlCl3,H2])
    balanced_rx = Reaction(1.0,[Al,HCl],[AlCl3,H2],[2,6], [2,3])
    brxs = balance_reaction(rx)
    @test length(brxs) == 1
    @test isequal(balanced_rx, first(brxs))
end

# @reaction k, Na2CO3 + 2HCl → 2NaCl + H2O + CO2
let
    @species Na(t) C(t) O(t) H(t) Cl(t)
    @compound Na2CO3 ~ 2Na + C + 3O
    @compound HCl ~ H + Cl
    @compound NaCl ~ Na + Cl
    @compound H2O ~ 2H + O
    @compound CO2 ~ C + 2O

    rx = Reaction(1.0,[Na2CO3,HCl],[NaCl,H2O,CO2])
    balanced_rx = Reaction(1.0,[Na2CO3,HCl],[NaCl,H2O,CO2], [1,2], [2,1,1])
    brxs = balance_reaction(rx)
    @test length(brxs) == 1
    @test isequal(balanced_rx, first(brxs))
end

# @reaction k, 2C7H6O2 + 15O2 → 14CO2 + 6H2O
let
    @species C(t) H(t) O(t)
    @compound C7H6O2 ~ 7C + 6H + 2O
    @compound O2 ~ 2O
    @compound CO2 ~ C + 2O
    @compound H2O ~ 2H + O

    rx = Reaction(1.0,[C7H6O2,O2],[CO2,H2O])
    balanced_rx = Reaction(1.0,[C7H6O2,O2],[CO2,H2O], [2,15], [14,6])
    brxs = balance_reaction(rx)
    @test length(brxs) == 1
    @test isequal(balanced_rx, first(brxs))
end

# @reaction k,Fe2(SO4)3 + 6KOH → 3K2SO4 + 2Fe(OH)3
let
    @species Fe(t) S(t) O(t) H(t) K(t)
    @compound Fe2S3O12 ~ 2Fe + 3S + 12O
    @compound KOH ~ K + O + H
    @compound K2SO4 ~ 2K + S + 4O
    @compound FeO3H3 ~ Fe + 3O + 3H

    rx = Reaction(1.0,[Fe2S3O12,KOH],[K2SO4,FeO3H3]) #5x4 matrix
    balanced_rx = Reaction(1.0,[Fe2S3O12,KOH],[K2SO4,FeO3H3], [1,6], [3,2])
    brxs = balance_reaction(rx)
    @test length(brxs) == 1
    @test isequal(balanced_rx, first(brxs))
end

# @reaction k, 2Ca3(PO4)2 + 6SiO2 → P4O10 + 6CaSiO3
let
    @species Ca(t) P(t) O(t) Si(t)
    @compound Ca3P2O8 ~ 3Ca + 2P + 8O
    @compound SiO2 ~ Si + 2O
    @compound P4O10 ~ 4P + 10O
    @compound CaSiO3 ~ Ca + Si + 3O

    rx = Reaction(1.0,[Ca3P2O8,SiO2],[P4O10,CaSiO3]) #5x4 matrix
    balanced_rx = Reaction(1.0,[Ca3P2O8,SiO2],[P4O10,CaSiO3], [2,6] , [1,6])
    brxs = balance_reaction(rx)
    @test length(brxs) == 1
    @test isequal(balanced_rx, first(brxs))
end

# @reaction k, 4KClO3 → 3KClO4 + KCl
let
    @species K(t) Cl(t) O(t)
    @compound KClO3 ~ K + Cl + 3O
    @compound KClO4 ~ K + Cl + 4O
    @compound KCl ~ K + Cl

    rx = Reaction(1.0,[KClO3],[KClO4,KCl])
    balanced_rx = Reaction(1.0,[KClO3],[KClO4,KCl], [4], [3,1])
    brxs = balance_reaction(rx)
    @test length(brxs) == 1
    @test isequal(balanced_rx, first(brxs))
end

# @reaction k, Al2(SO4)3 + 3Ca(OH)2 → 2Al(OH)3 + 3CaSO4
let
    @species Al(t) S(t) O(t) Ca(t) O(t) H(t)
    @compound Al2S3O12 ~ 2Al + 3S + 12O
    @compound CaO2H2 ~ Ca + 2O + 2H
    @compound AlO3H3 ~ Al + 3O + 3H
    @compound CaSO4 ~ Ca + S + 4O

    rx = Reaction(1.0,[Al2S3O12,CaO2H2],[AlO3H3,CaSO4])
    balanced_rx = Reaction(1.0,[Al2S3O12,CaO2H2],[AlO3H3,CaSO4], [1,3], [2,3])
    brxs = balance_reaction(rx)
    @test length(brxs) == 1
    @test isequal(balanced_rx, first(brxs))
end

# @reaction k, H2SO4 + 8HI → H2S + 4I2 + 4H2O
let
    @species H(t) S(t) O(t) I(t)
    @compound H2SO4 ~ 2H + S + 4O
    @compound HI ~ H + I
    @compound H2S ~ 2H + S
    @compound I2 ~ 2I
    @compound H2O ~ 2H + O

    rx = Reaction(1.0,[H2SO4,HI],[H2S,I2,H2O])
    balanced_rx = Reaction(1.0,[H2SO4,HI],[H2S,I2,H2O], [1,8], [1,4,4])
    brxs = balance_reaction(rx)
    @test length(brxs) == 1
    @test isequal(balanced_rx, first(brxs))
end

# @reaction k, C2H4 + 3O2 ↔ 2CO2 + 2H2O
let
    @species C(t) H(t) O(t)
    @compound C2H4 ~ 2C + 4H
    @compound O2 ~ 2O
    @compound CO2 ~ C + 2O
    @compound H2O ~ 2H + O

    rx = Reaction(1.0,[C2H4,O2],[CO2,H2O])
    balanced_rx = Reaction(1.0,[C2H4,O2],[CO2,H2O],[1,3],[2,2])
    brxs = balance_reaction(rx)
    @test length(brxs) == 1
    @test isequal(balanced_rx, first(brxs))
end

# Infinite solutions
let
    @species C(t) H(t) O(t)
    @compound CO ~ C + O
    @compound CO2 ~ C + 2O
    @compound H2 ~ 2H
    @compound CH4 ~ C + 4H
    @compound H2O ~ 2H + O

    rx = Reaction(1.0,[CO,CO2,H2],[CH4,H2O])
    brxs = balance_reaction(rx)
    @test_logs (:warn, r"Infinite balanced reactions*") match_mode=:any balance_reaction(rx)
    @test length(brxs) == 2
end

# No way to balance
let
    @species Fe(t) S(t) O(t) H(t) N(t)

    @compound FeS2 ~ Fe + 2S
    @compound HNO3 ~ H + N + 3O
    @compound Fe2S3O12 ~ 2Fe + 3S + 12O
    @compound NO ~ N + O
    @compound H2SO4 ~ 2H + S + 4O

    rx = Reaction(1.0,[FeS2,HNO3],[Fe2S3O12,NO,H2SO4])
    brxs = balance_reaction(rx)
    @test_logs (:warn, "Unable to balance reaction.") match_mode=:any balance_reaction(rx)
    @test isempty(brxs)
end

# test errors on compounds of compounds
let
    @species C(t) H(t) O(t)
    @compound CO ~ C + O
    @compound H2 ~ 2H
    @compound COH2 ~ CO + H2

    rx = Reaction(1.0, [CO, H2], [COH2])
    @test_throws Catalyst.COMPOUND_OF_COMPOUND_ERROR balance_reaction(rx)
end

# Checks that balancing work for a reaction from a reaction_network.
let
    rn = complete(@reaction_network begin
        @species C(t) H(t) O(t)
        @compounds begin
            O2 ~ 2O
            CO2 ~ 1C + 2O
            H2O ~ 2H + 1O
            C6H12O6 ~ 6C + 12H + 6O
        end
        k, CO2 + H2O --> C6H12O6 + O2
    end)
    
    brxs = balance_reaction(reactions(rn)[1])[1]
    
    @test isequal(rn.k, brxs.rate)
    @test isequal([rn.CO2,  rn.H2O], brxs.substrates)
    @test isequal([rn.C6H12O6,  rn.O2], brxs.products)
    @test isequal([6, 6], brxs.substoich)
    @test isequal([1,  6], brxs.prodstoich)    
end