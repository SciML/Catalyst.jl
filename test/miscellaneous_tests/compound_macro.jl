### Fetch Packages and Set Global Variables ###
using Catalyst

using ModelingToolkit, Test
const MT = ModelingToolkit

let
    @variables t
    @species C(t) H(t) O(t)
    @compound C6H12O2(t) 6C 12H 2O 
    #WIP
    # @test isequal([C, H, O], components(C6H12O2))
    # @test isequal([6, 12, 2], coefficients(C6H12O2))
end
