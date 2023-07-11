### Fetch Packages and Set Global Variables ###
using Catalyst, ModelingToolkit, Test
@variables t

let 
    @parameters k
    @species H(t) O(t)
    @compound H2O(t) 1H 2O
    @test typeof(H2O) == Num
    arr_components = [1H ,2O]
    @test all((string(c1) == string(c2) for (c1, c2) = zip(components(H2O), arr_components)))
    @test iscompound(H2O) == true
end

let 
    @parameters k
    @species C(t) H(t) O(t)
    @compound C6H12O6(t) 6C 12H 6O
    @test typeof(C6H12O6) == Num
    arr_components = [6C ,12H, 6O]
    @test all((string(c1) == string(c2) for (c1, c2) = zip(components(C6H12O6), arr_components)))
    @test iscompound(C6H12O6) == true
end