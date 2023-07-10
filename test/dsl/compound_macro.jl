### Fetch Packages and Set Global Variables ###
using Catalyst, ModelingToolkit
@variables t

let 
    @parameters k
    @species H(t) O(t)
    @test typeof(@compound H2O(t) 1H 2O) == Vector{Num}
    @test iscompound(H2O) == true
end

let 
    @parameters k
    @species C(t) H(t) O(t)
    @test typeof(@compound C6H12O6(t) 6C 12H 6O) == Vector{Num}
    @test typeof(C6H12O6) == Num
    @test iscompound(C) == false
    @test iscompound(C6H12O6) == true
end