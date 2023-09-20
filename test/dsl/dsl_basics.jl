#! format: off

### Fetch Packages and Set Global Variables ###

using Catalyst
@variables t

### Naming Tests ###

let
    @parameters k
    @species A(t)
    rx = Reaction(k, [A], nothing)
    function rntest(rn, name)
        @test nameof(rn) == name
        @test isequal(species(rn)[1], ModelingToolkit.unwrap(A))
        @test isequal(parameters(rn)[1], ModelingToolkit.unwrap(k))
        @test reactions(rn)[1] == rx
    end

    function emptyrntest(rn, name)
        @test nameof(rn) == name
        @test numreactions(rn) == 0
        @test numspecies(rn) == 0
        @test numreactionparams(rn) == 0
    end

    rn = @reaction_network name begin
        @parameters k
        k, A --> 0
    end
    rntest(rn, :name)

    name = :blah
    rn = @reaction_network $name begin
        @parameters k
        k, A --> 0
    end
    rntest(rn, :blah)

    rn = @reaction_network begin
        @parameters k
        k, A --> 0
    end
    rntest(rn, nameof(rn))

    function makern(; name)
        @reaction_network $name begin
        @parameters k
            k, A --> 0
        end
    end
    @named testnet = makern()
    rntest(testnet, :testnet)

    rn = @reaction_network name
    emptyrntest(rn, :name)

    rn = @reaction_network $name
    emptyrntest(rn, :blah)
end


### Test Interpolation Within the DSL ###

@parameters α k k1 k2
@species A(t) B(t) C(t) D(t)

let
    AA = A
    AAA = A^2 + B
    rn = @reaction_network rn begin
        @parameters k
        @species A(t) B(t) C(t) D(t)
        k*$AAA, C --> D
    end
    rn2 = ReactionSystem([Reaction(k*AAA, [C], [D])], t; name=:rn)
    @test rn == rn2

    rn = @reaction_network rn begin
        @parameters k
        @species A(t) C(t) D(t)
        k, $AA + C --> D
    end
    rn2 = ReactionSystem([Reaction(k, [AA,C], [D])], t; name=:rn)
    @test rn == rn2
end

let
    BB = B; A2 = A
    rn = @reaction_network rn begin
        @parameters k1 k2
        (k1,k2), C + $A2 + $BB + $A2 <--> $BB + $BB
    end
    rn2 = ReactionSystem([Reaction(k1, [C, A, B], [B], [1,2,1],[2]),
                        Reaction(k2, [B], [C, A, B], [2], [1,2,1])],
                        t; name=:rn)
    @test rn == rn2
end

let
    AA = A
    kk1 = k^2*A
    kk2 = k1+k2
    rn = @reaction_network rn begin
        @parameters α k k1 k2
        α+$kk1*$kk2*$AA, 2*$AA + B --> $AA
    end
    rn2 = ReactionSystem([Reaction(α+kk1*kk2*AA, [A, B], [A], [2, 1], [1])], t; name=:rn)
    @test rn == rn2
end

@testset "make_reaction_system can be called from another module" begin
    ex = quote
        (Ka, Depot --> Central)
        (CL / Vc, Central --> 0)
    end
    # Line number nodes aren't ignored so have to be manually removed
    Base.remove_linenums!(ex)
    @test eval(Catalyst.make_reaction_system(ex)) isa ReactionSystem
end

let
    rx = @reaction k*h, A + 2*B --> 3*C + D
    @parameters k h
    @species A(t) B(t) C(t) D(t)
    @test rx == Reaction(k*h,[A,B],[C,D],[1,2],[3,1])

    ex = k*A^2 + B
    V  = A
    rx = @reaction b+$ex, 2*$V + C--> ∅
    @parameters b
    @test rx == Reaction(b+ex, [A,C], nothing, [2,1], nothing)
end


### Other tests ###

# Test floating point stoichiometry work.
let
    @parameters k
    @species B(t) C(t) D(t)
    rx1 = Reaction(k,[B,C],[B,D], [2.5,1],[3.5, 2.5])
    rx2 = Reaction(2*k, [B], [D], [1], [2.5])
    rx3 = Reaction(2*k, [B], [D], [2.5], [2])
    @named mixedsys = ReactionSystem([rx1,rx2,rx3],t,[B,C,D],[k])
    osys = convert(ODESystem, mixedsys; combinatoric_ratelaws=false)
    rn = @reaction_network mixedsys begin
        @parameters k
        k, 2.5*B + C --> 3.5*B + 2.5*D
        2*k, B --> 2.5*D
        2*k, 2.5*B --> 2*D
    end
    @test rn == mixedsys
end

# Test variables that appear only in rates and aren't ps
# are categorized as species.
let
    rn = @reaction_network begin
        @parameters k k2 n
        @species A(t) B(t) C(t) D(t) H(t)
        π*k*D*hill(B,k2,B*D*H,n), 3*A  --> 2*C
    end
    @parameters k k2 n
    @variables t
    @species A(t) B(t) C(t) D(t) H(t)
    @test issetequal([A,B,C,D,H], species(rn))
    @test issetequal([k,k2,n], parameters(rn))
end

# Test species have the right metadata via the DSL.
let
    rn = @reaction_network begin
        k, 2*A + B --> C
    end
    @test issetequal(states(rn), species(rn))
    @test all(isspecies, species(rn))

    rn2 = @reaction_network begin
        @species A(t) = 1 B(t) = 2 [isbcspecies = true]
        k, A + 2*B --> 2*B
    end
    @variables t
    @unpack A,B = rn2
    D = Differential(t)
    eq = D(B) ~ -B
    @named osys = ODESystem([eq], t)
    @named rn2 = extend(osys, rn2)
    @test issetequal(states(rn2), species(rn2))
    @test all(isspecies, species(rn))
    @test Catalyst.isbc(ModelingToolkit.value(B))
    @test Catalyst.isbc(ModelingToolkit.value(A)) == false
    osys2 = convert(ODESystem, rn2)
    @test issetequal(states(osys2), states(rn2))
    @test length(equations(osys2)) == 2
end

# test @variables in DSL
let
    rn = @reaction_network tester begin
        @parameters k1
        @variables V1(t) V2(t) V3(t)
        @species B1(t) B2(t)
        (k1*k2 + V3), V1*A + 2*B1 --> V2*C + B2
    end

    @parameters k1 k2
    @variables t V1(t) V2(t) V3(t)
    @species A(t) B1(t) B2(t) C(t)
    rx = Reaction(k1*k2 + V3, [A, B1], [C, B2], [V1, 2], [V2, 1])
    @named tester = ReactionSystem([rx], t)
    @test tester == rn

    sts = (A, B1, B2, C, V1, V2, V3)
    spcs = (A, B1, B2, C)
    @test issetequal(states(rn), sts)
    @test issetequal(species(rn), spcs)

    @test_throws ArgumentError begin
        rn = @reaction_network begin
            @variables K
            k, K*A --> B
        end
    end
end

# ivs test
let
    rn = @reaction_network ivstest begin
        @ivs s x
        @parameters k2
        @variables D(x) E(s) F(s,x)
        @species A(s,x) B(s) C(x)
        k*k2*D, E*A +B --> F*C + C2
    end
    @parameters k k2
    @variables s x D(x) E(s) F(s,x)
    @species A(s,x) B(s) C(x) C2(s,x)
    rx = Reaction(k*k2*D, [A, B], [C, C2], [E, 1], [F, 1])
    @named ivstest = ReactionSystem([rx], s; spatial_ivs = [x])
    @test ivstest == rn
    @test issetequal(states(rn), [D, E, F, A, B, C, C2])
    @test issetequal(species(rn), [A, B, C, C2])
    @test isequal(ModelingToolkit.get_iv(rn), s)
    @test issetequal(Catalyst.get_sivs(rn), [x])
end

# array variables test
let
    rn = @reaction_network arrtest begin
        @parameters k[1:2] a
        @variables (V(t))[1:2] W(t)
        @species (X(t))[1:2] Y(t)
        k[1]*a+k[2], X[1] + V[1]*X[2] --> V[2]*W*Y + B*C
    end

    @parameters k[1:3] a B
    @variables t (V(t))[1:2] W(t)
    @species (X(t))[1:2] Y(t) C(t)
    rx = Reaction(k[1]*a+k[2], [X[1], X[2]], [Y, C], [1, V[1]], [V[2] * W, B])
    @named arrtest = ReactionSystem([rx], t)
    arrtest == rn
end
