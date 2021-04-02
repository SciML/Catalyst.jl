using MacroTools, Test
include("../src/Catalyst.jl")

systems = [ 
    # graphs.jl
    quote
        α, S + I --> 2I
        β, I --> R
        S^2, R --> 0
    end => (:α, :β),
            
    # custom_functions.jl
    quote
        hill(X,v1,K1,2), X + Y --> Z1
        mm(X,v2,K2), X + Y --> Z2
        p1*X^2+p2, X + Y --> Z3
        exp(-p3*Y), X + Y --> Z4
        hillR(X,v3,K3,2), X + Y --> Z5
        mmR(X,v4,K4), X + Y --> Z6
    end => (:v1, :K1, :v2, :K2, :p1, :p2, :p3, :v3, :K3, :v4, :K4),
         
    # higher_order_reactions.jl
    quote
        p, ∅ ⟼ X1
        r1, 2X1 ⟼ 3X2
        mm(X1,r2,K), 3X2 ⟼ X3 + 2X4
        r3, X3 + 2X4 ⟼ 3X5 + 3X6
        r4*X2, 3X5 + 3X6 ⟼ 3X5 + 2X7 + 4X8
        r5, 3X5 + 2X7 + 4X8 ⟼ 10X9
        r6, 10X9 ⟼ X10
        d, 2X10 ⟼ ∅
    end => (:p, :r1, :r2, :K, :r3, :r4, :r5, :r6, :d),
                
    quote
        p, ∅ ⟾ X1
        r1*X1^2/factorial(2), 2X1 ⟾ 3X2
        mm(X1,r2,K)*X2^3/factorial(3), 3X2 ⟾ X3 + 2X4
        r3*X3*X4^2/factorial(2), X3 + 2X4 ⟾ 3X5 + 3X6
        r4*X2*X5^3*X6^3/(factorial(3)*factorial(3)), 3X5 + 3X6 ⟾ 3X5 + 2X7 + 4X8
        r5*X5^3*X7^2*X8^4/(factorial(3)*factorial(2)*factorial(4)), 3X5 + 2X7 + 4X8 ⟾ 10X9
        r6*X9^10/factorial(10), 10X9 ⟾ X10
        d*X10^2/factorial(2), 2X10 ⟾ ∅
    end => (:p, :r1, :r2, :K, :r3, :r4, :r5, :r6, :d),
    
    # latexify.jl
    quote
        hillR(X2,v1,K1,n1)*hill(X4,v1,K1,n1), ∅ → X1
        hill(X5,v2,K2,n2), ∅ → X2
        hill(X3,v3,K3,n3), ∅ → X3
        hillR(X1,v4,K4,n4), ∅ → X4
        hill(X2,v5,K5,n5), ∅ → X5
        (k1,k2), X2 ⟷ X1 + 2X4
        (k3,k4), X4 ⟷ X3
        (k5,k6), 3X5 + X1 ⟷ X2
        (d1,d2,d3,d4,d5), (X1,X2,X3,X4,X5)  ⟶ ∅
    end => (:v1, :K1, :n1, :v2, :K2, :n2, :v3, :K3, :n3, :v4, :K4, :n4, :v5, 
            :K5, :n5, :k1, :k2, :k3, :k4, :k5, :k6, :d1, :d2, :d3, :d4, :d5),

    quote
        (hill(B, p_a, k, n), d_a), 0 ↔ A
        (p_b, d_b), 0 ↔ B
        (r_a, r_b), 3B ↔ A
    end => (:p_a, :k, :n, :d_a, :p_b, :d_b, :r_a, :r_b),

    # make_jacobian.jl
    quote
        (2.0,1.0),   ∅ ↔ X
        (3.0,1.0),   ∅ ↔ Y
        (5.0,2.0),   X + Y ↔ XY
    end => (),

    quote
        (p1,1.0),   ∅ ↔ X
        (p2,1.0),   ∅ ↔ Y
        (p3*X,1.0), X + Y ↔ XY
    end => (:p1, :p2, :p3),

    quote
        k1, 2A → B
        k2, B → 2A
        k3, A + B → C
        k4, C → A + B
        k5, 3C → 3A
        k6, 0 → 2B
        hill(A,k7,k8,2), ∅ → B
    end => (:k1, :k2, :k3, :k4, :k5, :k6, :k7, :k8),
    
    # test_networks.jl              
    quote
        (p1,p2,p3), ∅ → (X1,X2,X3)
        (k1,k2), X2 ⟷ X1 + 2X3
        (k3,k4), X1 ⟷ X3
        (d1,d2,d3), (X1,X2,X3) → ∅
    end => (:p1, :p2, :p3, :k1, :k2, :k3, :k4, :d1, :d2, :d3),

    quote
        mmR(X2,v1,K1), ∅ → X1
        mm(X1,v2,K2), ∅ → X2
        d, X1+X2 → ∅
    end => (:v1, :K1, :v2, :K2, :d),

    quote
        mm(X2,v1,K1), ∅ → X1
        mm(X3,v2,K2), ∅ → X2
        (k1,k2), X1 ⟷ X3
        (k3,k4), X3 + X2 ⟷ X4 + X1
        d, (X1,X2,X3,X4) → ∅
    end => (:v1, :K1, :v2, :K2, :k1, :k2, :k3, :k4, :d),

    quote
        mmR(X4,v1,K1), ∅ → X1
        mmR(X1,v2,K2), ∅ → X2
        mmR(X2,v3,K3), ∅ → X3
        mmR(X3,v4,K4), ∅ → X4
        (d1,d2,d3,d4), (X1,X2,X3,X4) → ∅
    end => (:v1, :K1, :v2, :K2, :v3, :K3, :v4, :K4, :d1, :d2, :d3, :d4),

    quote
        p, ∅ → X1
        (k1,k2), X1 ⟷ X2
        (k3,k4), X2 ⟷ X3
        (k5,k6), X3 ⟷ X4
        d, X4 → ∅
    end => (:p, :k1, :k2, :k3, :k4, :k5, :k6, :d),

    quote
        (p1,p2), ∅ → (X1,X2)
        (k1,k2), 2X1 ⟷ X3
        (k3,k4), X2 + X3 ⟷ 4X4
        (k5,k6), X4 + X1 ⟷ 2X3
        d, (X1,X2,X3,X4) → ∅
    end => (:p1, :p2, :k1, :k2, :k3, :k4, :k5, :k6, :d),

    ### Reaction networks that are actual models that have been used ###

    # Brusselator.
    quote
        A, ∅ → X
        1, 2X + Y → 3X
        B, X → Y
        1, X → ∅
    end => (:A, :B),

    # The B.subtilis σV Lysozyme stress response.
    quote
        v0 + hill(σ,v,K,n), ∅ → (σ+A)
        deg, (σ,A,Aσ) → ∅
        (kB,kD), A + σ ↔ Aσ
        S*kC, Aσ → σ
    end => (:v0, :v, :K, :n, :kD, :kB, :kC, :deg, :S),

    # A cell cycle model
    quote
      k1, 0 --> Y
      k2p, Y --> 0
      k2pp*P, Y --> 0
      (k3p+k3pp*A)/(J3+Po), Po-->P
      (k4*m)/(J4+P), Y + P --> Y + Po
    end => (:k1, :k2p, :k2pp, :k3p, :k3pp, :A, :J3, :k4, :m, :J4),

    # A bistable switch
    quote
        d,    (X,Y) → ∅
        hillR(Y,v1,K1,n1), ∅ → X
        hillR(X,v2,K2,n2), ∅ → Y
    end => (:d, :v1, :K1, :n1, :v2, :K2, :n2),

    ### Reaction networks that contain weird functions, stuff, and other oddities ###

    quote
        exp(p), ∅ → X
        d*exp(X), X → ∅
    end => (:p, :d),

    quote
        k1, ∅ → X
        k2*log(12+X), X → Y
        k3*log(3+Y), Y → Z
        log(5,6+k4), Z → ∅
    end => (:k1, :k2, :k3, :k4),

    quote
        d,    (X,Y) → ∅
        hillR(hill(Y,v2,K2,n2),v1,K1,n1), ∅ → X
        hillR(hill(X,v1,K1,n1),v2,K2,n2), ∅ → Y
    end => (:d, :v1, :K1, :n1, :v2, :K2, :n2),

    quote
        p, ∅ → X1
        (k1,X2^X3), X1 ⟷ X2
        (X2/X4,k4), X2 ⟷ X3
        (k5,X1*X2*X4), X3 ⟷ X4
        d*X2, X4 → ∅
    end => (:p, :k1, :k2, :k3, :k4, :k5, :k6, :d),

    quote
        (k1*tanh(X1),k2), X1 ↔ 2X1
        (2+sin(X1),3+sin(k3)), X1 + X2 ↔ X3
        (k5,4+cos(X3+X1)+sin(X2)), X3 ↔ X2
    end => (:k1, :k2, :k3, :k4, :k5, :k6),

    quote
        (p,d), ∅ ↔ X1
        (p,d), ∅ ↔ X1
        (p,d), ∅ ↔ X1
        (p,d), ∅ ↔ X1
        (p,d), ∅ ↔ X1
        (p,d), ∅ ↔ X1
    end => (:p, :d),

    quote
        k1+X3, X1 → X2
        k2*X4, X2 → X3
        k3/(X5+0.01), X3 → X4
        X6/k4, X4 → X5
        X7^k5, X5 → X6
        k6^X8, X6 → X7
        sqrt(abs(k7*X9)), X7 → X8
        cbrt(abs(k8+X1)), X8 → X9
        X2^3+2X2^2+k9, X9 → X1
    end => (:k1, :k2, :k3, :k4, :k5, :k6, :k7, :k8, :k9),

    quote
        p, ∅ → X1
        hill(hill(hill(hill(X1,v1,K1,n1),v2,K2,n2),v3,K3,n3),v4,K4,n4), X1 → ∅
    end => (:p, :v1, :K1, :n1, :v2, :K2, :n2, :v3, :K3, :n3, :v4, :K4, :n4),
                                                                                                            
    # solve_ODEs.jl
    quote
        (k1,k2), X1 ↔ X2
        (k3,k4), X3+ X4 ↔ X5
        (k5,k6), 2X6 ↔ 3X7
        (k7,k8), ∅ ↔ X8
    end => (:k1, :k2, :k3, :k4, :k5, :k6, :k7, :k8),


    quote
        (1.2,5), X1 ↔ X2
    end => (),
]
            
for (expr, params) in systems
    expr = MacroTools.striplines(expr)
    expr_spatial = Catalyst.make_spatial_reaction_network(expr, [])
            
    rn = Catalyst.make_reaction_system(expr, params)
    rn_copy = Catalyst.make_reaction_system(expr_spatial, params)

    @test rn == rn_copy
end
