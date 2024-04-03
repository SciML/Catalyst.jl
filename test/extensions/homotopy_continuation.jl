### Prepares Tests ###

# Fetch packages.
using Catalyst, Test
import HomotopyContinuation

# Fetch test functions.
include("../test_functions.jl")

### Basic Tests ###

# Tests for network without conservation laws.
# Tests for Symbol parameter input.
# Tests for Symbolics initial condition input.
# Tests for different types (Symbol/Symbolics) for parameters and initial conditions.
# Tests that attempts to find steady states of system with conservation laws, while u0 is not provided, gives an error.
let 
    # Creates the model.
    rs = @reaction_network begin
        (k1,k2), X1 <--> X2
        (k3,k4), 2X2 + X3 <--> X2_2X3
    end
    @unpack k1, k2, k3, k4 = rs
    ps = [k1 => 1.0, k2 => 2.0, k3 => 2.0, k4 => 2.0]
    u0 = [:X1 => 2.0, :X2 => 2.0, :X3 => 2.0, :X2_2X3 => 2.0]

    # Computes the single steady state, checks that when given to the ODE rhs, all are evaluated to 0.
    hc_ss = hc_steady_states(rs, ps; u0=u0, show_progress=false)
    hc_ss = Pair.(unknowns(rs), hc_ss[1])
    @test maximum(abs.(f_eval(rs, hc_ss, ps, 0.0))) ≈ 0.0 atol=1e-12

    # Checks that not giving a `u0` argument yields an error for systems with conservation laws.
    @test_throws Exception hc_steady_states(rs, ps; show_progress=false)
end

# Tests for network with multiple steady state.
# Tests for Symbol parameter input.
# Tests that passing kwargs to HC.solve does not error.
let 
    wilhelm_2009_model = @reaction_network begin
        k1, Y --> 2X
        k2, 2X --> X + Y
        k3, X + Y --> Y
        k4, X --> 0
    end
    ps = [:k3 => 1.0, :k2 => 2.0, :k4 => 1.5, :k1=>8.0]

    hc_ss_1 = hc_steady_states(wilhelm_2009_model, ps; seed=0x000004d1, show_progress=false)
    @test sort(hc_ss_1, by=sol->sol[1]) ≈ [[0.0, 0.0], [0.5, 2.0], [4.5, 6.0]]

    hc_ss_2 = hc_steady_states(wilhelm_2009_model, ps; seed=0x000004d2, show_progress=false)
    hc_ss_3 = hc_steady_states(wilhelm_2009_model, ps; seed=0x000004d2, show_progress=false)
    @test hc_ss_1 != hc_ss_2
    @test hc_ss_2 == hc_ss_3
end

# Tests that reordering is correct.
# Tests correctness in presence of default values.
# Tests where some default values are overwritten with other values.
# Tests where input ps/u0 are tuples with mixed types.
let
    rs_1 = @reaction_network begin
        @parameters kX1=1.0 kX2=2.0 kY1=12345.0 
        @species X1(t)=0.1 X2(t)=0.2 Y1(t)=12345.0
        (kX1,kX2), X1 <--> X2
        (kY1,kY2), Y1 <--> Y2
        (kZ1,kZ2), Z1 <--> Z2
    end
    ps = (:kY1 => 1.0, :kY2 => 3, :kZ1 => 1.0, :kZ2 => 4.0)
    u0_1 = (:Y1 => 1.0, :Y2 => 3, :Z1 => 10, :Z2 =>40.0)
    
    ss_1 = sort(hc_steady_states(rs_1, ps; u0=u0_1, show_progress=false), by=sol->sol[1])
    @test ss_1 ≈ [[0.2, 0.1, 3.0, 1.0, 40.0, 10.0]]
    
    rs_2 = @reaction_network begin
        @parameters kX1=1.0 kX2=2.0 kY1=12345.0 
        @species C2(t)=0.1 C1(t)=0.2 B2(t)=12345.0
        (kX1,kX2), C2 <--> C1
        (kY1,kY2), B2 <--> B1
        (kZ1,kZ2), A2 <--> A1
    end
    u0_2 = [:B2 => 1.0, :B1 => 3.0, :A2 => 10.0, :A1 =>40.0]
    
    ss_2 = sort(hc_steady_states(rs_2, ps; u0=u0_2, show_progress=false), by=sol->sol[1])
    @test ss_1 ≈ ss_2
end

# Tests that non-scalar reaction rates work.
# Tests that rational polynomial steady state systems work. 
# Tests that Hill function is correctly expanded even if nested.
# Test filter_negative=false works.
# Tests than non-integer exponents throws an error.
let 
    rs = @reaction_network begin
        v*(0.1/v + hill(X,1,K,n)), 0 --> X
        d, X --> 0
    end
    ps = [:v => 5.0, :K => 2.5, :n => 3, :d => 1.0]
    sss = hc_steady_states(rs, ps; filter_negative=false, show_progress=false)
    
    @test length(sss) == 4
    for ss in sss
        @test f_eval(rs,sss[1], last.(ps), 0.0)[1] ≈ 0.0 atol=1e-12
    end

    @test_throws Exception hc_steady_states(rs, [:v => 5.0, :K => 2.5, :n => 2.7, :d => 1.0]; show_progress=false)
end


### Other Tests ###

# Checks that `hc_steady_states` cannot be applied to non-complete `ReactionSystems`s.
let 
    # Create model.
    incomplete_network = @network_component begin
        (p, d), 0 <--> X
    end
    p_start = [:p => 1.0, :d => 0.2]
    
    # Computes bifurcation diagram.
    @test_throws Exception hc_steady_states(incomplete_network, p_start)
end