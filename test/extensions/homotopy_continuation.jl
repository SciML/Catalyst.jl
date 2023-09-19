### Fetch Packages ###
using Catalyst, OrdinaryDiffEq, Test
import HomotopyContinuation

### Run Tests ###

# Tests for network without conservation laws.
# Tests for Symbol parameter input.
# Tests for Symbolics initial condiiton input.
# Tests for different types (Symbol/Symbolics) for parameters and initial conditions.
let 
    rs = @reaction_network begin
        (k1,k2), X1 <--> X2
        (k3,k4), 2X2 + X3 <--> X2_2X3
    end
    @unpack k1, k2, k3, k4 = rs
    ps = [k1 => 1.0, k2 => 2.0, k3 => 2.0, k4 => 2.0]
    u0 = [:X1 => 2.0, :X2 => 2.0, :X3 => 2.0, :X2_2X3 => 2.0]

    sim_ss = solve(ODEProblem(rs, u0, (0.0,1000.0), ps), Tsit5(); abstol=1e-12, reltol=1e-12)[end]
    hc_ss = hc_steady_states(rs, ps; u0=u0)[1]
    @test sim_ss ≈ hc_ss
end

# Tests for network with multiple steady state.
# Tests for Symbol parameter input.
# Tests tha passing kwargs to HC.solve does not error.
let 
    wilhelm_2009_model = @reaction_network begin
        k1, Y --> 2X
        k2, 2X --> X + Y
        k3, X + Y --> Y
        k4, X --> 0
    end
    ps = [:k3 => 1.0, :k2 => 2.0, :k4 => 1.5, :k1=>8.0]

    hc_ss_1 = hc_steady_states(wilhelm_2009_model, ps, seed=0x000004d1)
    @test sort(hc_ss_1, by=sol->sol[1]) ≈ [[0.0, 0.0], [0.5, 2.0], [4.5, 6.0]]

    hc_ss_2 = hc_steady_states(wilhelm_2009_model, ps, seed=0x000004d2)
    hc_ss_3 = hc_steady_states(wilhelm_2009_model, ps, seed=0x000004d2)
    @test hc_ss_1 != hc_ss_2
    @test hc_ss_2 == hc_ss_3
end

# Tests that reordering is correct.
# Tests corectness in presence of default values.
# Tests where some defaul values are overwritten with other values
let    
    rs_1 = @reaction_network begin
        @parameters kX1=1.0 kX2=2.0 kY1=12345.0 
        @species X1(t)=0.1 X2(t)=0.2 Y1(t)=12345.0
        (kX1,kX2), X1 <--> X2
        (kY1,kY2), Y1 <--> Y2
        (kZ1,kZ2), Z1 <--> Z2
    end
    ps = [:kY1 => 1.0, :kY2 => 3.0, :kZ1 => 1.0, :kZ2 => 4.0]
    u0_1 = [:Y1 => 1.0, :Y2 => 3.0, :Z1 => 10.0, :Z2 =>40.0]
    
    ss_1 = sort(hc_steady_states(rs_1, ps; u0=u0_1), by=sol->sol[1])
    @test ss_1 ≈ [[0.2, 0.1, 3.0, 1.0, 40.0, 10.0]]
    
    rs_2 = @reaction_network begin
        @parameters kX1=1.0 kX2=2.0 kY1=12345.0 
        @species C2(t)=0.1 C1(t)=0.2 B2(t)=12345.0
        (kX1,kX2), C2 <--> C1
        (kY1,kY2), B2 <--> B1
        (kZ1,kZ2), A2 <--> A1
    end
    u0_2 = [:B2 => 1.0, :B1 => 3.0, :A2 => 10.0, :A1 =>40.0]
    
    ss_2 = sort(hc_steady_states(rs_2, ps; u0=u0_2), by=sol->sol[1])
    @test ss_1 ≈ ss_2
end