#! format: off

### Fetch Packages and Set Global Variables ###

# Fetch packages.
using Catalyst, JumpProcesses, OrdinaryDiffEqTsit5, StochasticDiffEq, Test
using Catalyst: isequivalent
using ModelingToolkitBase: get_brownians, get_jumps, Pre
import ModelingToolkitBase as MT

# Sets stable rng number.
using StableRNGs
rng = StableRNG(12345)

# Sets the default `t` to use as time variable.
t = default_t()
D = default_time_deriv()

### Test Brownians and Jumps in ReactionSystem ###

# Tests creating a ReactionSystem with user-provided brownians via the 5-arg constructor.
let
    @brownians B1 B2
    @variables V(t)
    @species S(t) P(t)
    @parameters k

    # Create a ReactionSystem with brownians as 5th positional argument.
    @named rn = ReactionSystem(
        [Reaction(k, [S], [P]), D(V) ~ -V + B1 + B2],
        t, [S, P, V], [k], [B1, B2]
    )

    # Test that brownians are stored correctly.
    @test length(get_brownians(rn)) == 2
    @test Set(get_brownians(rn)) == Set([B1, B2])

    # Test recursive accessor (no subsystems, so same as get_brownians).
    @test length(MT.brownians(rn)) == 2
end

# Tests creating a ReactionSystem with user-provided jumps.
let
    @species S(t)
    @parameters k d

    # Create a ConstantRateJump manually.
    rate = k
    affect = [S ~ Pre(S) + 1]
    user_jump = MT.ConstantRateJump(rate, affect)

    # Create a ReactionSystem with jumps as keyword argument.
    @named rn = ReactionSystem(
        [Reaction(d, [S], nothing)],
        t, [S], [k, d];
        jumps = [user_jump]
    )

    # Test that jumps are stored correctly.
    @test length(get_jumps(rn)) == 1
    @test get_jumps(rn)[1] === user_jump

    # Test recursive accessor.
    @test length(MT.jumps(rn)) == 1
end

# Tests auto-discovery of brownians via the two-argument constructor when equations are provided.
# Note: Auto-discovery works through MT.System which extracts brownians from equations.
let
    @brownians B
    @variables V(t)
    @species S(t) P(t)
    @parameters k

    # Two-argument constructor with equations containing a brownian.
    # The brownian should be auto-discovered via MT.System.
    rn = ReactionSystem([Reaction(k, [S], [P]), D(V) ~ -V + B], t; name = :rn)

    # Brownians should be auto-discovered from equations.
    @test length(get_brownians(rn)) == 1

    # V should be an unknown, not a brownian.
    @test any(isequal(V), unknowns(rn))
    @test !any(isequal(V), get_brownians(rn))
end

# Tests auto-discovery of variables from ConstantRateJump via the two-argument constructor.
let
    @species S(t) P(t)
    @parameters k d

    # Create a ConstantRateJump that references S and k.
    user_jump = MT.ConstantRateJump(k, [S ~ Pre(S) + 1])

    # Two-argument constructor - should auto-discover S, k from jump.
    # Only explicitly provide the reaction's variables (P, d).
    rn = ReactionSystem(
        [Reaction(d, [P], nothing)],
        t;
        jumps = [user_jump],
        name = :rn
    )

    # S from jump affect, P from reaction.
    @test issetequal(unknowns(rn), [S, P])
    # k from jump rate, d from reaction.
    @test issetequal(parameters(rn), [k, d])
end

# Tests auto-discovery of variables from VariableRateJump.
let
    @species X(t)
    @variables Y(t)
    @parameters α β

    # VariableRateJump with rate depending on Y and α.
    user_jump = MT.VariableRateJump(α * Y, [X ~ Pre(X) + 1])

    rn = ReactionSystem(
        [D(Y) ~ -β * Y],  # ODE for Y
        t;
        jumps = [user_jump],
        name = :rn
    )

    # X from jump affect, Y from equation and jump rate.
    @test issetequal(unknowns(rn), [X, Y])
    # α from jump rate, β from equation.
    @test issetequal(parameters(rn), [α, β])
end

# Tests auto-discovery of variables from MassActionJump.
let
    @species A(t) B(t)
    @parameters k_ma

    # MassActionJump: A -> B with rate k_ma.
    user_jump = MT.MassActionJump(k_ma, [A => 1], [A => -1, B => 1])

    rn = ReactionSystem([], t; jumps = [user_jump], name = :rn)

    # A and B from jump stoichiometry.
    @test issetequal(unknowns(rn), [A, B])
    # k_ma from jump rate.
    @test issetequal(parameters(rn), [k_ma])
end

### Test Composition with Brownians and Jumps ###

# Tests that brownians are collected from subsystems during flatten.
let
    @brownians B1 B2
    @variables V1(t) V2(t)
    @species S(t)
    @parameters k

    # Create two subsystems, each with a brownian.
    @named sub1 = ReactionSystem(
        [D(V1) ~ -V1 + B1],
        t, [V1], [], [B1]
    )
    @named sub2 = ReactionSystem(
        [D(V2) ~ -V2 + B2],
        t, [V2], [], [B2]
    )

    # Compose them under a parent system.
    @named parent = ReactionSystem(
        [Reaction(k, nothing, [S])],
        t, [S], [k];
        systems = [sub1, sub2]
    )

    # Before flattening, recursive accessor should collect all brownians.
    @test length(MT.brownians(parent)) == 2

    # After flattening, brownians should be merged.
    flat = Catalyst.flatten(parent)
    @test length(get_brownians(flat)) == 2
end

# Tests that jumps are collected from subsystems during flatten.
let
    @species S(t) P(t)
    @parameters k1 k2 d

    # Create a user jump.
    user_jump = MT.ConstantRateJump(k2, [P ~ Pre(P) + 1])

    # Create a subsystem with a user jump.
    @named sub = ReactionSystem(
        [Reaction(d, [P], nothing)],
        t, [P], [k2, d];
        jumps = [user_jump]
    )

    # Create a parent system.
    @named parent = ReactionSystem(
        [Reaction(k1, nothing, [S])],
        t, [S], [k1];
        systems = [sub]
    )

    # Recursive accessor should find the jump in the subsystem.
    @test length(MT.jumps(parent)) == 1

    # After flattening, jump should be present.
    flat = Catalyst.flatten(parent)
    @test length(get_jumps(flat)) == 1
end

# Tests extend with brownians and jumps (union semantics).
let
    @brownians B1 B2
    @variables V1(t) V2(t)
    @species S(t) P(t)
    @parameters k1 k2

    user_jump1 = MT.ConstantRateJump(k1, [S ~ Pre(S) + 1])
    user_jump2 = MT.ConstantRateJump(k2, [P ~ Pre(P) + 1])

    @named sys1 = ReactionSystem(
        [D(V1) ~ -V1 + B1],
        t, [V1, S], [k1], [B1];
        jumps = [user_jump1]
    )
    @named sys2 = ReactionSystem(
        [D(V2) ~ -V2 + B2],
        t, [V2, P], [k2], [B2];
        jumps = [user_jump2]
    )

    # Extend sys1 with sys2.
    extended = MT.extend(sys2, sys1; name = :extended)

    # Both brownians should be present.
    @test length(get_brownians(extended)) == 2

    # Both jumps should be present.
    @test length(get_jumps(extended)) == 2
end

# Tests compose (subsystem brownians/jumps accessible after flatten).
let
    @brownians B
    @variables V(t)
    @species S(t) P(t)
    @parameters k1 k2

    user_jump = MT.ConstantRateJump(k2, [P ~ Pre(P) + 1])

    @named sub = ReactionSystem(
        [D(V) ~ -V + B],
        t, [V], [k2], [B];
        jumps = [user_jump]
    )

    @named main = ReactionSystem(
        [Reaction(k1, [S], [P])],
        t, [S, P], [k1]
    )

    # Compose main with sub.
    composed = MT.compose(main, [sub]; name = :composed)

    # Recursive accessors should find brownians and jumps in subsystem.
    @test length(MT.brownians(composed)) == 1
    @test length(MT.jumps(composed)) == 1

    # After flatten, they should be in the top-level system.
    flat = Catalyst.flatten(composed)
    @test length(get_brownians(flat)) == 1
    @test length(get_jumps(flat)) == 1
end

### Test Conversion Error Checks ###

# Tests that ode_model errors on systems with brownians.
let
    @brownians B
    @variables V(t)
    @parameters k

    @named rn = ReactionSystem(
        [D(V) ~ -V + B],
        t, [V], [k], [B]
    )
    rn = complete(rn)

    @test_throws ErrorException Catalyst.ode_model(rn)
end

# Tests that ode_model errors on systems with user jumps.
let
    @species S(t)
    @parameters k

    user_jump = MT.ConstantRateJump(k, [S ~ Pre(S) + 1])

    @named rn = ReactionSystem(
        Reaction[], t, [S], [k];
        jumps = [user_jump]
    )
    rn = complete(rn)

    @test_throws ErrorException Catalyst.ode_model(rn)
end

# Tests that sde_model errors on systems with user jumps.
let
    @species S(t)
    @parameters k d

    user_jump = MT.ConstantRateJump(k, [S ~ Pre(S) + 1])

    @named rn = ReactionSystem(
        [Reaction(d, [S], nothing)],
        t, [S], [k, d];
        jumps = [user_jump]
    )
    rn = complete(rn)

    @test_throws ErrorException Catalyst.sde_model(rn)
end

# Tests that jump_model errors on systems with non-reaction equations.
let
    @variables V(t)
    @species S(t)
    @parameters k

    @named rn = ReactionSystem(
        [Reaction(k, nothing, [S]), D(V) ~ -V],
        t, [S, V], [k]
    )
    rn = complete(rn)

    @test_throws ErrorException Catalyst.jump_model(rn)
end

### Test Conversion Success with Merged Brownians/Jumps ###

# Tests that HybridProblem correctly handles user brownians (SDE path).
let
    @brownians B
    @variables V(t)
    @species S(t) P(t)
    @parameters k σ

    # System with ODE reaction and user brownian in a variable equation.
    @named rn = ReactionSystem(
        [Reaction(k, [S], [P]), D(V) ~ -V + σ * B],
        t, [S, P, V], [k, σ], [B]
    )
    rn = complete(rn)

    u0 = [:S => 10.0, :P => 0.0, :V => 1.0]
    ps = [:k => 1.0, :σ => 0.1]
    tspan = (0.0, 1.0)

    # Should create an SDEProblem (has brownians → SDE path).
    # Note: default_scale=ODE ensures the reaction is treated as ODE, not Jump.
    prob = HybridProblem(rn, u0, tspan, ps; default_scale = PhysicalScale.ODE, structural_simplify = true)
    @test prob isa SDEProblem

    # Should be solvable.
    sol = solve(prob, EM(); dt = 0.01)
    @test sol.retcode == ReturnCode.Success
end

# Tests that HybridProblem correctly handles user jumps (Jump path).
let
    @species S(t) P(t)
    @parameters k1 k2

    user_jump = MT.ConstantRateJump(k2, [P ~ Pre(P) + 1])

    # System with reaction and user jump (both Jump scale by default).
    @named rn = ReactionSystem(
        [Reaction(k1, [S], [P])],
        t, [S, P], [k1, k2];
        jumps = [user_jump]
    )
    rn = complete(rn)

    u0 = [:S => 10, :P => 0]  # Integer u0 for pure jump
    ps = [:k1 => 1.0, :k2 => 0.5]
    tspan = (0.0, 10.0)

    # Should create a JumpProblem (has user jumps → Jump path).
    prob = HybridProblem(rn, u0, tspan, ps)
    @test prob isa JumpProblem

    # Should be solvable with SSAStepper (pure jump).
    sol = solve(prob, SSAStepper())
    @test sol.retcode == ReturnCode.Success
end

# Tests that hybrid_model merges user brownians with reaction brownians.
let
    @brownians B_user
    @variables V(t)
    @species S(t) P(t)
    @parameters k σ

    # System with SDE-scale reaction (generates reaction brownian) and user brownian.
    @named rn = ReactionSystem(
        [Reaction(k, [S], [P]), D(V) ~ -V + σ * B_user],
        t, [S, P, V], [k, σ], [B_user]
    )
    rn = complete(rn)

    # Use hybrid_model with SDE scale for the reaction (index 1).
    sys = Catalyst.hybrid_model(rn;
        physical_scales = [1 => PhysicalScale.SDE],
        name = :test_sys
    )

    # Should have both reaction-generated and user brownians.
    # Reaction brownian + user brownian = at least 2.
    @test length(get_brownians(sys)) == 2
    @test any(isequal(B_user), get_brownians(sys))
end

# Tests that hybrid_model merges user jumps with reaction jumps.
let
    @species S(t) P(t)
    @parameters k1 k2

    user_jump = MT.ConstantRateJump(k2, [P ~ Pre(P) + 1])

    # System with Jump-scale reaction and user jump.
    @named rn = ReactionSystem(
        [Reaction(k1, [S], [P])],
        t, [S, P], [k1, k2];
        jumps = [user_jump]
    )
    rn = complete(rn)

    # Use hybrid_model with Jump scale for the reaction (index 1).
    sys = Catalyst.hybrid_model(rn;
        physical_scales = [1 => PhysicalScale.Jump],
        name = :test_sys
    )

    # Should have both reaction-generated and user jumps.
    @test length(get_jumps(sys)) == 2
end

### Test isequivalent with Brownians and Jumps ###

# Tests isequivalent correctly compares brownians.
let
    @brownians B1 B2
    @variables V(t)
    @parameters k

    @named rn1 = ReactionSystem([D(V) ~ -V + B1], t, [V], [k], [B1])
    @named rn2 = ReactionSystem([D(V) ~ -V + B1], t, [V], [k], [B1])
    @named rn3 = ReactionSystem([D(V) ~ -V + B2], t, [V], [k], [B2])

    # Same brownians → equivalent.
    @test isequivalent(rn1, rn2)

    # Different brownians → not equivalent.
    @test !isequivalent(rn1, rn3)
end

# Tests isequivalent correctly compares jumps.
let
    @species S(t)
    @parameters k1 k2

    jump1 = MT.ConstantRateJump(k1, [S ~ Pre(S) + 1])
    jump2 = MT.ConstantRateJump(k2, [S ~ Pre(S) + 2])

    @named rn1 = ReactionSystem(Reaction[], t, [S], [k1]; jumps = [jump1])
    @named rn2 = ReactionSystem(Reaction[], t, [S], [k1]; jumps = [jump1])
    @named rn3 = ReactionSystem(Reaction[], t, [S], [k1, k2]; jumps = [jump2])

    # Same jumps → equivalent.
    @test isequivalent(rn1, rn2)

    # Different jumps → not equivalent.
    @test !isequivalent(rn1, rn3)
end

### Test Subsystem Validation ###

# Tests that only ReactionSystems are allowed as subsystems.
let
    @variables V(t)
    @species S(t)
    @parameters k

    # Create a non-ReactionSystem (plain System).
    non_rs = MT.System([D(V) ~ -V], t; name = :non_rs)

    # Attempting to use it as a subsystem should error.
    @test_throws ErrorException ReactionSystem(
        [Reaction(k, nothing, [S])],
        t, [S], [k];
        systems = [non_rs],
        name = :parent
    )
end
