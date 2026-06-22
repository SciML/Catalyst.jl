# Tests for alias elimination functionality.
using Catalyst, Test
import ModelingToolkitBase as MT
import DiffEqBase
using JumpProcesses

const t = Catalyst.default_t()

### AliasClass and Validation Tests ###

@testset "AliasClass enum" begin
    @test AliasClass.Species isa AliasClass.T
    @test AliasClass.BCSpecies isa AliasClass.T
    @test AliasClass.ConstantSpecies isa AliasClass.T
    @test AliasClass.CompoundSpecies isa AliasClass.T
    @test AliasClass.Unknown isa AliasClass.T
    @test AliasClass.Parameter isa AliasClass.T
    @test AliasClass.Brownian isa AliasClass.T
    @test AliasClass.Poissonian isa AliasClass.T
    @test AliasClass.Bound isa AliasClass.T
    @test AliasClass.Observable isa AliasClass.T
    @test AliasClass.Unsupported isa AliasClass.T
end

@testset "alias_compatible" begin
    @test Catalyst.alias_compatible(AliasClass.Species, AliasClass.Species)
    @test Catalyst.alias_compatible(AliasClass.BCSpecies, AliasClass.BCSpecies)
    @test Catalyst.alias_compatible(AliasClass.ConstantSpecies, AliasClass.ConstantSpecies)
    @test Catalyst.alias_compatible(AliasClass.Unknown, AliasClass.Unknown)
    @test Catalyst.alias_compatible(AliasClass.Parameter, AliasClass.Parameter)
    @test !Catalyst.alias_compatible(AliasClass.Species, AliasClass.Parameter)
    @test !Catalyst.alias_compatible(AliasClass.Brownian, AliasClass.Brownian)
    @test !Catalyst.alias_compatible(AliasClass.CompoundSpecies, AliasClass.CompoundSpecies)
end

@testset "check_aliases (constructor checks)" begin
    @species A(t) B(t)
    @parameters k1 k2
    Catalyst.check_aliases(Equation[])
    Catalyst.check_aliases([A ~ B])
    Catalyst.check_aliases([A ~ B, k1 ~ k2])
    @test_throws ErrorException Catalyst.check_aliases([A ~ A])
    @test_throws ErrorException Catalyst.check_aliases([A ~ B, A ~ B])
end

### Alias Resolution Tests ###

@testset "Alias resolution" begin
    @species A(t) B(t) C(t) D(t)
    @parameters k1 k2

    # Basic species alias
    @named rn = ReactionSystem([Reaction(k1, [A], [B])], t, [A, B], [k1, k2];
        aliases = [A ~ B])
    rn_c = complete(rn)
    flat = Catalyst.flatten(rn_c)
    result = Catalyst.validate_and_resolve_aliases(Catalyst.get_aliases(flat), flat)
    @test isequal(result.unknown_submap[A], B)

    # Chain aliasing
    @named rn_chain = ReactionSystem(
        [Reaction(k1, [A], [C])], t, [A, B, C], [k1];
        aliases = [A ~ B, B ~ C])
    flat_chain = Catalyst.flatten(complete(rn_chain))
    result_chain = Catalyst.validate_and_resolve_aliases(Catalyst.get_aliases(flat_chain), flat_chain)
    @test isequal(result_chain.unknown_submap[A], C)
    @test isequal(result_chain.unknown_submap[B], C)

    # Fan-in + chain
    @named rn_fan = ReactionSystem(
        [Reaction(k1, [A], [C])], t, [A, B, C, D], [k1];
        aliases = [A ~ C, B ~ C, D ~ B])
    flat_fan = Catalyst.flatten(complete(rn_fan))
    result_fan = Catalyst.validate_and_resolve_aliases(Catalyst.get_aliases(flat_fan), flat_fan)
    @test isequal(result_fan.unknown_submap[A], C)
    @test isequal(result_fan.unknown_submap[B], C)
    @test isequal(result_fan.unknown_submap[D], C)

    # Cycle detection
    @named rn_cycle = ReactionSystem(
        [Reaction(k1, [A], [C])], t, [A, B, C], [k1];
        aliases = [A ~ B, B ~ A])
    flat_cycle = Catalyst.flatten(complete(rn_cycle))
    @test_throws ErrorException Catalyst.validate_and_resolve_aliases(
        Catalyst.get_aliases(flat_cycle), flat_cycle)

    # Parameter alias
    @named rn_param = ReactionSystem(
        [Reaction(k1, [A], [B])], t, [A, B], [k1, k2];
        aliases = [k1 ~ k2])
    flat_param = Catalyst.flatten(complete(rn_param))
    result_param = Catalyst.validate_and_resolve_aliases(
        Catalyst.get_aliases(flat_param), flat_param)
    @test isequal(result_param.param_submap[k1], k2)
end

### Alias Elimination Tests ###

@testset "Alias elimination" begin
    @species A(t) B(t) C(t)
    @parameters k1 k2

    # Basic species elimination
    @named rn = ReactionSystem(
        [Reaction(k1, [A], [B]), Reaction(k2, [B], [C])], t, [A, B, C], [k1, k2];
        aliases = [A ~ B])
    flat = Catalyst.flatten(complete(rn))
    elim = Catalyst.eliminate_aliases(flat)
    @test !any(isequal(A), Catalyst.get_unknowns(elim))
    @test any(eq -> isequal(eq.lhs, A), Catalyst.get_observed(elim))
    @test isempty(Catalyst.get_aliases(elim))

    # Parameter elimination
    @named rn_p = ReactionSystem(
        [Reaction(k1, [A], [B])], t, [A, B], [k1, k2];
        aliases = [k1 ~ k2])
    flat_p = Catalyst.flatten(complete(rn_p))
    elim_p = Catalyst.eliminate_aliases(flat_p)
    @test !any(isequal(k1), Catalyst.get_ps(elim_p))
    @test any(isequal(k2), Catalyst.get_ps(elim_p))

    # No-op reaction removal (A→B with A~B becomes B→B, dropped)
    @named rn_noop = ReactionSystem(
        [Reaction(k1, [A], [B]), Reaction(k2, [B], [C])], t, [A, B, C], [k1, k2];
        aliases = [A ~ B])
    elim_noop = Catalyst.eliminate_aliases(Catalyst.flatten(complete(rn_noop)))
    @test length(Catalyst.get_rxs(elim_noop)) == 1

    # Stoichiometry merging: 2A + B → C with A~B becomes 3B → C
    @named rn_merge = ReactionSystem(
        [Reaction(k1, [A, B], [C], [2, 1], [1])], t, [A, B, C], [k1];
        aliases = [A ~ B])
    elim_merge = Catalyst.eliminate_aliases(Catalyst.flatten(complete(rn_merge)))
    rx = first(Catalyst.get_rxs(elim_merge))
    @test length(rx.substrates) == 1
    @test isequal(rx.substrates[1], B)
    @test rx.substoich[1] == 3

    # Initial conditions transfer
    @named rn_ic = ReactionSystem(
        [Reaction(k1, [A], [B])], t, [A, B], [k1];
        aliases = [A ~ B], initial_conditions = Dict(A => 5.0))
    elim_ic = Catalyst.eliminate_aliases(Catalyst.flatten(complete(rn_ic)))
    ics = MT.initial_conditions(elim_ic)
    @test !haskey(ics, A)
    @test haskey(ics, B)

    # Initial conditions conflict
    @named rn_conflict = ReactionSystem(
        [Reaction(k1, [A], [B])], t, [A, B], [k1];
        aliases = [A ~ B], initial_conditions = Dict(A => 5.0, B => 10.0))
    @test_throws ErrorException Catalyst.eliminate_aliases(
        Catalyst.flatten(complete(rn_conflict)))

    # Alias maps in metadata
    @named rn_meta = ReactionSystem(
        [Reaction(k1, [A], [B])], t, [A, B], [k1, k2];
        aliases = [A ~ B, k1 ~ k2])
    elim_meta = Catalyst.eliminate_aliases(Catalyst.flatten(complete(rn_meta)))
    maps = Catalyst.get_alias_maps(elim_meta)
    @test maps !== nothing
    @test isequal(maps.unknown_submap[A], B)

    # Input remapping
    u0 = Dict(A => 1.0)
    remapped = Catalyst.remap_alias_inputs(u0, elim_meta)
    @test !haskey(remapped, A)
    @test haskey(remapped, B)
end

### DSL Tests ###

@testset "DSL @aliases" begin
    rn = @reaction_network begin
        @aliases A ~ B
        @species A(t) B(t) C(t)
        @parameters k
        k, A --> C
        k, B --> C
    end
    @test aliases_present(rn)
    @test length(Catalyst.get_aliases(rn)) == 1

    rn2 = @reaction_network begin
        @aliases begin
            A ~ B
            k1 ~ k2
        end
        @species A(t) B(t) C(t)
        @parameters k1 k2
        k1, A --> C
        k2, B --> C
    end
    @test length(Catalyst.get_aliases(rn2)) == 2
end

### Conversion Pipeline Tests ###

@testset "Conversion pipeline with aliases" begin
    @species A(t) B(t) C(t)
    @parameters k1 k2

    @named rn = ReactionSystem(
        [Reaction(k1, [A], [B]), Reaction(k2, [B], [C])], t, [A, B, C], [k1, k2];
        aliases = [A ~ B])
    rn_c = complete(rn)

    # Model conversions
    @test ode_model(rn_c) isa MT.System
    @test sde_model(rn_c) isa MT.System
    @test jump_model(rn_c) isa MT.System
    @test hybrid_model(rn_c; default_scale = PhysicalScale.ODE) isa MT.System

    # Problem construction
    prob = ODEProblem(rn_c, [B => 1.0, C => 0.0], (0.0, 1.0), [k1 => 0.5, k2 => 0.3])
    @test prob isa ODEProblem
end

### Compose / Extend / Flatten Tests ###

@testset "Compose/extend/flatten with aliases" begin
    @species A(t) B(t) C(t) D(t)
    @parameters k1 k2 k3

    # Flatten collects aliases
    @named sub1 = ReactionSystem([Reaction(k1, [A], [B])], t, [A, B], [k1];
        aliases = [A ~ B])
    @named sub2 = ReactionSystem([Reaction(k2, [C], [D])], t, [C, D], [k2])
    @named parent = ReactionSystem(Reaction[], t; systems = [sub1, sub2])
    flat = Catalyst.flatten(parent)
    @test length(Catalyst.get_aliases(flat)) == 1

    # Compose with aliases keyword
    @named sys1 = ReactionSystem([Reaction(k1, [A], [B])], t, [A, B], [k1])
    @named sys2 = ReactionSystem([Reaction(k2, [C], [D])], t, [C, D], [k2])
    composed = compose(sys1, [sys2]; aliases = [A ~ C])
    @test length(Catalyst.get_aliases(composed)) == 1

    # Extend unions aliases
    @named ext1 = ReactionSystem([Reaction(k1, [A], [B])], t, [A, B], [k1];
        aliases = [A ~ B])
    @named ext2 = ReactionSystem([Reaction(k2, [C], [D])], t, [C, D], [k2, k3];
        aliases = [k2 ~ k3])
    extended = extend(ext1, ext2)
    @test length(Catalyst.get_aliases(extended)) == 2
end

### isequivalent Tests ###

@testset "isequivalent with aliases" begin
    @species A(t) B(t)
    @parameters k1 k2

    @named rn1 = ReactionSystem([Reaction(k1, [A], [B])], t, [A, B], [k1, k2];
        aliases = [A ~ B])
    @named rn2 = ReactionSystem([Reaction(k1, [A], [B])], t, [A, B], [k1, k2];
        aliases = [A ~ B])
    @test Catalyst.isequivalent(rn1, rn2)

    @named rn3 = ReactionSystem([Reaction(k1, [A], [B])], t, [A, B], [k1, k2];
        aliases = [k1 ~ k2])
    @test !Catalyst.isequivalent(rn1, rn3)

    @named rn4 = ReactionSystem([Reaction(k1, [A], [B])], t, [A, B], [k1, k2])
    @test !Catalyst.isequivalent(rn1, rn4)
end

### Serialization Guard Test ###

@testset "Serialization guard" begin
    @species A(t) B(t)
    @parameters k1
    @named rn = ReactionSystem([Reaction(k1, [A], [B])], t, [A, B], [k1];
        aliases = [A ~ B])
    @test_throws ErrorException Catalyst.save_reactionsystem(tempname() * ".jl", rn)
end

### Input Remapping Tests ###

@testset "remap_alias_inputs formats" begin
    @species A(t) B(t) C(t)
    @parameters k1 k2

    @named rn = ReactionSystem(
        [Reaction(k1, [A], [B])], t, [A, B], [k1, k2];
        aliases = [A ~ B, k1 ~ k2])
    elim = Catalyst.eliminate_aliases(Catalyst.flatten(complete(rn)))

    # Vector{Pair} input
    u0_vec = [A => 1.0, C => 2.0]
    remapped_vec = Catalyst.remap_alias_inputs(u0_vec, elim)
    @test remapped_vec isa Vector
    @test any(p -> isequal(p[1], B), remapped_vec)
    @test any(p -> isequal(p[1], C), remapped_vec)

    # Dict input
    u0_dict = Dict(A => 1.0, C => 2.0)
    remapped_dict = Catalyst.remap_alias_inputs(u0_dict, elim)
    @test remapped_dict isa Dict
    @test haskey(remapped_dict, B)
    @test !haskey(remapped_dict, A)

    # NullParameters no-op
    @test Catalyst.remap_alias_inputs(
        DiffEqBase.NullParameters(), elim) isa DiffEqBase.NullParameters

    # Conflict: Dict with both eliminated and canonical keys (different values)
    @test_throws ErrorException Catalyst.remap_alias_inputs(Dict(A => 1.0, B => 5.0), elim)

    # Conflict: Vector{Pair} with both eliminated and canonical keys (different values)
    @test_throws ErrorException Catalyst.remap_alias_inputs([A => 1.0, B => 5.0], elim)
end

### Constructor-Level Remapping Tests ###

@testset "Constructor remapping with eliminated symbols" begin
    @species A(t) B(t) C(t)
    @parameters k1 k2

    @named rn = ReactionSystem(
        [Reaction(k1, [A], [B]), Reaction(k2, [B], [C])], t, [A, B, C], [k1, k2];
        aliases = [A ~ B])
    rn_c = complete(rn)

    # ODEProblem with eliminated symbol A in u0
    prob = ODEProblem(rn_c, [A => 1.0, C => 0.0], (0.0, 1.0), [k1 => 0.5, k2 => 0.3])
    @test prob isa ODEProblem

    # ODEProblem with NullParameters (system with no free parameters)
    @named rn_np = ReactionSystem(
        [Reaction(1.0, [A], [B])], t, [A, B], [];
        aliases = [A ~ B])
    prob_np = ODEProblem(complete(rn_np), [A => 1.0], (0.0, 1.0))
    @test prob_np isa ODEProblem
end

### Non-Symbolic Jump Guard ###

@testset "Non-symbolic jump constructor guard" begin
    @species A(t) B(t)
    @parameters k1
    # Imperative (function-based) CRJ should be rejected at construction
    crj = ConstantRateJump((u, p, t) -> p[1], integrator -> (integrator.u[1] -= 1))
    @test_throws ErrorException ReactionSystem(
        Reaction[], Catalyst.default_t(), [A, B], [k1]; jumps = [crj], name = :test)
end

### MassActionJump + Aliases Guard ###

@testset "MassActionJump + aliases error guard" begin
    @species A(t) B(t)
    @parameters k1
    maj = MassActionJump(k1, [A => 1], [B => 1])
    @test_throws ErrorException Catalyst.substitute_jump(maj, Dict(A => B))
end

### Additional Constructor-Level Remapping Tests ###

@testset "SDEProblem with eliminated-symbol p" begin
    @species A(t) B(t)
    @parameters k1 k2
    # Use parameter alias (no species alias avoids no-op reaction / noise matrix size issue)
    @named rn = ReactionSystem(
        [Reaction(k1, [A], [B]), Reaction(k2, [B], [A])], t, [A, B], [k1, k2];
        aliases = [k1 ~ k2])
    prob = SDEProblem(complete(rn), [A => 1.0, B => 0.0], (0.0, 1.0), [k1 => 0.5])
    @test prob isa SDEProblem
end

@testset "JumpProblem with eliminated-symbol u0" begin
    @species A(t) B(t) C(t)
    @parameters k1 k2
    @named rn = ReactionSystem(
        [Reaction(k1, [A], [B]), Reaction(k2, [B], [C])], t, [A, B, C], [k1, k2];
        aliases = [A ~ B])
    prob = JumpProblem(complete(rn), [A => 10, C => 0], (0.0, 1.0), [k1 => 0.5, k2 => 0.3])
    @test prob isa JumpProblem
end

@testset "HybridProblem with eliminated-symbol u0" begin
    @species A(t) B(t) C(t)
    @parameters k1 k2
    @named rn = ReactionSystem(
        [Reaction(k1, [A], [B]), Reaction(k2, [B], [C])], t, [A, B, C], [k1, k2];
        aliases = [A ~ B])
    prob = HybridProblem(complete(rn), [A => 1.0, C => 0.0], (0.0, 1.0),
        [k1 => 0.5, k2 => 0.3]; default_scale = PhysicalScale.ODE)
    @test prob isa ODEProblem
end

@testset "SteadyStateProblem with eliminated-symbol u0" begin
    @species A(t) B(t) C(t)
    @parameters k1 k2
    @named rn = ReactionSystem(
        [Reaction(k1, [A], [B]), Reaction(k2, [B], [C])], t, [A, B, C], [k1, k2];
        aliases = [A ~ B])
    prob = SteadyStateProblem(complete(rn), [A => 1.0, C => 0.0], [k1 => 0.5, k2 => 0.3])
    @test prob isa SteadyStateProblem
end

### eliminate_aliases=false Tests ###

@testset "eliminate_aliases=false" begin
    @species A(t) B(t) C(t)
    @parameters k1 k2
    @named rn = ReactionSystem(
        [Reaction(k1, [A], [B]), Reaction(k2, [B], [C])], t, [A, B, C], [k1, k2];
        aliases = [A ~ B])
    rn_c = complete(rn)

    # Without mtkcompile should error (alias-generated constraints need DAE handling)
    @test_throws ErrorException ODEProblem(rn_c, [B => 1.0, C => 0.0], (0.0, 1.0),
        [k1 => 0.5, k2 => 0.3]; eliminate_aliases = false)
end

### Binding Value Substitution Test ###

@testset "Binding value substitution through aliases" begin
    @species A(t) B(t)
    @parameters k1 k2

    # Test initial_conditions value substitution
    @named rn_ic = ReactionSystem(
        [Reaction(k1, [A], [B])], t, [A, B], [k1, k2];
        aliases = [k1 ~ k2],
        initial_conditions = Dict(A => 2 * k1))
    elim_ic = Catalyst.eliminate_aliases(Catalyst.flatten(complete(rn_ic)))
    ics = MT.initial_conditions(elim_ic)
    @test haskey(ics, A)
    ic_vars = Symbolics.get_variables(ics[A])
    @test !any(isequal(k1), ic_vars)
    @test any(isequal(k2), ic_vars)

    # Test actual bindings-map rebuilding: create a species with a symbolic default
    # (which MTKBase extracts into bindings via process_variables!).
    @species X(t) = 3 * k1
    @named rn_bind = ReactionSystem(
        [Reaction(k1, [X], [B])], t, [X, B], [k1, k2];
        aliases = [k1 ~ k2])
    flat_bind = Catalyst.flatten(complete(rn_bind))
    elim_bind = Catalyst.eliminate_aliases(flat_bind)
    bindings = MT.get_bindings(elim_bind)
    # X's binding should exist and reference k2 (canonical), not k1 (eliminated)
    @test haskey(bindings, X)
    bind_vars = Symbolics.get_variables(bindings[X])
    @test !any(isequal(k1), bind_vars)
    @test any(isequal(k2), bind_vars)
end

### Discrete Event Alias Substitution Test ###

@testset "Discrete event with aliased parameter" begin
    @species A(t) B(t)
    @parameters k1 k2

    # Create a system with a discrete event that references an aliased parameter.
    # The event affect should have k1 substituted to k2 after elimination.
    devt = MT.SymbolicDiscreteCallback(
        1.0, [A ~ A + k1])
    @named rn = ReactionSystem(
        [Reaction(k1, [A], [B])], t, [A, B], [k1, k2];
        aliases = [k1 ~ k2],
        discrete_events = [devt])
    flat = Catalyst.flatten(complete(rn))
    elim = Catalyst.eliminate_aliases(flat)

    # After elimination, discrete events should reference k2, not k1
    devents = MT.get_discrete_events(elim)
    @test length(devents) == 1
    # Check the affect equation: should be A ~ A + k2 (not k1)
    affect_eqs = devents[1].affect.affect
    affect_vars = union([Symbolics.get_variables(eq.rhs) for eq in affect_eqs]...)
    @test !any(isequal(k1), affect_vars)
    @test any(isequal(k2), affect_vars)
end

### Alias Autodiscovery Tests ###

@testset "Autodiscovery of species from aliases" begin
    # B only appears in the alias, not in any reaction or the explicit species list.
    # It should be autodiscovered as a species from the alias equation.
    @species A(t) B(t) C(t)
    @parameters k
    @named rn = ReactionSystem(
        [Reaction(k, [A], [C])];
        aliases = [A ~ B])
    rn = complete(Catalyst.flatten(rn))
    @test any(isequal(B), Catalyst.get_species(rn))
    elim = Catalyst.eliminate_aliases(rn)
    @test !any(isequal(A), Catalyst.get_species(elim))
    @test any(isequal(B), Catalyst.get_species(elim))
end

@testset "Autodiscovery of parameters from aliases" begin
    # k2 only appears in the alias, not in any reaction or the explicit parameter list.
    # It should be autodiscovered as a parameter from the alias equation.
    @species A(t) B(t)
    @parameters k1 k2
    @named rn = ReactionSystem(
        [Reaction(k1, [A], [B])];
        aliases = [k1 ~ k2])
    rn = complete(Catalyst.flatten(rn))
    @test any(isequal(k2), parameters(rn))
    elim = Catalyst.eliminate_aliases(rn)
    @test !any(isequal(k1), parameters(elim))
    @test any(isequal(k2), parameters(elim))
end
