### Prepare Tests ###

# Fetch packages.
using Catalyst, Test
using Catalyst: get_rxs
using ModelingToolkit: getdefault, getdescription, get_metadata

# Creates missing getters for MTK metadata (can be removed once added to MTK).
getmisc(x) = SymbolicUtils.getmetadata(Symbolics.unwrap(x), ModelingToolkit.VariableMisc, nothing)
getinput(x) = SymbolicUtils.getmetadata(Symbolics.unwrap(x), ModelingToolkit.VariableInput, nothing)

# Sets the default `t` and `D` to use.
t = default_t()
D = default_time_deriv()


### Basic Test ###

# Checks for a simple reaction network (containing variables, equations, and observables).
# Checks that declaration via DSL works.
# Checks annotated and non-annotated files against manually written ones.
let 
    # Creates and serialises the model.
    rn = @reaction_network rn begin
        @observables X2 ~ 2X
        @equations D(V) ~ 1 - V
        d, X --> 0
    end
    save_reactionsystem("test_serialisation_annotated.jl", rn; safety_check = false)
    save_reactionsystem("test_serialisation.jl", rn; annotate = false, safety_check = false)

    # Checks equivalence.
    file_string_annotated = read("test_serialisation_annotated.jl", String)
    file_string = read("test_serialisation.jl", String)
    file_string_annotated_real = """let

    # Independent variable:
    @variables t

    # Parameters:
    ps = @parameters d

    # Species:
    sps = @species X(t)

    # Variables:
    vars = @variables V(t)

    # Reactions:
    rxs = [Reaction(d, [X], nothing, [1], nothing)]

    # Equations:
    eqs = [Differential(t)(V) ~ 1 - V]

    # Observables:
    @variables X2(t)
    observed = [X2 ~ 2X]

    # Declares ReactionSystem model:
    rs = ReactionSystem([rxs; eqs], t, [sps; vars], ps; name = :rn, observed)
    complete(rs)

    end"""
    file_string_real = """let

    @variables t
    ps = @parameters d
    sps = @species X(t)
    vars = @variables V(t)
    rxs = [Reaction(d, [X], nothing, [1], nothing)]
    eqs = [Differential(t)(V) ~ 1 - V]
    @variables X2(t)
    observed = [X2 ~ 2X]

    rs = ReactionSystem([rxs; eqs], t, [sps; vars], ps; name = :rn, observed)
    complete(rs)

    end"""
    @test file_string_annotated == file_string_annotated_real
    @test file_string == file_string_real

    # Deletes the files.
    rm("test_serialisation_annotated.jl")
    rm("test_serialisation.jl")
end

# Tests for hierarchical system created programmatically.
# Checks that the species, variables, and parameters have their non-default types, default values,
# and metadata recorded correctly (these are not considered for system equality is tested).
# Checks that various types (processed by the `x_2_string` function) are serialised properly.
# Checks that `ReactionSystem` and `Reaction` metadata fields are recorded properly.
let 
    # Prepares various stuff to add as metadata.
    bool_md = false
    int_md = 3
    float_md = 1.2
    rat_md = 4//5
    sym_md = :sym
    c_md = 'c'
    str_md = "A string"
    nothing_md = nothing
    @parameters s r
    symb_md = s
    expr_md = 2s + r^3
    pair_md = rat_md => symb_md 
    tup_md = (float_md, str_md, expr_md)
    vec_md = [float_md, sym_md, tup_md]
    dict_md = Dict([c_md => str_md, symb_md => vec_md])
    mat_md = [rat_md sym_md; symb_md tup_md]

    # Creates parameters, variables, and species (with various metadata and default values).
    @parameters begin
        a, [input=bool_md]
        b, [misc=int_md]
        c = float_md, [misc=rat_md]
        d1, [misc=c_md]
        d2, [description=str_md]
        e1, [misc=nothing_md]
        e2, [misc=symb_md]
    end
    @variables begin
        A(t) = float_md
        B(t), [misc=expr_md]
        C(t), [misc=pair_md]
        D1(t), [misc=tup_md]
        D2(t), [misc=vec_md]
        E1(t), [misc=dict_md]
        E2(t), [misc=mat_md]
    end
    @species begin
        X(t), [input=bool_md]
        Y(t), [misc=int_md]
        Z(t), [misc=float_md]
        V1(t), [description=str_md]
        V2(t), [misc=dict_md]
        W1(t), [misc=mat_md]
        W2(t) = float_md
    end

    # Creates the hierarchical model. Adds metadata to both reactions and the systems.
    # First reaction has `+ s + r` as that is easier than manually listing all symbolics.
    # (These needs to be part of the system somehow, as they are only added through the `misc` metadata)
    rxs1 = [
        Reaction(a + A + s + r, [X], [], metadata = [:misc => bool_md])
        Reaction(b + B, [Y], [], metadata = [:misc => int_md])
        Reaction(c + C, [Z], [], metadata = [:misc => sym_md])
        Reaction(d1 + D1, [V1], [], metadata = [:misc => str_md])
        Reaction(e1 + E1, [W1], [], metadata = [:misc => nothing_md])
    ]
    rxs2 = [
        Reaction(a + A, [X], [], metadata = [:misc => expr_md])
        Reaction(b + B, [Y], [], metadata = [:misc => tup_md])
        Reaction(c + C, [Z], [], metadata = [:misc => vec_md])
        Reaction(d2 + D2, [V2], [], metadata = [:misc => dict_md])
        Reaction(e2 + E2, [W2], [], metadata = [:misc => mat_md])
    ]
    @named rs2 = ReactionSystem(rxs2, t; metadata = dict_md)
    @named rs1 = ReactionSystem(rxs1, t; systems = [rs2], metadata = mat_md)
    rs = complete(rs1)

    # Loads the model and checks that it is correct. Removes the saved file
    save_reactionsystem("serialised_rs.jl", rs; safety_check = false)
    rs_loaded = include("../serialised_rs.jl")
    @test rs == rs_loaded
    rm("serialised_rs.jl")

    # Checks that parameters/species/variables metadata fields are correct.
    @test isequal(getinput(rs_loaded.a), bool_md)
    @test isequal(getmisc(rs_loaded.b), int_md)
    @test isequal(getdefault(rs_loaded.c), float_md)
    @test isequal(getmisc(rs_loaded.c), rat_md)
    @test isequal(getmisc(rs_loaded.d1), c_md)
    @test isequal(getdescription(rs_loaded.rs2.d2), str_md)
    @test isequal(getmisc(rs_loaded.e1), nothing_md)
    @test isequal(getmisc(rs_loaded.rs2.e2), symb_md)

    @test isequal(getdefault(rs.A), float_md)
    @test isequal(getmisc(rs_loaded.B), expr_md)
    @test isequal(getmisc(rs_loaded.C), pair_md)
    @test isequal(getmisc(rs_loaded.D1), tup_md)
    @test isequal(getmisc(rs_loaded.rs2.D2), vec_md)
    @test isequal(getmisc(rs_loaded.E1), dict_md)
    @test isequal(getmisc(rs_loaded.rs2.E2), mat_md)

    @test isequal(getinput(rs_loaded.X), bool_md)
    @test isequal(getmisc(rs_loaded.Y), int_md)
    @test isequal(getmisc(rs_loaded.Z), float_md)
    @test isequal(getdescription(rs_loaded.V1), str_md)
    @test isequal(getmisc(rs_loaded.rs2.V2), dict_md)
    @test isequal(getmisc(rs_loaded.W1), mat_md)
    @test isequal(getdefault(rs_loaded.rs2.W2), float_md)

    # Checks that `Reaction` metadata fields are correct.
    @test isequal(getmetadata(get_rxs(rs_loaded)[1], :misc), bool_md)
    @test isequal(getmetadata(get_rxs(rs_loaded)[2], :misc), int_md)
    @test isequal(getmetadata(get_rxs(rs_loaded)[3], :misc), sym_md)
    @test isequal(getmetadata(get_rxs(rs_loaded)[4], :misc), str_md)
    @test isequal(getmetadata(get_rxs(rs_loaded)[5], :misc), nothing_md)
    @test isequal(getmetadata(get_rxs(rs_loaded.rs2)[1], :misc), expr_md)
    @test isequal(getmetadata(get_rxs(rs_loaded.rs2)[2], :misc), tup_md)
    @test isequal(getmetadata(get_rxs(rs_loaded.rs2)[3], :misc), vec_md)
    @test isequal(getmetadata(get_rxs(rs_loaded.rs2)[4], :misc), dict_md)
    @test isequal(getmetadata(get_rxs(rs_loaded.rs2)[5], :misc), mat_md)

    # Checks that `ReactionSystem` metadata fields are correct.
    @test isequal(get_metadata(rs_loaded), mat_md)
    @test isequal(get_metadata(rs_loaded.rs2), dict_md)
end

# Checks systems where parameters/species/variables have complicated interdependency are correctly
# serialised.
# Checks for system with non-default independent variable.
let 
    # Prepares parameters/variables/species with complicated dependencies.
    @variables τ
    @parameters begin
        b = 3.0
        c
        f
    end
    @variables begin
        A(τ) = c
        B(τ) = c + A + f
        C(τ) = 2.0
        D(τ) = C
        G(τ)
    end
    @species begin
        Y(τ) = f
        Z(τ)
        U(τ) = G + Z
        V(τ)
    end
    @parameters begin
        a = G + D
        e = U
    end
    @variables begin
        E(τ) = G + Z
        F(τ) = f + G
    end
    @species begin
        X(τ) = f + a
    end
    @parameters d = X
    @species W(τ) = d + e

    # Creates model and checks it against serialised version.
    @named rs = ReactionSystem([], τ, [X, Y, Z, U, V, W, A, B, C, D, E, F, G], [a, b, c, d, e, f])
    save_reactionsystem("serialised_rs.jl", rs; safety_check = false)
    @test rs == include("../serialised_rs.jl")
    rm("serialised_rs.jl")
end


# Tests for multi-layered hierarchical system. Tests with spatial independent variables,
# variables, (differential and algebraic) equations, observables (continuous and discrete) events, 
# and with various species/variables/parameter/reaction/system metadata.
# Tests for complete and incomplete system.
let
    # Prepares spatial independent variables (technically not used and only supplied to systems).
    sivs = @variables x y z [description="A spatial independent variable."]

    # Prepares parameters, species, and variables.
    @parameters p d k1_1 k2_1 k1_2 k2_2 k1_3 k2_3 k1_4 k2_4 a b_1 b_2 b_3 b_4 η
    @parameters begin
        t_1 = 2.0
        t_2::Float64
        t_3, [description="A parameter."]
        t_4::Float32 = p, [description="A parameter."]
    end 
    @species X(t) X2_1(t) X2_2(t) X2_3(t) X2_4(t)=p [description="A species."]
    @variables A(t)=p [description="A variable."] B_1(t) B_2(t) B_3(t) B_4(t)

    # Prepares all equations.
    eqs_1 = [
        Reaction(p, [], [X]; metadata = [:description => "A reaction"]),
        Reaction(d, [X], []; metadata = [:noise_scaling => η]),
        Reaction(k1_1, [X], [X2_1], [2], [1]),
        Reaction(k2_1, [X2_1], [X], [1], [2]),
        D(A) ~ a - A,
        A + 2B_1^3 ~ b_1 * X
    ]
    eqs_2 = [
        Reaction(p, [], [X]; metadata = [:description => "A reaction"]),
        Reaction(d, [X], []; metadata = [:noise_scaling => η]),
        Reaction(k1_2, [X], [X2_2], [2], [1]),
        Reaction(k2_2, [X2_2], [X], [1], [2]),
        D(A) ~ a - A,
        A + 2B_2^3 ~ b_2 * X
    ]
    eqs_3 = [
        Reaction(p, [], [X]; metadata = [:description => "A reaction"]),
        Reaction(d, [X], []; metadata = [:noise_scaling => η]),
        Reaction(k1_3, [X], [X2_3], [2], [1]),
        Reaction(k2_3, [X2_3], [X], [1], [2]),
        D(A) ~ a - A,
        A + 2B_3^3 ~ b_3 * X
    ]
    eqs_4 = [
        Reaction(p, [], [X]; metadata = [:description => "A reaction"]),
        Reaction(d, [X], []; metadata = [:noise_scaling => η]),
        Reaction(k1_4, [X], [X2_4], [2], [1]),
        Reaction(k2_4, [X2_4], [X], [1], [2]),
        D(A) ~ a - A,
        A + 2B_4^3 ~ b_4 * X
    ]

    # Prepares all events.
    continuous_events_1 = [(A ~ t_1) => [A ~ A + 2.0, X ~ X/2]]
    continuous_events_2 = [(A ~ t_2) => [A ~ A + 2.0, X ~ X/2]]
    continuous_events_3 = [(A ~ t_3) => [A ~ A + 2.0, X ~ X/2]]
    continuous_events_4 = [(A ~ t_4) => [A ~ A + 2.0, X ~ X/2]]
    discrete_events_1 = [
        10.0 => [X2_1 ~ X2_1 + 1.0]
        [5.0, 10.0] => [b_1 ~ 2 * b_1]
        (X > 5.0) => [X2_1 ~ X2_1 + 1.0, X ~ X - 1]
    ]
    discrete_events_2 = [
        10.0 => [X2_2 ~ X2_2 + 1.0]
        [5.0, 10.0] => [b_2 ~ 2 * b_2]
        (X > 5.0) => [X2_2 ~ X2_2 + 1.0, X ~ X - 1]
    ]
    discrete_events_3 = [
        10.0 => [X2_3 ~ X2_3 + 1.0]
        [5.0, 10.0] => [b_3 ~ 2 * b_3]
        (X > 5.0) => [X2_3 ~ X2_3 + 1.0, X ~ X - 1]
    ]
    discrete_events_4 = [
        10.0 => [X2_4 ~ X2_4 + 1.0]
        [5.0, 10.0] => [b_4 ~ 2 * b_4]
        (X > 5.0) => [X2_4 ~ X2_4 + 1.0, X ~ X - 1]
    ]

    # Creates the systems.
    @named rs_4 = ReactionSystem(eqs_4, t; continuous_events = continuous_events_4,
                                discrete_events = discrete_events_4, spatial_ivs = sivs, 
                                metadata = "System 4", systems = [])
    @named rs_2 = ReactionSystem(eqs_2, t; continuous_events = continuous_events_2,
                                discrete_events = discrete_events_2, spatial_ivs = sivs, 
                                metadata = "System 2", systems = [])
    @named rs_3 = ReactionSystem(eqs_3, t; continuous_events = continuous_events_3,
                                discrete_events = discrete_events_3, spatial_ivs = sivs, 
                                metadata = "System 3", systems = [rs_4])
    @named rs_1 = ReactionSystem(eqs_1, t; continuous_events = continuous_events_1,
                                discrete_events = discrete_events_1, spatial_ivs = sivs, 
                                metadata = "System 1", systems = [rs_2, rs_3])
    rs = complete(rs_1)

    # Checks that the correct system is saved (both complete and incomplete ones).
    save_reactionsystem("serialised_rs_incomplete.jl", rs_1; safety_check = false)
    @test isequal(rs_1, include("../serialised_rs_incomplete.jl"))
    save_reactionsystem("serialised_rs_complete.jl", rs; safety_check = false)
    @test isequal(rs, include("../serialised_rs_complete.jl"))
    rm("serialised_rs_incomplete.jl")
    rm("serialised_rs_complete.jl")
end

# Tests for (slightly more) complicate system created via the DSL.
# Tests for cases where the number of input is untested (i.e. multiple observables and continuous
# events, but single equations and discrete events).
# Tests with and without `safety_check`.
let 
    # Declares the model.
    rs = @reaction_network begin
        @equations D(V) ~ 1 - V
        @continuous_events begin
            [X ~ 5.0] => [X ~ X + 1.0]
            [X ~ 20.0] => [X ~ X - 1.0]
        end
        @discrete_events 5.0 => [d ~ d/2]
        d, X --> 0
    end

    # Checks that serialisation works.
    save_reactionsystem("serialised_rs_1.jl", rs)
    save_reactionsystem("serialised_rs_2.jl", rs; safety_check = false)
    isequal(rs, include("../serialised_rs_1.jl"))
    isequal(rs, include("../serialised_rs_2.jl"))
    rm("serialised_rs_1.jl")
    rm("serialised_rs_2.jl")
end

# Tests for system where species depends on multiple independent variables.
# Tests for system where variables depends on multiple independent variables.
let
    rs = @reaction_network begin
        @ivs t x y z
        @parameters p
        @species X(t,x,y) Y(t,x,y) XY(t,x,y) Z(t,x,y)
        @variables V(t,x,z)
        (kB,kD), X + Y <--> XY
    end
    save_reactionsystem("serialised_rs.jl", rs)
    @test ModelingToolkit.isequal(rs, include("../serialised_rs.jl"))
    rm("serialised_rs.jl")
end


### Other Tests ###

# Tests that an error is generated when non-`ReactionSystem` subs-systems are used.
let
    @variables V(t)
    @species X(t)
    @parameters p d V_max

    rxs = [
        Reaction(p, [], [X]),
        Reaction(d, [X], [])
    ]
    eq = D(V) ~ V_max - V

    @named osys = ODESystem([eq], t)
    @named rs = ReactionSystem(rxs, t; systems = [osys])
    @test_throws Exception save_reactionsystem("failed_serialisation.jl", rs)
end

# Checks that completeness is recorded correctly.
# Checks without turning off the `safety_check` option.
let 
    # Checks for complete system.
    rs_complete = @reaction_network begin
        (p,d), 0 <--> X
    end
    save_reactionsystem("serialised_rs_complete.jl", rs_complete)
    rs_complete_loaded = include("../serialised_rs_complete.jl")
    @test ModelingToolkit.iscomplete(rs_complete_loaded)
    rm("serialised_rs_complete.jl")

    # Checks for non-complete system.
    rs_incomplete = @network_component begin
        (p,d), 0 <--> X
    end
    save_reactionsystem("serialised_rs_incomplete.jl", rs_incomplete)
    rs_incomplete_loaded = include("../serialised_rs_incomplete.jl")
    @test !ModelingToolkit.iscomplete(rs_incomplete_loaded)
    rm("serialised_rs_incomplete.jl")
end