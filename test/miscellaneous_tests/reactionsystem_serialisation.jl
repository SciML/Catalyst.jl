### Prepare Tests ###

# Fetch packages.
using Catalyst

# Sets the default `t` and `D` to use.
t = default_t()
D = default_time_deriv()


### Basic Test ###

# Checks for a simple reaction network (containing variables, equations, and observables).
# Checks that declaration via DSL works.
# Checks annotated and non-annotated files against manually written ones.
let 
    # Creates and serialises the model.
    rn = @reaction_system begin
        @observables X2 ~ 2X
        @equations D(V) ~ 1 - V
        (p,d), 0 <--> X
    end
    save_reaction_network("test_serialisation_annotated.jl", rn)
    save_reaction_network("test_serialisation.jl", rn; annotate = false)

    # Checks equivalence.
    file_string_annotated = read("test_serialisation_annotated.jl", String)
    file_string = read("test_serialisation.jl", String)
    file_string_annotated_real = ""
    file_string_real = ""
    @test file_string_annotated == file_string_annotated_real
    @test file == file_string_real

    # Deletes the files.
    rm("test_serialisation_annotated.jl")
    rm("test_serialisation.jl")
end

# Tests for hierarchical system created programmatically.
# Checks that the species, variables, and parameters have their non-default types, default values,
# and metadata recorded correctly (these are not considered for system equality is tested).
let 
    
end

# Tests for complicated hierarchical system. Tests with non-default independent variable, 
# spatial independent variables, variables, (differential and algebraic) equations, observables
# (continuous and discrete) events, and with various species/variables/parameter metadata.
let

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
    @test_throws Exception save_reaction_network("failed_serialisation.jl", rs)
end

# Checks that completeness is recorded correctly.
let 
    # Checks for complete system.
    rs_complete = @reaction_network begin
        (p,d), 0 <--> X
    end
    save_reaction_network("serialised_rs_complete.jl", rs_complete)
    rs_complete_loaded = include("serialised_rs_complete.jl")
    @test ModelingToolkit.iscomplete(rs_complete_loaded)
    rn("serialised_rs_complete.jl")

    # Checks for non-complete system.
    rs_incomplete = @network_component begin
        (p,d), 0 <--> X
    end
    save_reaction_network("serialised_rs_incomplete.jl", rs_incomplete)
    rs_incomplete_loaded = include("serialised_rs_incomplete.jl")
    @test !ModelingToolkit.iscomplete(rs_incomplete_loaded)
    rn("serialised_rs_incomplete.jl")
end