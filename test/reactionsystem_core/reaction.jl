### Preparations ###

# Fetch packages.
using Catalyst, Test
using Catalyst: get_symbolics
using ModelingToolkit: value, get_variables!

# Sets the default `t` to use.
t = default_t()

### Test Basic Accessors ###

# Tests the `get_variables` function.
let
    # Declare symbolic variables.
    @parameters k1 k2 n1 n2 η1 η2 p
    @species X(t) Y(t) Z(t)
    @variables A(t)
    
    # Create `Reaction`s.
    rx1 = Reaction(k1, [], [X])
    rx2 = Reaction(k1 + k2, [X], [Y], [1], [n1]; metadata = [:noise_scaling => η1])
    rx3 = Reaction(k1 + k2 + A, [X], [X, Y, Z], [1], [n1 + n2, 2, 1])
    rx4 = Reaction(X + t, [], [Y]; metadata = [:noise_scaling => η1 + η2])
    
    # Test `get_variables!`.
    @test issetequal(get_variables!([value(p)], rx1), [k1, X, p])
    @test issetequal(get_variables!([value(p)], rx2), [k1, k2, X, Y, n1, η1, p])
    @test issetequal(get_variables!([value(p)], rx3), [k1, k2, A, X, Y, Z, n1, n2, p])
    @test issetequal(get_variables!([value(p)], rx4), [X, t, Y, η1, η2, p])
    
    # Test `get_symbolics`.
    @test issetequal(get_symbolics(rx1), [k1, X])
    @test issetequal(get_symbolics(rx2), [k1, k2, X, Y, n1, η1])
    @test issetequal(get_symbolics(rx3), [k1, k2, A, X, Y, Z, n1, n2])
    @test issetequal(get_symbolics(rx4), [X, t, Y, η1, η2])
end

### Tests Metadata ###

# Tests creation.
# Tests basic accessor functions.
# Tests that repeated metadata entries are not permitted.
let
    @variables t
    @parameters k
    @species X(t) X2(t)

    metadata = [:noise_scaling => 0.0]
    r = Reaction(k, [X], [X2], [2], [1]; metadata=metadata)

    @test Catalyst.getmetadata_dict(r) == [:noise_scaling => 0.0]
    @test Catalyst.hasmetadata(r, :noise_scaling)
    @test !Catalyst.hasmetadata(r, :nonexisting_metadata)
    @test Catalyst.getmetadata(r, :noise_scaling) == 0.0

    metadata_repeated = [:noise_scaling => 0.0, :noise_scaling => 1.0, :metadata_entry => "unused"]
    @test_throws Exception Reaction(k, [X], [X2], [2], [1]; metadata=metadata_repeated)
end

# Tests accessors for system without metadata.
let
    @variables t
    @parameters k
    @species X(t) X2(t)

    metadata = Pair{Symbol,Any}[]
    r1 = Reaction(k, [X], [X2], [2], [1])
    r2 = Reaction(k, [X], [X2], [2], [1]; metadata=metadata)

    @test isequal(r1, r2)
    @test Catalyst.getmetadata_dict(r1) == Pair{Symbol,Any}[]
    @test !Catalyst.hasmetadata(r1, :md)
end

# Tests creation.
# Tests basic accessor functions.
# Tests various metadata types.
let
    @variables t
    @parameters k
    @species X(t) X2(t)

    metadata = Pair{Symbol,Any}[]
    push!(metadata, :md_1 => 1.0)
    push!(metadata, :md_2 => false)
    push!(metadata, :md_3 => "Hello world")
    push!(metadata, :md_4 => :sym)
    push!(metadata, :md_5 => X + X2^k -1)
    push!(metadata, :md_6 => (0.1, 2.0))
    r = Reaction(k, [X], [X2], [2], [1]; metadata=metadata)

    @test Catalyst.getmetadata_dict(r) isa Vector{Pair{Symbol,Any}}
    @test Catalyst.hasmetadata(r, :md_1)
    @test Catalyst.hasmetadata(r, :md_2)
    @test Catalyst.hasmetadata(r, :md_3)
    @test Catalyst.hasmetadata(r, :md_4)
    @test Catalyst.hasmetadata(r, :md_5)
    @test Catalyst.hasmetadata(r, :md_6)
    @test !Catalyst.hasmetadata(r, :md_8)
    
    @test isequal(Catalyst.getmetadata(r, :md_1), 1.0)
    @test isequal(Catalyst.getmetadata(r, :md_2), false)
    @test isequal(Catalyst.getmetadata(r, :md_3), "Hello world")
    @test isequal(Catalyst.getmetadata(r, :md_4), :sym)
    @test isequal(Catalyst.getmetadata(r, :md_5), X + X2^k -1)
    @test isequal(Catalyst.getmetadata(r, :md_6), (0.1, 2.0))
end

# Tests the noise scaling metadata.
let
    @variables t
    @parameters k η  
    @species X(t) X2(t)

    metadata = Pair{Symbol,Any}[]
    push!(metadata, :noise_scaling => η)
    r1 = Reaction(k, [X], [X2], [2], [1])
    r2 = Reaction(k, [X], [X2], [2], [1]; metadata=metadata)

    @test !Catalyst.has_noise_scaling(r1)
    @test Catalyst.has_noise_scaling(r2)
    @test_throws Exception Catalyst.get_noise_scaling(r1)
    @test isequal(Catalyst.get_noise_scaling(r2), η)
end