### Preparations ###

# Fetch packages.
using Catalyst, Test

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

    @test get_metadata_dict(r) == [:noise_scaling => 0.0]
    @test has_metadata(r, :noise_scaling)
    @test !has_metadata(r, :nonexisting_metadata)
    @test get_metadata(r, :noise_scaling) == 0.0

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
    @test get_metadata_dict(r1) == Pair{Symbol,Any}[]
    @test !has_metadata(r1, :md)
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

    @test get_metadata_dict(r) isa Vector{Pair{Symbol,Any}}
    @test has_metadata(r, :md_1)
    @test has_metadata(r, :md_2)
    @test has_metadata(r, :md_3)
    @test has_metadata(r, :md_4)
    @test has_metadata(r, :md_5)
    @test has_metadata(r, :md_6)
    @test !has_metadata(r, :md_8)
    
    @test isequal(get_metadata(r, :md_1), 1.0)
    @test isequal(get_metadata(r, :md_2), false)
    @test isequal(get_metadata(r, :md_3), "Hello world")
    @test isequal(get_metadata(r, :md_4), :sym)
    @test isequal(get_metadata(r, :md_5), X + X2^k -1)
    @test isequal(get_metadata(r, :md_6), (0.1, 2.0))
end