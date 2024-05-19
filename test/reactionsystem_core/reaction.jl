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

    r1 = Reaction(k, [X], [X2], [2], [1])
    r2 = Reaction(k, [X], [X2], [2], [1]; metadata=[:noise_scaling => η])

    @test !Catalyst.has_noise_scaling(r1)
    @test Catalyst.has_noise_scaling(r2)
    @test_throws Exception Catalyst.get_noise_scaling(r1)
    @test isequal(Catalyst.get_noise_scaling(r2), η)
end

# Tests the description metadata.
let
    @variables t
    @parameters k η  
    @species X(t) X2(t)

    r1 = Reaction(k, [X], [X2], [2], [1])
    r2 = Reaction(k, [X], [X2], [2], [1]; metadata=[:description => "A reaction"])

    @test !Catalyst.has_description(r1)
    @test Catalyst.has_description(r2)
    @test_throws Exception Catalyst.get_description(r1)
    @test isequal(Catalyst.get_description(r2), "A reaction")
end

# Tests the misc metadata.
let
    @variables t
    @parameters k η  
    @species X(t) X2(t)

    r1 = Reaction(k, [X], [X2], [2], [1])
    r2 = Reaction(k, [X], [X2], [2], [1]; metadata=[:misc => ('M', :M)])

    @test !Catalyst.has_misc(r1)
    @test Catalyst.has_misc(r2)
    @test_throws Exception Catalyst.get_misc(r1)
    @test isequal(Catalyst.get_misc(r2), ('C', :C))
end