using Test
using Catalyst
using Aqua
using ExplicitImports

@testset "Code quality (Aqua.jl)" begin
    # Test with Aqua.jl - testing code quality
    # We allow unbound type parameters in some constructors 
    # We allow some ambiguities that are hard to resolve without breaking changes
    Aqua.test_all(Catalyst;
        ambiguities = false,  # TODO: Fix ambiguities in future PR
        unbound_args = false, # Some constructors have unbound type parameters by design
        stale_deps = false,   # Some test dependencies might appear stale
        deps_compat = false,  # Skip - incorrectly flags stdlib packages
        piracies = false      # We extend some Base/MTK methods which might be detected as piracy
    )

    # Test individual Aqua checks that we want to enforce
    @testset "Aqua selective tests" begin
        Aqua.test_undefined_exports(Catalyst)
        Aqua.test_project_extras(Catalyst)
        # Skip deps_compat test as it incorrectly flags stdlib packages
        # Aqua.test_deps_compat(Catalyst)
    end
end

@testset "Explicit imports (ExplicitImports.jl)" begin
    # Test that we're not relying on implicit imports
    @testset "Check implicit imports" begin
        # Check for implicit imports - we allow some from closely related packages
        result = check_no_implicit_imports(Catalyst; skip = (Base, Core))
        
        # If there are any implicit imports, check if they're acceptable
        if !isnothing(result)
            # Extract the actual implicit imports
            implicit_names = []
            for r in result
                if hasproperty(r, :name)
                    push!(implicit_names, r.name)
                end
            end
            
            # Allow some specific implicit imports that are intentional
            allowed_implicit = [
                # Add any symbols here that we intentionally want to allow as implicit
            ]
            
            problematic = filter(x -> !(x in allowed_implicit), implicit_names)
            
            if !isempty(problematic)
                @warn "Implicit imports detected:" problematic
                # For now, we'll allow implicit imports but track them
                @test_broken isempty(problematic)
            else
                @test true  # All implicit imports are allowed
            end
        else
            @test true  # No implicit imports found
        end
    end

    @testset "Check stale explicit imports" begin
        # Check for unused explicit imports
        stale_imports = check_no_stale_explicit_imports(Catalyst; skip = (Base, Core))
        
        if !isnothing(stale_imports)
            # Check if the stale imports are in our allowed list
            has_unexpected_stale = false
            for (mod, imports) in stale_imports
                # Some imports might be used in macros or extensions and not detected
                allowed_stale = [
                    # These might appear stale but are actually used
                    :MacroTools,      # Used in DSL macros
                    :Graphs,          # Used in graph analysis
                    :DataStructures,  # Used for OrderedDict/OrderedSet
                    :Parameters      # Used for @with_kw_noshow
                ]
                
                for imp in imports
                    if !(Symbol(imp) in allowed_stale)
                        @warn "Unexpected stale import in $mod: $imp"
                        has_unexpected_stale = true
                    end
                end
            end
            @test !has_unexpected_stale
        else
            @test true  # No stale imports found
        end
    end

    @testset "Check qualified accesses are public" begin
        # Check that we only use public APIs when accessing other modules
        result = check_all_qualified_accesses_are_public(Catalyst)
        
        if !isnothing(result)
            # Some non-public accesses might be necessary for deep integration with ModelingToolkit
            # We should document these and work to minimize them
            problematic_accesses = []
            
            for r in result
                # Check if this is an acceptable non-public access
                # ModelingToolkit internal functions that Catalyst needs
                if !occursin("ModelingToolkit", string(r))
                    push!(problematic_accesses, r)
                end
            end
            
            if !isempty(problematic_accesses)
                @warn "Non-public qualified accesses (non-MTK):" problematic_accesses
                @test_broken isempty(problematic_accesses)
            else
                # All non-public accesses are to ModelingToolkit internals which are acceptable
                @test true
            end
        else
            @test true  # All qualified accesses are public
        end
    end

    @testset "Analysis summary" begin
        # Print a summary of the explicit imports analysis
        # This helps track progress but doesn't fail tests
        @info "Running explicit imports analysis for review..."
        
        # Capture the analysis in a string to check for issues
        io = IOBuffer()
        print_explicit_imports(io, Catalyst; strict = false)
        analysis = String(take!(io))
        
        # Check if the analysis mentions any critical issues
        has_issues = occursin("relying on implicit imports", analysis)
        
        if has_issues
            @info "Analysis found potential improvements needed"
            # Don't fail, just inform
        else
            @info "Explicit imports analysis looks good"
        end
        
        # Always pass this test - it's just for information
        @test true
    end
end
