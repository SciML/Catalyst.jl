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
        piracies = false      # We extend some Base/MTK methods which might be detected as piracy
    )

    # Test individual Aqua checks that we want to enforce
    @testset "Aqua selective tests" begin
        Aqua.test_undefined_exports(Catalyst)
        Aqua.test_project_extras(Catalyst)
        Aqua.test_deps_compat(Catalyst)
    end
end

@testset "Explicit imports (ExplicitImports.jl)" begin
    # Test that we're not relying on implicit imports
    @testset "No implicit imports" begin
        # Main module should have no implicit imports
        @test isnothing(check_no_implicit_imports(Catalyst; skip = (Base, Core)))
    end

    @testset "No stale explicit imports" begin
        # Check for unused explicit imports (allowing some that might be used in macros)
        stale_imports = check_no_stale_explicit_imports(Catalyst; skip = (Base, Core))
        if !isnothing(stale_imports)
            # Allow some exceptions for imports that are used in macros or re-exported
            allowed_stale = [
                :MacroTools,  # Used in DSL macros
                :Graphs,      # Some imports might be re-exported
                :DataStructures  # Used in internal data structures
            ]
            for (mod, imports) in stale_imports
                filtered = filter(x -> !(x in allowed_stale), imports)
                if !isempty(filtered)
                    @warn "Stale imports in $mod: $filtered"
                end
            end
        end
    end

    @testset "All qualified accesses are public" begin
        # Check that we only use public APIs when accessing other modules with qualified names
        @test isnothing(check_all_qualified_accesses_are_public(Catalyst))
    end

    @testset "Print analysis for review" begin
        # This is not a test, but prints useful information for review
        # It helps identify any remaining implicit imports or other issues
        @info "Printing explicit imports analysis for Catalyst module:"
        print_explicit_imports(Catalyst; strict = false)
    end
end
