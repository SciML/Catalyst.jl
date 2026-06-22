# Tests for import hygiene using ExplicitImports.jl
# Ensures Catalyst imports symbols from their owner modules and has no stale imports.
using Catalyst, ExplicitImports, Test

@testset "Explicit Imports" begin
    # Test that there are no implicit imports
    # allow_unanalyzable is needed because EnumX-generated submodules cause parsing issues
    @test check_no_implicit_imports(Catalyst;
        allow_unanalyzable=(Catalyst, Catalyst.PhysicalScale, Catalyst.AliasClass)) === nothing

    # Test that there are no stale explicit imports
    @test check_no_stale_explicit_imports(Catalyst;
        allow_unanalyzable=(Catalyst, Catalyst.PhysicalScale, Catalyst.AliasClass)) === nothing

    # Test that all explicit imports are from owner modules
    @test check_all_explicit_imports_via_owners(Catalyst) === nothing
end
