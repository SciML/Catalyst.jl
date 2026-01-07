### Tests ExplicitImports.jl Compliance ###

# Fetch packages
using Catalyst, ExplicitImports, Test

@testset "ExplicitImports" begin
    # Test that there are no implicit imports (names used without being imported)
    @test check_no_implicit_imports(Catalyst; allow_unanalyzable=(Catalyst, Catalyst.PhysicalScale)) === nothing

    # Test that there are no stale explicit imports (imports that are not used)
    @test check_no_stale_explicit_imports(Catalyst; allow_unanalyzable=(Catalyst, Catalyst.PhysicalScale)) === nothing

    # Test that all explicit imports are from their owner modules (not re-exported)
    @test check_all_explicit_imports_via_owners(Catalyst) === nothing
end
