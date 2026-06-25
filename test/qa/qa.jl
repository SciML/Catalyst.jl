using SciMLTesting, Catalyst, Test

# JET is a SciMLTesting weak dependency: `using JET` registers it and turns the JET
# check on. JET 0.11 crashes (UndefRefError in `collect_callee_reports!`) when run
# under the Julia 1.13 prerelease `Compiler.jl`, so only load it on the Julia versions
# JET supports. Aqua and ExplicitImports still run on every version. Re-enable JET on
# 1.13 once a JET release analyses cleanly there. Tracked in SciML/Catalyst.jl#1496.
@static if VERSION < v"1.13-"
    using JET
end

# ExplicitImports cannot fully analyze `Catalyst` (the `PhysicalScale` EnumX enum
# parses as an unanalyzable submodule), so the per-module checks are told to allow it.
const EI_ALLOW_UNANALYZABLE = (Catalyst, Catalyst.PhysicalScale)

run_qa(
    Catalyst;
    explicit_imports = true,
    # Test-only [extras] in the root Project.toml intentionally carry no [compat]
    # entries (they are pinned via the resolver, not declared bounds).
    aqua_kwargs = (; deps_compat = (; check_extras = false)),
    # Pre-existing Aqua findings tracked in SciML/Catalyst.jl#1496:
    #   - ambiguities:       55 method ambiguities across the API surface
    #   - unbound_args:      GridLattice `Union{Array{Bool,N}, CartesianGridRej{N,T}}`
    #                        methods leave `T` unbound on the Array branch
    #   - undefined_exports: `Variable` is reexported (via `@reexport using
    #                        ModelingToolkitBase`) but no longer defined upstream
    aqua_broken = (:ambiguities, :unbound_args, :undefined_exports),
    ei_kwargs = (;
        no_implicit_imports = (; allow_unanalyzable = EI_ALLOW_UNANALYZABLE),
        no_stale_explicit_imports = (; allow_unanalyzable = EI_ALLOW_UNANALYZABLE),
    ),
    # Pre-existing ExplicitImports public-API findings tracked in
    # SciML/Catalyst.jl#1496: Catalyst pervasively imports/accesses internal (not
    # `public`/exported) names of ModelingToolkitBase, Symbolics, SymbolicUtils,
    # SciMLBase, DiffEqBase, JumpProcesses, Latexify, Graphs and Base. These go green
    # as the upstream packages mark the names public; tracked-broken until then.
    ei_broken = (
        :all_qualified_accesses_via_owners,
        :all_qualified_accesses_are_public,
        :all_explicit_imports_are_public,
    ),
)
