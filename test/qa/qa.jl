using SciMLTesting, Catalyst, Test

# ExplicitImports only checks an extension module once it exists, and an extension module
# only exists once every one of its trigger weakdeps has been loaded. Without these
# `using`s `Base.get_extension` returns `nothing` for all of them and QA scans no
# extension code at all. Every Catalyst weakdep is loadable here, so all five extensions
# are covered.
using BifurcationKit, CairoMakie, DynamicPolynomials, GraphMakie, HomotopyContinuation,
    Makie, NetworkLayout, StructuralIdentifiability, TermInterface

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

# JET typo-mode surfaces 5 pre-existing latent findings (see the `jet_broken` comment
# below). They are only reported by JET 0.11 (the version resolved on Julia 1.12+); the
# lts lane resolves JET 0.9, which reports nothing, so the JET check must stay a hard
# (un-broken) `test_package` there. `jet_broken` is therefore only enabled on the Julia
# versions whose JET surfaces the findings, otherwise the lts lane would flip to an
# Unexpected Pass. Gating on `VERSION` (not on a JET version query) keeps this a
# compile-time `@static` choice that matches the JET-load gate above.
const JET_BROKEN = VERSION >= v"1.12-"

# Catalyst intentionally re-exports only the dependency-owned names that its
# documentation and tests exercise through `using Catalyst` (see the `export`
# block in src/Catalyst.jl). Keep the exact approved surface version-controlled:
# changing it is a public-API decision, not a consequence of a dependency adding
# a new public binding.
const LEGACY_DEPENDENCY_REEXPORTS = (
    Symbol("@brownians"), Symbol("@discretes"), Symbol("@independent_variables"), Symbol("@mtkcomplete"), Symbol("@named"), Symbol("@nonamespace"),
    Symbol("@parameters"), Symbol("@poissonians"), Symbol("@register_symbolic"), Symbol("@unpack"), Symbol("@variables"), :CartesianGrid,
    :CartesianGridRej, :Differential, :Equation, :GlobalScope, :Initial, :JumpProblem,
    :ModelingToolkitBase, :NonlinearProblem, :ODEProblem, :ParentScope, :Pre, :SDEProblem,
    :SteadyStateProblem, :SymbolicUtils, :Symbolics, :System, :complete, :compose,
    :continuous_events, :equations, :extend, :full_equations, :getmetadata, :hasmetadata,
    :independent_variables, :isinitial, :mtkcompile, :observed, :operation, :parameters, :setmetadata,
    :simplify, :substitute, :unknowns,
)

# ExplicitImports silently skips an extension that fails to load, so assert the
# extension modules actually exist rather than trusting a green `run_qa`.
@testset "Extensions loaded" begin
    for ext in (
            :CatalystBifurcationKitExtension,
            :CatalystCairoMakieExtension,
            :CatalystGraphMakieExtension,
            :CatalystHomotopyContinuationExtension,
            :CatalystStructuralIdentifiabilityExtension,
        )
        @test Base.get_extension(Catalyst, ext) !== nothing
    end
end

run_qa(
    Catalyst;
    reexports_allow = LEGACY_DEPENDENCY_REEXPORTS,
    api_docs_kwargs = (;
        ignore = LEGACY_DEPENDENCY_REEXPORTS,
        rendered_ignore = LEGACY_DEPENDENCY_REEXPORTS,
    ),
    # Test-only [extras] in the root Project.toml intentionally carry no [compat]
    # entries (they are pinned via the resolver, not declared bounds).
    aqua_kwargs = (; deps_compat = (; check_extras = false)),
    # Pre-existing Aqua findings tracked in SciML/Catalyst.jl#1496:
    #   - ambiguities:       55 method ambiguities across the API surface
    #   - unbound_args:      GridLattice `Union{Array{Bool,N}, CartesianGridRej{N,T}}`
    #                        methods leave `T` unbound on the Array branch
    aqua_broken = (:ambiguities, :unbound_args),
    # Pre-existing JET typo-mode findings tracked in SciML/Catalyst.jl#1496. They are
    # latent on master (reproduce byte-identical on the unmodified base) and were never
    # caught before because Catalyst had no JET check. JET surfaces them only on the
    # Julia 1 lane (JET 0.11 on Julia 1.12); the lts lane runs JET 0.9 and does not
    # flag them. `report_package(Catalyst; mode = :typo)` reports 5 names:
    #   - `metadata_entries` (src/dsl.jl:436) and `reaction`
    #     (src/chemistry_functionality.jl:400): wrong variable names interpolated into
    #     `error(...)` strings on already-error paths (should be `metadata_i` / `rx`).
    #   - `clipboard` x2 (src/latexify_recipes.jl:65, 161): `clipboard` (InteractiveUtils)
    #     is not imported; only reachable when `Latexify.COPY_TO_CLIPBOARD` is set.
    #   - `V` (src/reactionsystem.jl:78): the unparameterised `NetworkProperties()`
    #     `@kwdef` keyword constructor references the type parameter `V` in its
    #     `Set{V}()`/`Dict{V,Int}()` defaults; that bare constructor is never called
    #     (only `NetworkProperties{Int, SymbolicT}()` is), so it is a latent dead path.
    # `jet_broken` runs `report_package` and asserts the reports are non-empty, so it
    # auto-flags an Unexpected Pass once the source typos are fixed and the bare
    # constructor is removed/parameterised.
    jet_broken = JET_BROKEN,
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
