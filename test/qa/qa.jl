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

# JET typo-mode surfaces 5 pre-existing latent findings (see the `jet_broken` comment
# below). They are only reported by JET 0.11 (the version resolved on Julia 1.12+); the
# lts lane resolves JET 0.9, which reports nothing, so the JET check must stay a hard
# (un-broken) `test_package` there. `jet_broken` is therefore only enabled on the Julia
# versions whose JET surfaces the findings, otherwise the lts lane would flip to an
# Unexpected Pass. Gating on `VERSION` (not on a JET version query) keeps this a
# compile-time `@static` choice that matches the JET-load gate above.
const JET_BROKEN = VERSION >= v"1.12-"

# Catalyst's only blanket reexport is `@reexport using ModelingToolkitBase`
# (src/Catalyst.jl). 312 of the 314 public names on Catalyst's surface that Catalyst does
# not own arrive through it -- from ModelingToolkitBase itself and from everything
# MTKBase re-exports in turn (Symbolics, SymbolicUtils, SciMLBase, TermInterface, UnPack,
# JumpProcesses, BipartiteGraphs, CommonSolve).
#
# Approve that surface by provenance rather than by name. Enumerating it re-breaks QA
# whenever any of those packages adds a public binding, through no change of Catalyst's
# own: Symbolics 7.35 adding the Laplace helpers is exactly that, and is what prompted
# this. Keying on the reexport states the actual intent -- Catalyst deliberately
# re-exposes ModelingToolkitBase's public API -- and still fails if Catalyst starts
# exposing something that reexport does not account for.
const MTKBASE = Catalyst.ModelingToolkitBase

function _from_mtkbase(name)
    Base.ispublic(MTKBASE, name) || return false
    isdefined(MTKBASE, name) || return false
    return try
        getfield(MTKBASE, name) === getfield(Catalyst, name)
    catch
        false
    end
end

const MTKBASE_REEXPORTS = Tuple(filter(_from_mtkbase, public_reexports(Catalyst)))

# The only two public names Catalyst exposes from outside its own module that do *not*
# come from that reexport. Catalyst imports them from JumpProcesses and exports them
# itself, which is a deliberate API decision rather than a side effect, so they stay
# pinned by name.
const DIRECT_REEXPORTS = (:CartesianGrid, :CartesianGridRej)

const APPROVED_REEXPORTS = (MTKBASE_REEXPORTS..., DIRECT_REEXPORTS...)

run_qa(
    Catalyst;
    reexports_allow = APPROVED_REEXPORTS,
    api_docs_kwargs = (;
        ignore = APPROVED_REEXPORTS,
        rendered_ignore = APPROVED_REEXPORTS,
    ),
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
