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

# Legacy Catalyst releases expose this compatibility facade through
# `@reexport using ModelingToolkitBase` and direct imports from its owner packages.
# Keep the exact approved surface version-controlled: changing it is a public-API
# decision, not a consequence of a dependency adding a new public binding.
const LEGACY_DEPENDENCY_REEXPORTS = (
    Symbol("@acrule"), Symbol("@arrayop"), Symbol("@brownian"), Symbol("@brownians"), Symbol("@component"), Symbol("@connector"), Symbol("@constants"), Symbol("@derivative_rule"), Symbol("@derivatives"), Symbol("@discretes"), Symbol("@independent_variables"), Symbol("@makearray"),
    Symbol("@mtkbuild"), Symbol("@mtkcompile"), Symbol("@mtkcomplete"), Symbol("@named"), Symbol("@namespace"), Symbol("@nonamespace"), Symbol("@pack!"), Symbol("@parameters"), Symbol("@poissonians"), Symbol("@register_array_symbolic"), Symbol("@register_derivative"), Symbol("@register_discontinuity"),
    Symbol("@register_inverse"), Symbol("@register_symbolic"), Symbol("@rule"), Symbol("@symbolic_wrap"), Symbol("@syms"), Symbol("@symstruct"), Symbol("@unpack"), Symbol("@variables"), Symbol("@wrapped"), :AbstractCollocation, :AbstractDynamicOptProblem, :AbstractNonlinearProblem,
    :AnalysisPoint, :AssignmentAffect, :BS, :BipartiteGraph, :CartesianGrid, :CartesianGridRej, :CasADiCollocation, :CasADiDynamicOptProblem, :Clock, :CompilerOptions, :Connection, :Differential,
    :DiscreteFunction, :DiscreteProblem, :DiscreteSystem, :DynamicOptSolution, :Equation, :EvalAt, :Flow, :Girsanov_transform, :GlobalScope, :Hold, :HomotopyContinuationProblem, :IRStructure,
    :ImplicitDiscreteFunction, :ImplicitDiscreteProblem, :ImplicitDiscreteSystem, :Inequality, :InfiniteOptCollocation, :InfiniteOptDynamicOptProblem, :Initial, :InitializationProblem, :Integral, :IntervalNonlinearFunction, :IntervalNonlinearProblem, :JuMPCollocation,
    :JuMPDynamicOptProblem, :JumpProblem, :JumpSystem, :LocalScope, :MTKParameters, :MiscSystemData, :MissingGuessValue, :ModelingToolkitBase, :NonlinearFunction, :NonlinearProblem, :NonlinearSystem, :Num,
    :ODEFunction, :ODEProblem, :ODESystem, :OptimizationProblem, :OptimizationSystem, :PDESystem, :ParentScope, :Pre, :PyomoCollocation, :PyomoDynamicOptProblem, :Rewriters, :RuleSet,
    :SDEFunction, :SDEProblem, :SDESystem, :SafeReal, :Sample, :SampleTime, :Shift, :ShiftIndex, :SolverStepClock, :SteadyStateProblem, :Stream, :SymReal,
    :SymScope, :SymStruct, :SymbolicLinearODE, :SymbolicMassActionJump, :SymbolicUtils, :Symbolics, :SymbolicsSparsityDetector, :System, :Term, :TimeDomain, :TreeReal, :UnPack,
    :add_accumulations, :alg_equations, :analytically_integrated, :approximation_function, :arguments, :asdigraph, :asgraph, :bindings, :bound_parameters, :brownians, :build_explicit_observed_function, :build_function, :calculate_control_jacobian,
    :calculate_cost_gradient, :calculate_cost_hessian, :calculate_hessian, :calculate_jacobian, :calculate_massmatrix, :calculate_tgrad, :change_independent_variable, :change_of_variables, :complete, :compose, :connect, :constraints,
    :continuous_events, :convert_system_indepvar, :cost, :debug_system, :diff_equations, :discrete_events, :domain_connect, :eqeq_dependencies, :equation_dependencies, :equations, :expand, :expand_connections,
    :expand_derivatives, :extend, :factors, :flatten, :flatten_fractions, :fractional_to_ordinary, :full_equations, :gather_factor, :generate_W, :generate_control_jacobian, :generate_cost, :generate_cost_gradient,
    :generate_control_function, :generate_cost_hessian, :generate_custom_function, :generate_diffusion_function, :generate_initializesystem, :generate_jacobian, :generate_rhs, :generate_trajectory, :generate_tgrad, :get_alg_eqs, :get_canonical_expr, :get_diff_eqs, :get_reachability, :get_variables,
    :getbounds, :getconnect, :getdist, :getguess, :getmetadata, :getnominal, :getunit, :groebner_basis, :guesses, :has_alg_eqs, :has_alg_equations, :has_diff_eqs,
    :has_diff_equations, :has_inverse, :has_left_inverse, :has_right_inverse, :hasbounds, :hasconnect, :hasdist, :hasguess, :hasmetadata, :hasnominal, :hasunit, :hierarchy,
    :homotopy, :ifelse_branching, :ifelse_eager, :independent_variable, :independent_variables, :infimum, :initial_conditions, :initialization_equations, :instream, :inverse, :inverse_laplace, :irreducibles, :is_derivative,
    :is_groebner_basis, :iscall, :isdisturbance, :isinitial, :isinput, :isirreducible, :isoutput, :istree, :istunable, :jumps, :left_continuous_function, :left_inverse,
    :laplace, :laplace_solve_ode, :limit, :linear_fractional_to_ordinary, :liouville_transform, :majorization_function, :maybe_zeros, :minorization_function, :modelingtoolkitize, :mtkcompile, :noise_to_brownians, :observables, :observed, :open_loop,
    :operation, :parameters, :parse_expr_to_symbolic, :partial_frac_decomposition, :polynomial_coeffs, :populate_ir!, :print_ir, :quick_cancel, :reorder_dimension_by_tunables, :reorder_dimension_by_tunables!, :respecialize, :right_continuous_function, :right_inverse,
    :rootfunction, :semilinear_form, :semipolynomial_form, :semiquadratic_form, :series, :set_defaults, :setmetadata, :setnominal, :simplify, :simplify_fractions, :solve, :solve_for,
    :solve_linear_ode_system, :solve_symbolic_IVP, :sorted_arguments, :state_priorities, :state_priority, :stochastic_integral_transform, :structural_simplify, :subset_tunables, :substitute, :substitute_in_deriv, :substitute_in_deriv_and_depvar, :supremum,
    :symbolic_linear_solve, :symbolic_solve, :symbolic_solve_ode, :symbolics_to_sympy, :symbolics_to_sympy_pythoncall, :sympy_algebraic_solve, :sympy_integrate, :sympy_limit, :sympy_linear_solve, :sympy_ode_solve, :sympy_pythoncall_algebraic_solve, :sympy_pythoncall_integrate,
    :sympy_pythoncall_limit, :sympy_pythoncall_linear_solve, :sympy_pythoncall_ode_solve, :sympy_pythoncall_simplify, :sympy_pythoncall_to_symbolics, :sympy_simplify, :sympy_to_symbolics, :taylor, :taylor_coeff, :term, :terms, :toexpr,
    :toggle_namespacing, :tosymbol, :tunable_parameters, :unknowns, :unwrap_const, :variable_dependencies, :vartype, :varvar_dependencies, :≲, :≳,
)

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
