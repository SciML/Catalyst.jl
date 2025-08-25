pages = Any[
"Home" => "index.md",
"Introduction to Catalyst" => Any[
    "introduction_to_catalyst/catalyst_for_new_julia_users.md",
    "introduction_to_catalyst/introduction_to_catalyst.md",
    "introduction_to_catalyst/math_models_intro.md"
],
"Model creation and properties" => Any[
    "model_creation/dsl_basics.md",
    "model_creation/dsl_advanced.md",
    "model_creation/programmatic_CRN_construction.md",
    "model_creation/compositional_modeling.md",
    "model_creation/constraint_equations.md",
    "model_creation/conservation_laws.md",
    "model_creation/parametric_stoichiometry.md",
    "model_creation/functional_parameters.md",
    "model_creation/model_file_loading_and_export.md",
    "model_creation/model_visualisation.md",
    "model_creation/reactionsystem_content_accessing.md",
    "model_creation/chemistry_related_functionality.md",
    "Examples" => Any[
        "model_creation/examples/basic_CRN_library.md",
        "model_creation/examples/programmatic_generative_linear_pathway.md",
        "model_creation/examples/hodgkin_huxley_equation.md",
        "model_creation/examples/smoluchowski_coagulation_equation.md",
        "model_creation/examples/noise_modelling_approaches.md"
    ]
],
"Model simulation and visualization" => Any[
    "model_simulation/simulation_introduction.md",
    "model_simulation/simulation_plotting.md",
    "model_simulation/simulation_structure_interfacing.md",
    "model_simulation/ensemble_simulations.md",
    "model_simulation/ode_simulation_performance.md",
    "model_simulation/sde_simulation_performance.md",
    "model_simulation/finite_state_projection_simulation.md",
    "Examples" => Any[
        "model_simulation/examples/periodic_events_simulation.md",
        "model_simulation/examples/activation_time_distribution_measurement.md",
        "model_simulation/examples/interactive_brusselator_simulation.md"
    ]
],
"Network Analysis" => Any[
    "network_analysis/odes.md",
    "network_analysis/crn_theory.md",
    "network_analysis/network_properties.md"
],
"Steady state analysis" => Any[
    "steady_state_functionality/homotopy_continuation.md",
    "steady_state_functionality/nonlinear_solve.md",
    "steady_state_functionality/steady_state_stability_computation.md",
    "steady_state_functionality/bifurcation_diagrams.md",
    "steady_state_functionality/dynamical_systems.md",
    "Examples" => Any[
        "steady_state_functionality/examples/nullcline_plotting.md",
        "steady_state_functionality/examples/bifurcationkit_periodic_orbits.md",
        "steady_state_functionality/examples/bifurcationkit_codim2.md"
    ]
],
"Inverse problems" => Any[
    "inverse_problems/petab_ode_param_fitting.md",
    "inverse_problems/optimization_ode_param_fitting.md",
    "inverse_problems/behaviour_optimisation.md",
    "inverse_problems/structural_identifiability.md",
    "inverse_problems/global_sensitivity_analysis.md"    #"Examples" => Any[    #    "inverse_problems/examples/ode_fitting_oscillation.md"    #]
],
"Spatial modelling" => Any[
    "spatial_modelling/lattice_reaction_systems.md",
    "spatial_modelling/lattice_simulation_structure_ interaction.md",
    "spatial_modelling/lattice_simulation_plotting.md",
    "spatial_modelling/spatial_ode_simulations.md",
    "spatial_modelling/spatial_jump_simulations.md"
],
"FAQs" => "faqs.md",
"API" => Any[
    "api/core_api.md",
    "api/network_analysis_api.md"
],
"Developer Documentation" => "devdocs/dev_guide.md"
]
