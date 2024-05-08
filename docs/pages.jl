pages = Any[
    "Home" => "home.md",
    "Introduction to Catalyst" => Any[
        "01_introduction_to_catalyst/01_catalyst_for_new_julia_users.md",
        "01_introduction_to_catalyst/02_introduction_to_catalyst.md"
        # Advanced introduction.
    ],
    "Model Creation" => Any[
        "02_model_creation/01_dsl_description.md",
        # DSL - Advanced
        "02_model_creation/03_programmatic_CRN_construction.md",
        "02_model_creation/04_constraint_equations.md",
        # Events.
        "02_model_creation/06_compositional_modeling.md",
        # Distributed parameters, rates, and initial conditions.
        # Loading and writing models to files.
        "Model creation examples" => Any[
            "02_model_creation/examples/basic_CRN_examples.md",
            "02_model_creation/examples/hodgkin_huxley_equation.md",
            "02_model_creation/examples/smoluchowski_coagulation_equation.md"
        ]
    ],
    "Reaction network-related functionalities" => Any[
        # Model visualisation.
        "03_reaction_network_functionality/network_analysis.md",
        "03_reaction_network_functionality/chemistry_related_functionality.md",
    ],
    "Model simulation" => Any[
        # Simulation introduction.
        # Simulation plotting.
        "04_model_simulation/03_simulation_structure_interfacing.md",
        # Monte Carlo/Ensemble simulations.
        # Stochastic simulation statistical analysis.
        # ODE Performance considerations/advice.
        # SDE Performance considerations/advice.
        # Jump Performance considerations/advice.
        # Finite state projection
    ],
    "Spatial modelling" => Any[
        # Intro.
        # Lattice ODEs.
        # Lattice Jumps.
    ],
    "Steady state analysis" => Any[
        "06_steady_state_functionality/01_homotopy_continuation.md",
        "06_steady_state_functionality/02_nonlinear_solve.md",
        # Stability analysis.
        "06_steady_state_functionality/04_bifurcation_diagrams.md"
        # Dynamical systems analysis.
    ],
    "Inverse Problems" => Any[
        # Inverse problems introduction.
        "07_inverse_problems/02_optimization_ode_param_fitting.md",
        # "inverse_problems/03_petab_ode_param_fitting.md",
        # ODE parameter fitting using Turing.
        # SDE/Jump fitting.
        # Non-parameter fitting optimisation.
        "07_inverse_problems/07_structural_identifiability.md",
        # Practical identifiability.
        # GLobal and local sensitivity analysis.
        "Inverse problem examples" => Any[
            "07_inverse_problems/examples/ode_fitting_oscillation.md"
        ]
    ],
    # "Developer Documentation" => Any[
    #     # Contributor's guide.
    #     # Repository structure.
    # ],
    "FAQs" => "faqs.md",
    "API" => "api.md"
]