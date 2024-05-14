pages = Any[
    "Home" => "home.md",
    "Introduction to Catalyst" => Any[
        "01_introduction_to_catalyst/01_catalyst_for_new_julia_users.md",
        "01_introduction_to_catalyst/02_introduction_to_catalyst.md"
        # Advanced introduction.
    ],
    "Model Creation and Properties" => Any[
        "02_model_creation/01_dsl_description.md",
        # DSL - Advanced
        "02_model_creation/03_programmatic_CRN_construction.md",
        "02_model_creation/04_constraint_equations.md",
        "02_model_creation/05_compositional_modeling.md",
        # Events.
        # Distributed parameters, rates, and initial conditions.
        # Loading and writing models to files.
        # Model visualisation.
        "02_reaction_network_functionality/10_network_analysis.md",
        "02_reaction_network_functionality/11_chemistry_related_functionality.md",
        "Model creation examples" => Any[
            "02_model_creation/examples/basic_CRN_examples.md",
            "02_model_creation/examples/hodgkin_huxley_equation.md",
            "02_model_creation/examples/smoluchowski_coagulation_equation.md"
        ]
    ],
    "Model simulation" => Any[
        # Simulation introduction.
        # Simulation plotting.
        "03_model_simulation/03_simulation_structure_interfacing.md",
        # Monte Carlo/Ensemble simulations.
        # Stochastic simulation statistical analysis.
        # ODE Performance considerations/advice.
        # SDE Performance considerations/advice.
        # Jump Performance considerations/advice.
        # Finite state projection
    ],
    "Steady state analysis" => Any[
        "04_steady_state_functionality/01_homotopy_continuation.md",
        "04_steady_state_functionality/02_nonlinear_solve.md",
        # Stability analysis.
        "04_steady_state_functionality/04_bifurcation_diagrams.md"
        # Dynamical systems analysis.
    ],
    "Inverse Problems" => Any[
        # Inverse problems introduction.
        "05_inverse_problems/02_optimization_ode_param_fitting.md",
        # "inverse_problems/03_petab_ode_param_fitting.md",
        # ODE parameter fitting using Turing.
        # SDE/Jump fitting.
        # Non-parameter fitting optimisation.
        "05_inverse_problems/07_structural_identifiability.md",
        # Practical identifiability.
        # GLobal and local sensitivity analysis.
        "Inverse problem examples" => Any[
            "05_inverse_problems/examples/ode_fitting_oscillation.md"
        ]
    ],
    "Spatial modelling" => Any[
        # Intro.
        # Lattice ODEs.
        # Lattice Jumps.
    ],
    # "Developer Documentation" => Any[
    #     # Contributor's guide.
    #     # Repository structure.
    # ],
    "FAQs" => "faqs.md",
    "API" => "api.md"
]