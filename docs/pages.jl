pages = Any[
    "Home" => "home.md",
    "Introduction to Catalyst" => Any[
        "introduction_to_catalyst/catalyst_for_new_julia_users.md",
        "introduction_to_catalyst/introduction_to_catalyst.md"
        # Advanced introduction.
    ],
    "Model Creation and Properties" => Any[
        "model_creation/dsl_description.md",
        # DSL - Advanced
        "model_creation/programmatic_CRN_construction.md",
        "model_creation/coupled_crn_equations.md",
        "model_creation/compositional_modeling.md",
        # Events.
        # Distributed parameters, rates, and initial conditions.
        # Loading and writing models to files.
        # Model visualisation.
        "model_creation/network_analysis.md",
        "model_creation/chemistry_related_functionality.md",
        "Model creation examples" => Any[
            "model_creation/examples/basic_CRN_examples.md",
            "model_creation/examples/hodgkin_huxley_equation.md",
            "model_creation/examples/smoluchowski_coagulation_equation.md"
        ]
    ],
    "Model simulation" => Any[
        # Simulation introduction.
        # Simulation plotting.
        "model_simulation/simulation_structure_interfacing.md",
        # Monte Carlo/Ensemble simulations.
        # Stochastic simulation statistical analysis.
        # ODE Performance considerations/advice.
        # SDE Performance considerations/advice.
        # Jump Performance considerations/advice.
        # Finite state projection
    ],
    "Steady state analysis" => Any[
        "steady_state_functionality/homotopy_continuation.md",
        "steady_state_functionality/nonlinear_solve.md",
        # Stability analysis.
        "steady_state_functionality/bifurcation_diagrams.md"
        # Dynamical systems analysis.
    ],
    "Inverse Problems" => Any[
        # Inverse problems introduction.
        "inverse_problems/optimization_ode_param_fitting.md",
        # "inverse_problems/petab_ode_param_fitting.md",
        # ODE parameter fitting using Turing.
        # SDE/Jump fitting.
        # Non-parameter fitting optimisation.
        "inverse_problems/07_structural_identifiability.md",
        # Practical identifiability.
        # GLobal and local sensitivity analysis.
        "Inverse problem examples" => Any[
            "inverse_problems/examples/ode_fitting_oscillation.md"
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