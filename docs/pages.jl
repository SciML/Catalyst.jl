using Catalyst, ModelingToolkit

pages = Any[
        "Home" => "index.md",
        "Tutorials" => Any[
            "tutorials/using_catalyst.md",
            "tutorials/dsl.md",
            "tutorials/reaction_systems.md",
            "tutorials/basic_examples.md",
            "tutorials/compositional_modeling.md",
            "tutorials/symbolic_stoich.md",
            "tutorials/bifurcation_diagram.md",
            "tutorials/parameter_estimation.md",
            "tutorials/generating_reactions_programmatically.md"
        ],
        "FAQs" => "faqs.md",
        "API" => Any[
            "api/catalyst_api.md"
        ]
    ]
