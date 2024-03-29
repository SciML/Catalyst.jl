"""
    save_reaction_network(filename::String, rn::ReactionSystem; annotate = true)

Save a `ReactionSystem` model to a file.

Work in progress, currently missing features:
- Problems with ordering of declarations of species/variables/parameters that have defaults that are other species/variables/parameters.
- Saving of the `sivs` field has not been fully implemented.
- Saving of the `observed` field has not been fully implemented.
- Saving of the `continuous_events` field has not been fully implemented.
- Saving of the `discrete_events` field has not been fully implemented.
- Saving of the `systems` field has not been fully implemented.
- Saving of the `connection_type` field has not been fully implemented.
"""
function save_reaction_network(filename::String, rn::ReactionSystem; annotate = true)
    # Initiates the file string.
    file_text = ""

    # Goes through each type of system component, potentially adding it to the string.
    file_text, _ = push_field(file_text, rn, annotate, IV_FS)
    file_text, has_sivs = push_field(file_text, rn, annotate, SIVS_FS)
    file_text, has_species = push_field(file_text, rn, annotate, SPECIES_FS)
    file_text, has_variables = push_field(file_text, rn, annotate, VARIABLES_FS)
    file_text, has_parameters = push_field(file_text, rn, annotate, PARAMETERS_FS)
    file_text, has_reactions = push_field(file_text, rn, annotate, REACTIONS_FS)
    file_text, has_equations = push_field(file_text, rn, annotate, EQUATIONS_FS)
    file_text, has_observed = push_field(file_text, rn, annotate, OBSERVED_FS)
    file_text, has_continuous_events = push_field(file_text, rn, annotate, CONTINUOUS_EVENTS_FS)
    file_text, has_discrete_events = push_field(file_text, rn, annotate, DISCRETE_EVENTS_FS)   
    file_text, has_systems = push_field(file_text, rn, annotate, SYSTEMS_FS)
    file_text, has_connection_type = push_field(file_text, rn, annotate, CONNECTION_TYPE_FS)

    # Finalises the system. Creates the final `ReactionSystem` call.
    rs_creation_code = make_reaction_system_call(rn, file_text, annotate,
                                                 has_sivs, has_species, has_variables, has_parameters, 
                                                 has_reactions, has_equations, has_observed, 
                                                 has_discrete_events, has_continuous_events,
                                                 has_systems, has_connection_type)
    annotate || (@string_prepend! "\n" file_text) 
    @string_prepend! "let" file_text 
    @string_append! file_text "\n\n" rs_creation_code "\n\nend"

    # Writes the model to a file. Then, returns nothing.
    open(filename, "w") do file
        write(file, file_text)
    end
    return nothing
end

# Takes the actual text which creates the model, and wraps it in a `let ... end` statement and a
# ReactionSystem call. This creates the finalised text that is written to a file.
function make_reaction_system_call(rs::ReactionSystem, file_text, annotate, has_sivs, has_species, 
                                   has_variables, has_parameters, has_reactions, has_equations, 
                                   has_observed, has_continuous_events, has_discrete_events,
                                   has_systems, has_connection_type)

    # Gets the independent variable input.
    iv = x_2_string(get_iv(rs))

    # Gets the equations (reactions + equations) input.
    if has_reactions && has_equations
        eqs = "[rxs; eqs]"
    elseif has_reactions
        eqs = "rxs"
    elseif has_equations
        eqs = "eqs"
    else 
        eqs = "[]"
    end

    # Gets the unknowns (species + variables) input.
    if has_species && has_variables
        unknowns = "[sps; vars]"
    elseif has_species
        unknowns = "sps"
    elseif has_variables
        unknowns = "vars"
    else 
        unknowns = "[]"
    end

    # Gets the parameters input.
    if has_parameters
        ps = "ps"
    else 
        ps = "[]"
    end

    # Initiates the ReactionSystem call with the mandatory inputs.
    reaction_system_string = "ReactionSystem($eqs, $iv, $unknowns, $ps"

    # Appends the reaction system name. Also initiates the optional argument part of the call.
    if Base.isidentifier(Catalyst.getname(rs))
        rs_name = ":$(Catalyst.getname(rs))"
    else
        rs_name = "Symbol(\"$(Catalyst.getname(rs))\")"
    end
    @string_append! reaction_system_string "; name = $(rs_name)"

    # Goes through various fields that might exists, and if so, adds them to the string.
    has_sivs && (@string_append! reaction_system_string ", spatial_ivs")
    has_observed && (@string_append! reaction_system_string ", observed")
    has_continuous_events && (@string_append! reaction_system_string ", continuous_events")
    has_discrete_events && (@string_append! reaction_system_string ", discrete_events")
    has_systems && (@string_append! reaction_system_string ", systems")
    has_connection_type && (@string_append! reaction_system_string ", connection_type")

    # Potentially appends a combinatorial combinatoric_ratelaws statement.
    Symbolics.unwrap(rs.combinatoric_ratelaws) || (@string_append! reaction_system_string ", combinatoric_ratelaws = false")

    # Potentially appends `ReactionSystem` metadata value(s). Weird composite types are not supported.
    isnothing(rs.metadata) || (@string_append! reaction_system_string ", metadata = $(x_2_string(rs.metadata))")

    # Finalises the call. Appends potential annotation. If the system is complete, add a call for this. 
    @string_append! reaction_system_string ")"
    if !ModelingToolkit.iscomplete(rs)
        @string_append! reaction_system_string "rs = $(reaction_system_string)\ncomplete(rs)"
    end
    if annotate 
        @string_prepend! "# Declares ReactionSystem model:\n" reaction_system_string
    end
    return reaction_system_string
end