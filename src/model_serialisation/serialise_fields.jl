### Handles Independent Variables ###

# Checks if the reaction system have any independent variable. True for all valid reaction systems.
function has_iv(rn::ReactionSystem)
    return true
end

# Extract a string which declares the system's independent variable.
function get_iv_string(rn::ReactionSystem)
    iv_dec = ModelingToolkit.get_iv(rn)
    return "@variables $(iv_dec)"
end

# Creates an annotation for the system's independent variable.
function get_iv_annotation(rn::ReactionSystem)
    return "Independent variable:"
end

# Combines the 3 independent variable-related functions in a constant tuple.
IV_FS = (has_iv, get_iv_string, get_iv_annotation)


### Handles Spatial Independent Variables ###

# Checks if the reaction system have any spatial independent variables.
function has_sivs(rn::ReactionSystem)
    return !isempty(get_sivs(rn))
end

# Extract a string which declares the system's spatial independent variables.
function get_sivs_string(rn::ReactionSystem)
    return "spatial_ivs = @variables$(get_species_string(get_sivs(rn)))"
end

# Creates an annotation for the system's spatial independent variables.
function get_sivs_annotation(rn::ReactionSystem)
    return "Spatial independent variables:"
end

# Combines the 3 independent variables-related functions in a constant tuple.
SIVS_FS = (has_sivs, get_sivs_string, get_sivs_annotation)


### Handles Species ###

# Checks if the reaction system have any species.
function has_species(rn::ReactionSystem)
    return !isempty(get_species(rn))
end

# Extract a string which declares the system's species.
function get_species_string(rn::ReactionSystem)
    return "sps = @species$(syms_2_declaration_string(get_species(rn)))"
end

# Creates an annotation for the system's species.
function get_species_annotation(rn::ReactionSystem)
    return "Species:"
end

# Combines the 3 species-related functions in a constant tuple.
SPECIES_FS = (has_species, get_species_string, get_species_annotation)


### Handles Variables ###

# Checks if the reaction system have any variables.
function has_variables(rn::ReactionSystem)
    return length(get_unknowns(rn)) > length(get_species(rn))
end

# Extract a string which declares the system's variables.
function get_variables_string(rn::ReactionSystem)
    variables = filter(!isspecies, get_unknowns(rn))
    return "vars = @variables$(syms_2_declaration_string(variables))"
end

# Creates an annotation for the system's .
function get_variables_annotation(rn::ReactionSystem)
    return "Variables:"
end

# Combines the 3 variables-related functions in a constant tuple.
VARIABLES_FS = (has_variables, get_variables_string, get_variables_annotation)


### Handles Parameters ###

# Checks if the reaction system have any parameters.
function has_parameters(rn::ReactionSystem)
    return length(get_ps(rn)) != 0
end

# Extract a string which declares the system's parameters.
function get_parameters_string(rn::ReactionSystem)
    return "ps = @parameters$(syms_2_declaration_string(get_ps(rn)))"
end

# Creates an annotation for the system's parameters.
function get_parameters_annotation(rn::ReactionSystem)
    return "Parameters:"
end

# Combines the 3 parameters-related functions in a constant tuple.
PARAMETERS_FS = (has_parameters, get_parameters_string, get_parameters_annotation)


### Handles Reactions ###

# Checks if the reaction system have any reactions.
function has_reactions(rn::ReactionSystem)
    return length(reactions(rn)) != 0
end

# Extract a string which declares the system's reactions.
function get_reactions_string(rn::ReactionSystem)
    # Creates a dictionary for converting symbolics to their call-stripped form (e.g. X(t) to X).
    strip_call_dict = make_strip_call_dict(rn)

    # Handles the case with one reaction separately. Only effect is nicer formatting.
    (length(get_rxs(rn)) == 1) && (return "rxs = [$(reaction_string(rx, strip_call_dict))]")

    # Creates the string corresponding to the code which generates the system's reactions. 
    rxs_string = "rxs = ["
    for rx in get_rxs(rn)
        @string_append! rxs_string "\n\t" * reaction_string(rx, strip_call_dict) ","
    end

    # Updates the string (including removing the last `,`) and returns it.
    return rxs_string[1:end-1] * "\n]"
end

# Creates a string that corresponds to the declaration of a single `Reaction`.
function reaction_string(rx::Reaction, strip_call_dict)
    # Prepares the `Reaction` declaration components.
    rate = expression_2_string(rx.rate; strip_call_dict)
    substrates = isempty(rx.substrates) ? "nothing" : x_2_string(rx.substrates)    
    products = isempty(rx.products) ? "nothing" : x_2_string(rx.products)
    substoich = isempty(rx.substoich) ? "nothing" : x_2_string(rx.substoich)
    prodstoich = isempty(rx.prodstoich) ? "nothing" : x_2_string(rx.prodstoich)

    # Creates the full expression, including adding kwargs (`only_use_rate` and `metadata`).
    rx_string = "Reaction($rate, $(substrates), $(products), $(substoich), $(prodstoich)"
    if rx.only_use_rate
        @string_append! rx_string "; only_use_rate = true"
        isempty(getmetadata_dict(rx)) || (rx_string = rx_string * ", ")
    end
    if !isempty(getmetadata_dict(rx))
        rx.only_use_rate || (@string_append! rx_string "; ")
        @string_append! rx_string "metadata = ["
        for entry in getmetadata_dict(rx)
            metadata_entry = "$(x_2_string(entry)), "
            @string_append! rx_string metadata_entry
        end
        rx_string = rx_string[1:end-2] * "]"
    end

    # Returns the Reaction string.
    return rx_string * ")"
end

# Creates an annotation for the system's reactions.
function get_reactions_annotation(rn::ReactionSystem)
    return "Reactions:"
end

# Combines the 3 reactions-related functions in a constant tuple.
REACTIONS_FS = (has_reactions, get_reactions_string, get_reactions_annotation)


### Handles Equations ###

# Checks if the reaction system have any equations.
function has_equations(rn::ReactionSystem)
    return length(get_eqs(rn)) > length(get_rxs(rn))
end

# Extract a string which declares the system's equations.
function get_equations_string(rn::ReactionSystem)
    # Creates a dictionary for converting symbolics to their call-stripped form (e.g. X(t) to X).
    strip_call_dict = make_strip_call_dict(rn)

    # Handles the case with one equation separately. Only effect is nicer formatting.
    if length(get_eqs(rn)) - length(get_rxs(rn)) == 1
        return "eqs = [$(expression_2_string(get_eqs(rn)[end]; strip_call_dict))]"
    end

    # Creates the string corresponding to the code which generates the system's reactions. 
    eqs_string = "rxs = ["
    for eq in get_eqs(rn)[length(get_rxs(rn)) + 1:end]
        @string_append! eqs_string "\n\t" expression_2_string(eq; strip_call_dict) ","
    end

    # Updates the string (including removing the last `,`) and returns it.
    return eqs_string[1:end-1] * "\n]"
end

# Creates an annotation for the system's equations.
function get_equations_annotation(rn::ReactionSystem)
    return "Equations:"
end

# Combines the 3 equations-related functions in a constant tuple.
EQUATIONS_FS = (has_equations, get_equations_string, get_equations_annotation)


### Handles Observables ###

# Checks if the reaction system have any observables.
function has_observed(rn::ReactionSystem)
    return false
end

# Extract a string which declares the system's observables.
function get_observed_string(rn::ReactionSystem)
    get_unsupported_comp_string("observables")
end

# Creates an annotation for the system's observables.
function get_observed_annotation(rn::ReactionSystem)
    return "Observables:"
end

# Combines the 3 -related functions in a constant tuple.
OBSERVED_FS = (has_observed, get_observed_string, get_observed_annotation)


### Handles Continuous Events ###

# Checks if the reaction system have any continuous events.
function has_continuous_events(rn::ReactionSystem)
    return length(rn.continuous_events) > 0
end

# Extract a string which declares the system's continuous events.
function get_continuous_events_string(rn::ReactionSystem)
    # Creates a dictionary for converting symbolics to their call-stripped form (e.g. X(t) to X).
    strip_call_dict = make_strip_call_dict(rn)

    # Handles the case with one event separately. Only effect is nicer formatting.
    if length(rn.continuous_events) == 1
        return "continuous_events = [$(continuous_event_string(rn.continuous_events.value[1], strip_call_dict))]"
    end

    # Creates the string corresponding to the code which generates the system's reactions. 
    continuous_events_string = "continuous_events = ["
    for continuous_event in rn.continuous_events.value
        @string_append! continuous_events_string "\n\t" continuous_event_string(continuous_event, strip_call_dict) ","
    end

    # Updates the string (including removing the last `,`) and returns it.
    return continuous_events_string[1:end-1] * "\n]"
end

# Creates a string that corresponds to the declaration of a single continuous event.
function continuous_event_string(continuous_event, strip_call_dict)
    # Creates the string corresponding to the equations (i.e. conditions).
    eqs_string = "["
    for eq in continuous_event.eqs
        @string_append! eqs_string expression_2_string(eq; strip_call_dict) ", "
    end
    eqs_string = eqs_string[1:end-2] * "]"

    # Creates the string corresponding to the affects.
    # Continuous events' `affect` field should probably be called `affects`. Likely the `s` was
    # dropped by mistake in MTK.
    affects_string = "["
    for affect in continuous_event.affect
        @string_append! affects_string expression_2_string(affect; strip_call_dict) ", "
    end
    affects_string = affects_string[1:end-2] * "]"

    return eqs_string * " => " * affects_string
end

# Creates an annotation for the system's continuous events.
function get_continuous_events_annotation(rn::ReactionSystem)
    return "Continuous events:"
end

# Combines the 3 -related functions in a constant tuple.
CONTINUOUS_EVENTS_FS = (has_continuous_events, get_continuous_events_string, get_continuous_events_annotation)


### Handles Discrete Events ###

# Checks if the reaction system have any discrete events.
function has_discrete_events(rn::ReactionSystem)
    return length(rn.discrete_events) > 0
end

# Extract a string which declares the system's discrete events.
function get_discrete_events_string(rn::ReactionSystem)
    # Creates a dictionary for converting symbolics to their call-stripped form (e.g. X(t) to X).
    strip_call_dict = make_strip_call_dict(rn)

    # Handles the case with one event separately. Only effect is nicer formatting.
    if length(rn.discrete_events) == 1
        return "discrete_events = [$(discrete_event_string(rn.discrete_events.value[1], strip_call_dict))]"
    end

    # Creates the string corresponding to the code which generates the system's reactions. 
    discrete_events_string = "discrete_events = ["
    for discrete_event in rn.discrete_events.value
        @string_append! discrete_events_string "\n\t" discrete_event_string(discrete_event, strip_call_dict) ","
    end

    # Updates the string (including removing the last `,`) and returns it.
    return discrete_events_string[1:end-1] * "\n]"
end

# Creates a string that corresponds to the declaration of a single discrete event.
function discrete_event_string(discrete_event, strip_call_dict)
    # Creates the string corresponding to the conditions. The special check is if the condition is
    # an expression like `X > 5.0`. Here, "(...)" is added for purely aesthetic reasons.
    condition_string = x_2_string(discrete_event.condition)
    if discrete_event.condition isa SymbolicUtils.BasicSymbolic
        @string_prepend! "(" condition_string
        @string_append! condition_string ")"
    end

    # Creates the string corresponding to the affects.
    affects_string = "["
    for affect in discrete_event.affects
        @string_append! affects_string expression_2_string(affect; strip_call_dict) ", "
    end
    affects_string = affects_string[1:end-2] * "]"

    return condition_string * " => " * affects_string
end

# Creates an annotation for the system's discrete events.
function get_discrete_events_annotation(rn::ReactionSystem)
    return "Discrete events:"
end

# Combines the 3 -related functions in a constant tuple.
DISCRETE_EVENTS_FS = (has_discrete_events, get_discrete_events_string, get_discrete_events_annotation)


### Handles Systems ###

# Checks if the reaction system have any systems.
function has_systems(rn::ReactionSystem)
    return false
end

# Extract a string which declares the system's systems.
function get_systems_string(rn::ReactionSystem)
    get_unsupported_comp_string("systems")
end

# Creates an annotation for the system's systems.
function get_systems_annotation(rn::ReactionSystem)
    get_unsupported_comp_annotation("Systems:")
end

# Combines the 3 systems-related functions in a constant tuple.
SYSTEMS_FS = (has_systems, get_systems_string, get_systems_annotation)


### Handles Connection Types ###

# Checks if the reaction system have any connection types.
function has_connection_type(rn::ReactionSystem)
    return false
end

# Extract a string which declares the system's connection types.
function get_connection_type_string(rn::ReactionSystem)
    get_unsupported_comp_string("connection types")
end

# Creates an annotation for the system's connection types.
function get_connection_type_annotation(rn::ReactionSystem)
    get_unsupported_comp_annotation("Connection types:")
end

# Combines the 3 connection types-related functions in a constant tuple.
CONNECTION_TYPE_FS = (has_connection_type, get_connection_type_string, get_connection_type_annotation)