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


### Handles Species, Variables, and Parameters ###

# Function which handles the addition of species, variable, and parameter declarations to the file
# text. These must be handled as a unity in case there are default value dependencies between these.
function handle_us_n_ps(file_text::String, rn::ReactionSystem, annotate::Bool)
    # Fetches the systems parameters, species, and variables. Computes the `has_` `Bool`s.
    ps_all = get_ps(rn)
    sps_all = get_species(rn)
    vars_all = filter(!isspecies, get_unknowns(rn))
    has_ps = has_parameters(rn)
    has_sps = has_species(rn)
    has_vars = has_variables(rn)

    # Checks which sets have dependencies which requires managing.
    p_deps = any(depends_on(p, [ps_all; sps_all; vars_all]) for p in ps_all)
    sp_deps = any(depends_on(sp, [sps_all; vars_all]) for sp in sps_all)
    var_deps = any(depends_on(var, vars_all) for var in vars_all)

    # Makes the initial declaration.
    if !p_deps && has_ps
        annotate && (@string_append! file_text "\n\n# " get_parameters_annotation(rn))
        @string_append! file_text "\nps = " get_parameters_string(ps_all)
    end
    if !sp_deps && has_sps
        annotate && (@string_append! file_text "\n\n# " get_species_annotation(rn))
        @string_append! file_text "\nsps = " get_species_string(sps_all)
    end
    if !var_deps && has_vars
        annotate && (@string_append! file_text "\n\n# " get_variables_annotation(rn))
        @string_append! file_text "\nvars = " get_variables_string(vars_all)
    end

    # If any set have dependencies, handle these.
    # There are cases where the dependent syms come after their dependencies in the vector
    # (e.g. corresponding to `@parameters p1 p2=p1`)
    # which would not require this special treatment. However, this is currently not considered.
    # Considering it would make the written code prettier, but would also require additional
    # work in these functions to handle these cases (can be sorted out in the future).
    if p_deps || sp_deps || var_deps
        # Builds an annotation mentioning specially handled stuff.
        if annotate
            @string_append! file_text "\n\n# Some "
            p_deps && (@string_append! file_text "parameters, ")
            sp_deps && (@string_append! file_text "species, ")
            var_deps && (@string_append! file_text "variables, ")
            file_text = file_text[1:end-2]
            @string_append! file_text " depends on the declaration of other parameters, species, and/or variables.\n# These are specially handled here.\n"
        end

        # Pre-declares the sets with written/remaining parameters/species/variables.
        # Whenever all/none are written depends on whether there were any initial dependencies.
        remaining_ps = (p_deps ? ps_all : [])
        remaining_sps = (sp_deps ? sps_all : [])
        remaining_vars = (var_deps ? vars_all : [])

        # Iteratively loops through all parameters, species, and/or variables. In each iteration, 
        # adds the declaration of those that can still be declared.
        while !(isempty(remaining_ps) && isempty(remaining_sps) && isempty(remaining_vars))
            # Checks which parameters/species/variables can be written.
            writable_ps, nonwritable_ps = dependency_split([remaining_ps; remaining_sps; remaining_vars], remaining_ps)
            writable_sps, nonwritable_sps = dependency_split([remaining_sps; remaining_vars], remaining_sps)
            writable_vars, nonwritable_vars = dependency_split(remaining_vars, remaining_vars)
            
            # Writes those that can be written.
            isempty(writable_ps) || @string_append! file_text get_parameters_string(writable_ps) "\n"
            isempty(writable_sps) || @string_append! file_text get_species_string(writable_sps) "\n"
            isempty(writable_vars) || @string_append! file_text get_variables_string(writable_vars) "\n"

            # Updates the remaining parameters/species/variables sets.
            remaining_ps = nonwritable_ps
            remaining_sps = nonwritable_sps
            remaining_vars = nonwritable_vars
        end

        # For parameters, species, and/or variables with dependencies, creates final vectors.
        p_deps && (@string_append! file_text "ps = " syms_2_strings(ps_all) "\n")
        sp_deps && (@string_append! file_text "sps = " syms_2_strings(sps_all) "\n")
        var_deps && (@string_append! file_text "vars = " syms_2_strings(vars_all) "\n")
        file_text = file_text[1:end-1]
    end

    # Returns the finalised output.
    return file_text, has_ps, has_sps, has_vars
end


### Handles Parameters ###
# Unlike most other fields, there are not called via `push_field`, but rather via `handle_us_n_ps`.
# Hence they work slightly differently.

# Checks if the reaction system have any parameters.
function has_parameters(rn::ReactionSystem)
    return !isempty(get_ps(rn))
end

# Extract a string which declares the system's parameters. Uses multiline declaration (a 
# `begin ... end` block) if more than 3 parameters have a "complicated" declaration (if they
# have metadata, default value, or type designation).
function get_parameters_string(ps)
    multiline_format = count(complicated_declaration(p) for p in ps) > 3
    return "@parameters$(syms_2_declaration_string(ps; multiline_format))"
end

# Creates an annotation for the system's parameters.
function get_parameters_annotation(rn::ReactionSystem)
    return "Parameters:"
end


### Handles Species ###
# Unlike most other fields, there are not called via `push_field`, but rather via `handle_us_n_ps`.
# Hence they work slightly differently.

# Checks if the reaction system have any species.
function has_species(rn::ReactionSystem)
    return !isempty(get_species(rn))
end

# Extract a string which declares the system's species. Uses multiline declaration (a 
# `begin ... end` block) if more than 3 species have a "complicated" declaration (if they
# have metadata, default value, or type designation).
function get_species_string(sps)
    multiline_format = count(complicated_declaration(sp) for sp in sps) > 3
    return "@species$(syms_2_declaration_string(sps; multiline_format))"
end

# Creates an annotation for the system's species.
function get_species_annotation(rn::ReactionSystem)
    return "Species:"
end


### Handles Variables ###
# Unlike most other fields, there are not called via `push_field`, but rather via `handle_us_n_ps`.
# Hence they work slightly differently.

# Checks if the reaction system have any variables.
function has_variables(rn::ReactionSystem)
    return length(get_unknowns(rn)) > length(get_species(rn))
end

# Extract a string which declares the system's variables. Uses multiline declaration (a 
# `begin ... end` block) if more than 3 variables have a "complicated" declaration (if they
# have metadata, default value, or type designation).
function get_variables_string(vars)
    multiline_format = count(complicated_declaration(var) for var in vars) > 3
    return "@variables$(syms_2_declaration_string(vars; multiline_format))"
end

# Creates an annotation for the system's .
function get_variables_annotation(rn::ReactionSystem)
    return "Variables:"
end

# Combines the 3 variables-related functions in a constant tuple.
VARIABLES_FS = (has_variables, get_variables_string, get_variables_annotation)


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
    return !isempty(observed(rn))
end

# Extract a string which declares the system's observables.
function get_observed_string(rn::ReactionSystem)
    # Finds the observable species and variables.
    observed_unknowns = [obs_eq.lhs for obs_eq in observed(rn)]
    observed_species = filter(isspecies, observed_unknowns)
    observed_variables = filter(!isspecies, observed_unknowns)

    # Creates a dictionary for converting symbolics to their call-stripped form (e.g. X(t) to X).
    strip_call_dict = make_strip_call_dict([get_unknowns(rn); observed_unknowns])

    # Initialises the observables string with declaring the observable species/variables.
    observed_string = ""
    if !isempty(observed_species)
        @string_append! observed_string "@species$(syms_2_declaration_string(observed_species))\n"
    end
    if !isempty(observed_variables)
        @string_append! observed_string "@variables$(syms_2_declaration_string(observed_variables))\n"
    end

    # Handles the case with one observable separately. Only effect is nicer formatting.
    if length(observed(rn)) == 1
        @string_append! observed_string "observed = [$(expression_2_string(observed(rn)[1]; strip_call_dict))]"
        return observed_string
    end

    # Appends with the code which will generate observables equations.
    @string_append! observed_string "observed = ["
    for obs in observed(rn)
        @string_append! observed_string "\n\t" expression_2_string(obs, strip_call_dict) ","
    end

    # Updates the string (including removing the last `,`) and returns it.
    return observed_string[1:end-1] * "\n]"
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

# Specific `push_field` function, which is used for the system field (where the annotation option
# must be passed to the `get_component_string` function). Since non-ReactionSystem systems cannot be 
# written to file, this functions throws an error if any such systems are encountered.
function push_systems_field(file_text::String, rn::ReactionSystem, annotate::Bool)
    # Checks whther there are any subsystems, and if these are ReactionSystems.
    has_systems(rn) || (return (file_text, false))
    if any(!(system isa ReactionSystem) for system in ModelingToolkit.get_systems(rn)) 
        error("Tries to write a ReactionSystem to file which have non-ReactionSystem subs-systems. This is currently not possible.")
    end

    # Adds the system declaration string to the file string.
    write_string = "\n" * get_systems_string(rn, annotate)
    annotate && (@string_prepend! "\n\n# " get_systems_annotation(rn) write_string)
    return (file_text * write_string, true)
end

# Checks if the reaction system have any systems.
function has_systems(rn::ReactionSystem)
    return !isempty(ModelingToolkit.get_systems(rn))
end

# Extract a string which declares the system's systems.
function get_systems_string(rn::ReactionSystem, annotate::Bool)
    # Initiates the `systems` string. It is pre-declared vector, into which the systems are added.
    systems_string = "systems = Vector(undef, $(length(ModelingToolkit.get_systems(rn))))"

    # Loops through all systems, adding their declaration to the system string.
    for (idx, system) in enumerate(ModelingToolkit.get_systems(rn))
        annotate && (@string_append! systems_string "\n\n# Declares subsystem: $(getname(system))")

        # Manipulates the subsystem declaration to make it nicer.
        subsystem_string = get_full_system_string(system, annotate)
        subsystem_string = replace(subsystem_string, "\n" => "\n\t")
        subsystem_string = "let\n" * subsystem_string[7:end-6] * "end"
        @string_append! systems_string "\nsystems[$idx] = " subsystem_string
    end

    return systems_string
end

# Creates an annotation for the system's systems.
function get_systems_annotation(rn::ReactionSystem)
    return "Subystems:"
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