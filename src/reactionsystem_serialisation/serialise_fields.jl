### Handles Independent Variables ###

# Checks if the reaction system has any independent variable. True for all valid reaction systems.
function seri_has_iv(rn::ReactionSystem)
    return true
end

# Extract a string which declares the system's independent variable.
function get_iv_string(rn::ReactionSystem)
    iv_dec = MT.get_iv(rn)
    return "@independent_variables $(iv_dec)"
end

# Creates an annotation for the system's independent variable.
function get_iv_annotation(rn::ReactionSystem)
    return "Independent variable:"
end

# Combines the 3 independent variable-related functions in a constant tuple.
IV_FS = (seri_has_iv, get_iv_string, get_iv_annotation)

### Handles Spatial Independent Variables ###

# Checks if the reaction system has any spatial independent variables.
function seri_has_sivs(rn::ReactionSystem)
    return !isempty(get_sivs(rn))
end

# Extract a string which declares the system's spatial independent variables.
function get_sivs_string(rn::ReactionSystem)
    return "spatial_ivs = @variables$(syms_2_declaration_string(get_sivs(rn)))"
end

# Creates an annotation for the system's spatial independent variables.
function get_sivs_annotation(rn::ReactionSystem)
    return "Spatial independent variables:"
end

# Combines the 3 independent variables-related functions in a constant tuple.
SIVS_FS = (seri_has_sivs, get_sivs_string, get_sivs_annotation)

### Handles Species, Variables, and Parameters ###

# Function which handles the addition of species, variable, and parameter declarations to the file
# text. These must be handled as a unity in case there are default value dependencies between these.
function handle_us_n_ps(
        file_text::String, rn::ReactionSystem, annotate::Bool,
        top_level::Bool
    )
    # Fetches the system's parameters, species, and variables. Computes the `has_` `Bool`s.
    ps_all = filter(!_is_discrete, get_ps(rn))
    discs_all = filter(_is_discrete, get_ps(rn))
    sps_all = get_species(rn)
    vars_all = filter(!isspecies, get_unknowns(rn))
    has_ps = seri_has_parameters(rn)
    has_discs = seri_has_discretes(rn)
    has_sps = seri_has_species(rn)
    has_vars = seri_has_variables(rn)

    # Checks which sets have dependencies which require managing.
    p_deps = any(depends_on(p, [ps_all; discs_all; sps_all; vars_all]) for p in ps_all)
    disc_deps = any(depends_on(d, [discs_all; sps_all; vars_all]) for d in discs_all)
    sp_deps = any(depends_on(sp, [sps_all; vars_all]) for sp in sps_all)
    var_deps = any(depends_on(var, vars_all) for var in vars_all)

    # Makes the initial declaration.
    us_n_ps_string = ""
    if !p_deps && has_ps
        annotate && (@string_append! us_n_ps_string "\n\n# " get_parameters_annotation(rn))
        @string_append! us_n_ps_string "\nps = " get_parameters_string(ps_all)
    end
    if !disc_deps && has_discs
        annotate && (@string_append! us_n_ps_string "\n\n# " get_discretes_annotation(rn))
        @string_append! us_n_ps_string "\ndiscs = " get_discretes_string(discs_all)
    end
    if !sp_deps && has_sps
        annotate && (@string_append! us_n_ps_string "\n\n# " get_species_annotation(rn))
        @string_append! us_n_ps_string "\nsps = " get_species_string(sps_all)
    end
    if !var_deps && has_vars
        annotate && (@string_append! us_n_ps_string "\n\n# " get_variables_annotation(rn))
        @string_append! us_n_ps_string "\nvars = " get_variables_string(vars_all)
    end

    # If any set have dependencies, handle these.
    # There are cases where the dependent syms come after their dependencies in the vector
    # (e.g. corresponding to `@parameters p1 p2=p1`)
    # which would not require this special treatment. However, this is currently not considered.
    # Considering it would make the written code prettier, but would also require additional
    # work in these functions to handle these cases (can be sorted out in the future).
    if p_deps || disc_deps || sp_deps || var_deps
        # Builds an annotation mentioning specially handled stuff.
        if annotate
            @string_append! us_n_ps_string "\n\n# Some "
            p_deps && (@string_append! us_n_ps_string "parameters, ")
            disc_deps && (@string_append! us_n_ps_string "discretes, ")
            sp_deps && (@string_append! us_n_ps_string "species, ")
            var_deps && (@string_append! us_n_ps_string "variables, ")
            us_n_ps_string = get_substring_end(us_n_ps_string, 1, -2)
            @string_append! us_n_ps_string " depends on the declaration of other parameters, species, and/or variables.\n# These are specially handled here.\n"
        end

        # Pre-declares the sets with written/remaining parameters/species/variables.
        # Whenever all/none are written depends on whether there were any initial dependencies.
        # `deepcopy` is required as these get mutated by `dependency_split!`.
        remaining_ps = (p_deps ? deepcopy(ps_all) : [])
        remaining_discs = (disc_deps ? deepcopy(discs_all) : [])
        remaining_sps = (sp_deps ? deepcopy(sps_all) : [])
        remaining_vars = (var_deps ? deepcopy(vars_all) : [])

        # Iteratively loops through all parameters, species, and/or variables. In each iteration,
        # adds the declaration of those that can still be declared.
        while !(isempty(remaining_ps) && isempty(remaining_discs) && isempty(remaining_sps) && isempty(remaining_vars))
            # Checks which parameters/species/variables can be written. The `dependency_split`
            # function updates the `remaining_` input.
            writable_ps = dependency_split!(remaining_ps, [remaining_ps; remaining_discs; remaining_sps; remaining_vars])
            writable_discs = dependency_split!(remaining_discs, [remaining_ps; remaining_discs; remaining_sps; remaining_vars])
            writable_sps = dependency_split!(remaining_sps, [remaining_ps; remaining_discs; remaining_sps; remaining_vars])
            writable_vars = dependency_split!(remaining_vars, [remaining_ps; remaining_discs; remaining_sps; remaining_vars])

            # Writes those that can be written.
            isempty(writable_ps) ||
                @string_append! us_n_ps_string get_parameters_string(writable_ps) "\n"
            isempty(writable_discs) ||
                @string_append! us_n_ps_string get_discretes_string(writable_discs) "\n"
            isempty(writable_sps) ||
                @string_append! us_n_ps_string get_species_string(writable_sps) "\n"
            isempty(writable_vars) ||
                @string_append! us_n_ps_string get_variables_string(writable_vars) "\n"
        end

        # For parameters, species, and/or variables with dependencies, create final vectors.
        p_deps && (@string_append! us_n_ps_string "ps = " syms_2_strings(ps_all) "\n")
        disc_deps && (@string_append! us_n_ps_string "discs = " syms_2_strings(discs_all) "\n")
        sp_deps && (@string_append! us_n_ps_string "sps = " syms_2_strings(sps_all) "\n")
        var_deps && (@string_append! us_n_ps_string "vars = " syms_2_strings(vars_all) "\n")
        us_n_ps_string = get_substring_end(us_n_ps_string, 1, -1)
    end

    # If this is not a top-level system, `local ` must be added to all declarations.
    if !top_level
        us_n_ps_string = replace(us_n_ps_string, "\nps = " => "\nlocal ps = ")
        us_n_ps_string = replace(us_n_ps_string, "\ndiscs = " => "\nlocal discs = ")
        us_n_ps_string = replace(us_n_ps_string, "\nsps = " => "\nlocal sps = ")
        us_n_ps_string = replace(us_n_ps_string, "\nvars = " => "\nlocal vars = ")
    end

    # Merges the file text with `us_n_ps_string` and returns the final outputs.
    return file_text * us_n_ps_string, has_ps, has_discs, has_sps, has_vars
end

### Handles Parameters ###
# Unlike most other fields, these are not called via `push_field`, but rather via `handle_us_n_ps`.
# Hence they work slightly differently.

# Checks if the reaction system has any parameters.
function seri_has_parameters(rn::ReactionSystem)
    return any(!_is_discrete(p) for p in get_ps(rn))
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

### Handles Discretes ###
# Unlike most other fields, these are not called via `push_field`, but rather via `handle_us_n_ps`.
# Hence they work slightly differently.
# Discretes are basically time-dependent parameters.

# Checks if the reaction system has any discretes.
function seri_has_discretes(rn::ReactionSystem)
    return any(_is_discrete(p) for p in get_ps(rn))
end
# Temporary function, replace with better after reply from Aayush.
_is_discrete(param) = occursin('(', "$(param)")

# Extract a string which declares the system's discretes. Uses multiline declaration (a
# `begin ... end` block) if more than 3 discretes have a "complicated" declaration (if they
# have metadata, default value, or type designation).
function get_discretes_string(discs)
    multiline_format = count(complicated_declaration(d) for d in discs) > 3
    return "@discretes$(syms_2_declaration_string(discs; multiline_format))"
end

# Creates an annotation for the system's discretes.
function get_discretes_annotation(rn::ReactionSystem)
    return "Discretes:"
end

### Handles Species ###
# Unlike most other fields, these are not called via `push_field`, but rather via `handle_us_n_ps`.
# Hence they work slightly differently.

# Checks if the reaction system has any species.
function seri_has_species(rn::ReactionSystem)
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
# Unlike most other fields, these are not called via `push_field`, but rather via `handle_us_n_ps`.
# Hence they work slightly differently.

# Checks if the reaction system has any variables.
function seri_has_variables(rn::ReactionSystem)
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
VARIABLES_FS = (seri_has_variables, get_variables_string, get_variables_annotation)

### Handles Reactions ###

# Checks if the reaction system has any reactions.
function seri_has_reactions(rn::ReactionSystem)
    return length(reactions(rn)) != 0
end

# Extract a string which declares the system's reactions.
function get_reactions_string(rn::ReactionSystem)
    # Creates a dictionary for converting symbolics to their call-stripped form (e.g. X(t) to X).
    strip_call_dict = make_strip_call_dict(rn)

    # Handles the case with one reaction separately. Only effect is nicer formatting.
    if length(get_rxs(rn)) == 1
        return "rxs = [$(reaction_string(get_rxs(rn)[1], strip_call_dict))]"
    end

    # Creates the string corresponding to the code which generates the system's reactions.
    rxs_string = "rxs = ["
    for rx in get_rxs(rn)
        @string_append! rxs_string "\n\t" * reaction_string(rx, strip_call_dict) ","
    end

    # Updates the string (including removing the last `,`) and returns it.
    return get_substring_end(rxs_string, 1, -1) * "\n]"
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
        rx_string = get_substring_end(rx_string, 1, -2) * "]"
    end

    # Returns the Reaction string.
    return rx_string * ")"
end

# Creates an annotation for the system's reactions.
function get_reactions_annotation(rn::ReactionSystem)
    return "Reactions:"
end

# Combines the 3 reaction-related functions in a constant tuple.
REACTIONS_FS = (seri_has_reactions, get_reactions_string, get_reactions_annotation)

### Handles Equations ###

# Checks if the reaction system has any equations.
function seri_has_equations(rn::ReactionSystem)
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
    eqs_string = "eqs = ["
    for eq in get_eqs(rn)[(length(get_rxs(rn)) + 1):end]
        @string_append! eqs_string "\n\t" expression_2_string(eq; strip_call_dict) ","
    end

    # Updates the string (including removing the last `,`) and returns it.
    return get_substring_end(eqs_string, 1, -1) * "\n]"
end

# Creates an annotation for the system's equations.
function get_equations_annotation(rn::ReactionSystem)
    return "Equations:"
end

# Combines the 3 equations-related functions in a constant tuple.
EQUATIONS_FS = (seri_has_equations, get_equations_string, get_equations_annotation)

### Handles Observables ###

# Checks if the reaction system has any observables.
function seri_has_observed(rn::ReactionSystem)
    return !isempty(get_observed(rn))
end

# Extract a string which declares the system's observables.
function get_observed_string(rn::ReactionSystem)
    # Finds the observable species and variables.
    observed_unknowns = [obs_eq.lhs for obs_eq in MT.get_observed(rn)]
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
    if length(MT.get_observed(rn)) == 1
        @string_append! observed_string "observed = [$(expression_2_string(MT.get_observed(rn)[1]; strip_call_dict))]"
        return observed_string
    end

    # Appends with the code which will generate observables equations.
    @string_append! observed_string "observed = ["
    for obs in MT.get_observed(rn)
        @string_append! observed_string "\n\t" expression_2_string(obs; strip_call_dict) ","
    end

    # Updates the string (including removing the last `,`) and returns it.
    return observed_string[1:(end - 1)] * "\n]"
end

# Creates an annotation for the system's observables.
function get_observed_annotation(rn::ReactionSystem)
    return "Observables:"
end

# Combines the 3 -related functions in a constant tuple.
OBSERVED_FS = (seri_has_observed, get_observed_string, get_observed_annotation)

### Handles Observables ###

# Checks if the reaction system has any defaults.
function seri_has_defaults(rn::ReactionSystem)
    return !isempty(initial_conditions(rn))
end

# Extract a string which declares the system's defaults.
function get_defaults_string(rn::ReactionSystem)
    defaults_string = "initial_conditions = " * x_2_string(MT.get_initial_conditions(rn))
    return defaults_string
end

# Creates an annotation for the system's defaults.
function get_defaults_annotation(rn::ReactionSystem)
    return "Defaults:"
end

# Combines the 3 defaults-related functions in a constant tuple.
DEFAULTS_FS = (seri_has_defaults, get_defaults_string, get_defaults_annotation)

### Handles Continuous Events ###

# Checks if the reaction system have has continuous events.
function seri_has_continuous_events(rn::ReactionSystem)
    return !isempty(MT.get_continuous_events(rn))
end

# Extract a string which declares the system's continuous events.
function get_continuous_events_string(rn::ReactionSystem)
    # Creates a dictionary for converting symbolics to their call-stripped form (e.g. X(t) to X).
    strip_call_dict = make_strip_call_dict(rn)

    # Handles the case with one event separately. Only effect is nicer formatting.
    if length(MT.get_continuous_events(rn)) == 1
        return "continuous_events = [$(continuous_event_string(MT.get_continuous_events(rn)[1], strip_call_dict))]"
    end

    # Creates the string corresponding to the code which generates the system's reactions.
    continuous_events_string = "continuous_events = ["
    for continuous_event in MT.get_continuous_events(rn)
        @string_append! continuous_events_string "\n\t" continuous_event_string(
            continuous_event, strip_call_dict
        ) ","
    end

    # Updates the string (including removing the last `,`) and returns it.
    return get_substring_end(continuous_events_string, 1, -1) * "\n]"
end

# Creates a string that corresponds to the declaration of a single continuous event.
function continuous_event_string(continuous_event, strip_call_dict)
    # Creates the string corresponding to the equations (i.e. conditions).
    eqs_string = "["
    for eq in continuous_event.conditions
        @string_append! eqs_string expression_2_string(eq; strip_call_dict) ", "
    end
    eqs_string = get_substring_end(eqs_string, 1, -2) * "]"

    # Creates the string corresponding to the affects.
    # Continuous events' `affect` field should probably be called `affects`. Likely the `s` was
    # dropped by mistake in MTK.
    affects_string = "["
    for affect in continuous_event.affect.affect
        @string_append! affects_string expression_2_string(affect; strip_call_dict) ", "
    end
    affects_string = get_substring_end(affects_string, 1, -2) * "]"

    # If there are discrete parameters, these have to be added and the full constructor used.
    event = eqs_string * " => " * affects_string
    if !isempty(continuous_event.affect.discrete_parameters)
        dp_string = expression_2_string(continuous_event.affect.discrete_parameters; strip_call_dict)
        dp_string = dp_string[findfirst('[', dp_string):end]
        event = "ModelingToolkitBase.SymbolicContinuousCallback($event; discrete_parameters = $dp_string)"
    end
    return event
end

# Creates an annotation for the system's continuous events.
function get_continuous_events_annotation(rn::ReactionSystem)
    return "Continuous events:"
end

# Combines the 3 -related functions in a constant tuple.
CONTINUOUS_EVENTS_FS = (seri_has_continuous_events, get_continuous_events_string, get_continuous_events_annotation)

### Handles Discrete Events ###

# Checks if the reaction system has any discrete events.
function seri_has_discrete_events(rn::ReactionSystem)
    return !isempty(MT.get_discrete_events(rn))
end

# Extract a string which declares the system's discrete events.
function get_discrete_events_string(rn::ReactionSystem)
    # Creates a dictionary for converting symbolics to their call-stripped form (e.g. X(t) to X).
    strip_call_dict = make_strip_call_dict(rn)

    # Handles the case with one event separately. Only effect is nicer formatting.
    if length(MT.get_discrete_events(rn)) == 1
        return "discrete_events = [$(discrete_event_string(MT.get_discrete_events(rn)[1], strip_call_dict))]"
    end

    # Creates the string corresponding to the code which generates the system's reactions.
    discrete_events_string = "discrete_events = ["
    for discrete_event in MT.get_discrete_events(rn)
        @string_append! discrete_events_string "\n\t" discrete_event_string(
            discrete_event, strip_call_dict
        ) ","
    end

    # Updates the string (including removing the last `,`) and returns it.
    return get_substring_end(discrete_events_string, 1, -1) * "\n]"
end

# Creates a string that corresponds to the declaration of a single discrete event.
function discrete_event_string(discrete_event, strip_call_dict)
    # Creates the string corresponding to the conditions. The special check is if the condition is
    # an expression like `X > 5.0`. Here, "(...)" is added for purely aesthetic reasons.
    condition_string = x_2_string(discrete_event.conditions)
    if discrete_event.conditions isa SymbolicT
        @string_prepend! "(" condition_string
        @string_append! condition_string ")"
    end

    # Creates the string corresponding to the affects.
    affects_string = "["
    for affect in discrete_event.affect.affect
        @string_append! affects_string expression_2_string(affect; strip_call_dict) ", "
    end
    affects_string = get_substring_end(affects_string, 1, -2) * "]"

    # If there are discrete parameters, these have to be added and the full constructor used.
    event = condition_string * " => " * affects_string
    if !isempty(discrete_event.affect.discrete_parameters)
        dp_string = expression_2_string(discrete_event.affect.discrete_parameters; strip_call_dict)
        dp_string = dp_string[findfirst('[', dp_string):end]
        event = "ModelingToolkitBase.SymbolicDiscreteCallback($event; discrete_parameters = $dp_string)"
    end
    return event
end

# Creates an annotation for the system's discrete events.
function get_discrete_events_annotation(rn::ReactionSystem)
    return "Discrete events:"
end

# Combines the 3 -related functions in a constant tuple.
DISCRETE_EVENTS_FS = (seri_has_discrete_events, get_discrete_events_string, get_discrete_events_annotation)

### Handles Brownian Types ###

# Checks if the reaction system has any brownian types.
function seri_has_brownian_type(rn::ReactionSystem)
    return false
end

# Extract a string which declares the system's brownian types.
function get_brownian_type_string(rn::ReactionSystem)
    return get_unsupported_comp_string("brownian types")
end

# Creates an annotation for the system's brownian types.
function get_brownian_type_annotation(rn::ReactionSystem)
    return get_unsupported_comp_annotation("Brownian types:")
end

# Combines the 3 brownian types-related functions in a constant tuple.
BROWNIAN_TYPE_FS = (seri_has_brownian_type, get_brownian_type_string, get_brownian_type_annotation)

### Handles Jump Types ###

# Checks if the reaction system has any jump types.
function seri_has_jump_type(rn::ReactionSystem)
    return false
end

# Extract a string which declares the system's jump types.
function get_jump_type_string(rn::ReactionSystem)
    return get_unsupported_comp_string("jump types")
end

# Creates an annotation for the system's jump types.
function get_jump_type_annotation(rn::ReactionSystem)
    return get_unsupported_comp_annotation("Jump types:")
end

# Combines the 3 jump types-related functions in a constant tuple.
JUMP_TYPE_FS = (seri_has_jump_type, get_jump_type_string, get_jump_type_annotation)

### Handles Systems ###

# Specific `push_field` function, which is used for the system field (where the annotation option
# must be passed to the `get_component_string` function). Since non-ReactionSystem systems cannot be
# written to file, this function throws an error if any such systems are encountered.
function push_systems_field(
        file_text::String, rn::ReactionSystem, annotate::Bool,
        top_level::Bool
    )
    # Checks whether there are any subsystems, and if these are ReactionSystems.
    seri_has_systems(rn) || (return (file_text, false))
    if any(!(system isa ReactionSystem) for system in MT.get_systems(rn))
        error("Tries to write a ReactionSystem to file which have non-ReactionSystem subs-systems. This is currently not possible.")
    end

    # Adds the system declaration string to the file string.
    write_string = "\n"
    top_level || (@string_append! write_string "local ")
    @string_append! write_string get_systems_string(rn, annotate)
    annotate && (@string_prepend! "\n\n# " get_systems_annotation(rn) write_string)
    return (file_text * write_string, true)
end

# Checks if the reaction system has any systems.
function seri_has_systems(rn::ReactionSystem)
    return !isempty(MT.get_systems(rn))
end

# Extract a string which declares the system's systems.
function get_systems_string(rn::ReactionSystem, annotate::Bool)
    # Initiates the `systems` string. It is pre-declared vector, into which the systems are added.
    systems_string = "systems = Vector(undef, $(length(MT.get_systems(rn))))"

    # Loops through all systems, adding their declaration to the system string.
    for (idx, system) in enumerate(MT.get_systems(rn))
        if annotate
            @string_append! systems_string "\n\n# Declares subsystem: $(getname(system))"
        end

        # Manipulates the subsystem declaration to make it nicer.
        subsystem_string = get_full_system_string(system, annotate, false)
        subsystem_string = replace(subsystem_string, "\n" => "\n\t")
        subsystem_string = "let\n" * get_substring_end(subsystem_string, 7, -6) * "end"
        @string_append! systems_string "\nsystems[$idx] = " subsystem_string
    end

    return systems_string
end

# Creates an annotation for the system's systems.
function get_systems_annotation(rn::ReactionSystem)
    return "Subsystems:"
end

# Combines the 3 systems-related functions in a constant tuple.
SYSTEMS_FS = (seri_has_systems, get_systems_string, get_systems_annotation)
