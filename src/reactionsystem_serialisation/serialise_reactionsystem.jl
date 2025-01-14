"""
    save_reactionsystem(filename::String, rn::ReactionSystem; annotate = true, safety_check = true)

Save a `ReactionSystem` model to a file. The `ReactionSystem` is saved as runnable Julia code. This
can both be used to save a `ReactionSystem` model, but also to write it to a file for easy inspection.

Arguments:
- `filename`: The name of the file to which the `ReactionSystem` is saved.
- `rn`: The `ReactionSystem` which should be saved to a file.
- `annotate = true`: Whether annotation should be added to the file.
- `safety_check = true`: After serialisation, Catalyst will automatically load the serialised
  `ReactionSystem` and check that it is equal to `rn`. If it is not, an error will be thrown. For
  models without the `connection_type` field, this should not happen. If performance is required
  (i.e. when saving a large number of models), this can be disabled by setting `safety_check = false`.

Example:
```julia
rn = @reaction_network begin
    (p,d), 0 <--> X
end
save_reactionsystem("rn.jls", rn)
```
The model can now be loaded using
```julia
rn = include("rn.jls")
```

Notes:
- `ReactionSystem`s with the `connection_type` field has this ignored (saving of this field has not
  been implemented yet).
- `ReactionSystem`s with non-`ReactionSystem` sub-systems (e.g. `ODESystem`s) cannot be saved.
- Reaction systems with components that have units cannot currently be saved.
- The `ReactionSystem` is saved using *programmatic* (not DSL) format for model creation.
"""
function save_reactionsystem(filename::String, rn::ReactionSystem;
        annotate = true, safety_check = true)
    # Error and warning checks.
    reactionsystem_uptodate_check()
    if !isempty(get_networkproperties(rn))
        @warn "The serialised network has cached network properties (e.g. computed conservation laws). This will not be saved as part of the network, and must be recomputed when it is loaded."
    end

    # Write model to file and performs a safety check.
    open(filename, "w") do file
        write(file, get_full_system_string(rn, annotate, true))
    end
    if safety_check
        if !isequal(rn, include(joinpath(pwd(), filename)))
            rm(filename)
            error("The serialised `ReactionSystem` is not equal to the original one. Please make a report (including the full system) at https://github.com/SciML/Catalyst.jl/issues. To disable this behaviour, please pass the `safety_check = false` argument to `save_reactionsystem` (warning, this will permit the serialisation of an erroneous system).")
        end
    end
    return nothing
end

# Gets the full string which corresponds to the declaration of a system. Might be called recursively
# for systems with subsystems.
function get_full_system_string(rn::ReactionSystem, annotate::Bool, top_level::Bool)
    # MTK automatically performs flattening when `complete` is carried out. In case of a
    # hierarchical system, we must undo this an process the non-flattened system.
    iscomplete(rn) && (rn = MT.get_parent(rn))

    # Initiates the file string.
    file_text = ""

    # Goes through each type of system component, potentially adding it to the string.
    # Species, variables, and parameters must be handled differently in case there are default values
    # dependencies between them.
    # Systems use custom `push_field` function as these require the annotation `Bool`to be passed
    # to the function that creates the next sub-system declarations.
    file_text, _ = push_field(file_text, rn, annotate, top_level, IV_FS)
    file_text, has_sivs = push_field(file_text, rn, annotate, top_level, SIVS_FS)
    file_text, has_parameters, has_species, has_variables = handle_us_n_ps(
        file_text, rn, annotate, top_level)
    file_text, has_reactions = push_field(file_text, rn, annotate, top_level, REACTIONS_FS)
    file_text, has_equations = push_field(file_text, rn, annotate, top_level, EQUATIONS_FS)
    file_text, has_observed = push_field(file_text, rn, annotate, top_level, OBSERVED_FS)
    file_text, has_defaults = push_field(file_text, rn, annotate, top_level, DEFAULTS_FS)
    file_text, has_continuous_events = push_field(file_text, rn, annotate,
        top_level, CONTINUOUS_EVENTS_FS)
    file_text, has_discrete_events = push_field(file_text, rn, annotate,
        top_level, DISCRETE_EVENTS_FS)
    file_text, has_systems = push_systems_field(file_text, rn, annotate, top_level)
    file_text, has_connection_type = push_field(file_text, rn, annotate,
        top_level, CONNECTION_TYPE_FS)

    # Finalise the system. Creates the final `ReactionSystem` call.
    # Enclose everything in a `let ... end` block.
    rs_creation_code = make_reaction_system_call(
        rn, annotate, top_level, has_sivs, has_species,
        has_variables, has_parameters, has_reactions,
        has_equations, has_observed, has_defaults, has_continuous_events,
        has_discrete_events, has_systems, has_connection_type)
    annotate || (@string_prepend! "\n" file_text)
    @string_prepend! "let" file_text
    @string_append! file_text "\n\n" rs_creation_code "\n\nend"

    return file_text
end

# Creates a ReactionSystem call for creating the model. Adds all the correct inputs to it. The input
# `has_` `Bool`s described which inputs are used. If the model is `complete`, this is handled here.
function make_reaction_system_call(rs::ReactionSystem, annotate, top_level, has_sivs,
        has_species, has_variables, has_parameters, has_reactions, has_equations,
        has_observed, has_defaults, has_continuous_events, has_discrete_events, has_systems,
        has_connection_type)

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
    has_defaults && (@string_append! reaction_system_string ", defaults")
    has_continuous_events && (@string_append! reaction_system_string ", continuous_events")
    has_discrete_events && (@string_append! reaction_system_string ", discrete_events")
    has_systems && (@string_append! reaction_system_string ", systems")
    has_connection_type && (@string_append! reaction_system_string ", connection_type")

    # Potentially appends a combinatoric_ratelaws statement.
    if !Symbolics.unwrap(combinatoric_ratelaws(rs))
        @string_append! reaction_system_string ", combinatoric_ratelaws = false"
    end

    # Potentially appends `ReactionSystem` metadata value(s). Weird composite types are not supported.
    if !isnothing(MT.get_metadata(rs))
        @string_append! reaction_system_string ", metadata = $(x_2_string(MT.get_metadata(rs)))"
    end

    # Finalises the call. Appends potential annotation. If the system is complete, add a call for this.
    @string_append! reaction_system_string ")"
    if ModelingToolkit.iscomplete(rs)
        @string_prepend! "rs = " reaction_system_string
        top_level || (@string_prepend! "local " reaction_system_string)
        @string_append! reaction_system_string "\ncomplete(rs)"
    end
    if annotate
        @string_prepend! "# Declares ReactionSystem model:\n" reaction_system_string
    end
    return reaction_system_string
end
