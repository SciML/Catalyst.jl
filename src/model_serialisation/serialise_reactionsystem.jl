function save_reaction_network(filename::String, rn::ReactionSystem; annotate = true)
    # Initiates the file string.
    file_text = ""

    # Goes through each type of system component, potentially adding it to the string.
    file_text, has_iv = push_component(file_text, rn, annotate, IV_FS)
    file_text, has_sivs = push_component(file_text, rn, annotate, SIVS_FS)
    file_text, has_species = push_component(file_text, rn, annotate, SPECIES_FS)
    file_text, has_variables = push_component(file_text, rn, annotate, VARIABLES_FS)
    file_text, has_parameters = push_component(file_text, rn, annotate, PARAMETERS_FS)
    file_text, has_reactions = push_component(file_text, rn, annotate, REACTIONS_FS)
    file_text, has_equations = push_component(file_text, rn, annotate, EQUATIONS_FS)
   
    # Finalises the system. Creates the final `ReactionSystem` call.
    rs_creation_code = make_reaction_system_call(rn, file_text, annotate,
                                                 has_sivs, has_species, has_variables, 
                                                 has_parameters, has_reactions, has_equations)
    annotate || (file_text = "\n" * file_text)
    file_text = "let" * file_text * "\n\n" * rs_creation_code * "\n\nend"

    # Writes the model to a file. Then, returns nothing.
    open(filename, "w") do file
        write(file, file_text)
    end
    return nothing
end

# Takes the actual text which creates the model, and wraps it in a `let ... end` statement and A
# ReactionSystem call, to create the final file text.
function make_reaction_system_call(rs::ReactionSystem, file_text, annotate, has_sivs, has_species, has_variables, has_parameters, has_reactions, has_equations)
    # Creates base call.
    iv = Catalyst.get_iv(rs)
    if has_reactions && has_equations
        eqs = "[rxs; eqs]"
    elseif has_reactions
        eqs = "rxs"
    elseif has_equations
        eqs = "eqs"
    else 
        eqs = "[]"
    end
    if has_species && has_variables
        unknowns = "[sps; vars]"
    elseif has_species
        unknowns = "sps"
    elseif has_variables
        unknowns = "vars"
    else 
        unknowns = "[]"
    end
    if has_parameters
        ps = "ps"
    else 
        ps = "[]"
    end
    reaction_system_string = "ReactionSystem($eqs, $iv, $unknowns, $ps"

    # Appends additional (optional) arguments.
    if Base.isidentifier(Catalyst.getname(rs))
        rs_name = ":$(Catalyst.getname(rs))"
    else
        rs_name = "Symbol(\"$(Catalyst.getname(rs))\")"
    end
    reaction_system_string = reaction_system_string * "; name = $(rs_name)"
    reaction_system_string = reaction_system_string * ")"

    # Returns the full call.
    if !ModelingToolkit.iscomplete(rs)
        reaction_system_string = "rs = $(reaction_system_string)\ncomplete(rs)"
    end
    if annotate 
        reaction_system_string = "# Declares ReactionSystem model:\n" * reaction_system_string
    end
    return reaction_system_string
end