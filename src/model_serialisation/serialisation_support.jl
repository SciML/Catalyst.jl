### Generic Functions ###

# Function which handles the addition of a single component to the file string.
function push_component(file_text::String, rn::ReactionSystem, annotate::Bool, comp_funcs::Tuple)
    has_component, get_comp_string, get_comp_annotation = comp_funcs
    has_component(rn) || (return (file_text, false))
    write_string = "\n" * get_comp_string(rn)
    annotate && (write_string = "\n\n# " * get_comp_annotation(rn) * write_string)
    return (file_text * write_string, true)
end

# Converts a Num to a String.
function wrapped_num_2_string(num) 
    return String(Symbol(num))
end

# Converts a vector of Nums to a single string with all (e.g. [X(t), Y(t)] becomes " X(t) Y(t)").
function wrapped_nums_2_string(nums) 
    return prod([" " * wrapped_num_2_string(num) for num in nums])
end

# Converts a vector of Symbolics to a (uncalled) vector (e.g. [X(t), Y(t), p] becomes " [X, Y, p]").
function nums_2_uncalled_strin_vec(nums)
    return "$(convert(Vector{Any}, uncall.(nums)))"[4:end]
end

# Function for converting a Symbolics of the form `X(t)` to the form `X`.
function uncall(var) 
    istree(var) ? Sym{Real}(Symbolics.getname(var)) : var
end

function nums_2_uncalled_strin_vec(nums)
    "$(convert(Vector{Any}, uncall.(nums)))"[4:end]
end

# Converts a numeric expression to a string. `uncall_dict` is required to convert stuff like 
# `p1*X(t) + p2` to `p1*X + p2`.
function num_2_string(num; uncall_dict = make_uncall_dict(Symbolics.get_variables(num)))
    uncalled_num = substitute(num, uncall_dict)
    return repr(uncalled_num)
end

# For a list of symbolics, makes an "uncall dict". It permits the conversion of e.g. `X(t)` to `X`.
function make_uncall_dict(syms)
    return Dict([sym => uncall(Symbolics.unwrap(sym)) for sym in syms])
end

### Generic Unsupported Component FUnctions ###

# Generic function for creating an string for an unsupported argument.
function get_unsupported_comp_string(component::String)
    @warn "Writing ReactionSystem models with $(component) is currently not supported. This field is not written to the file."
    return ""
end

# Generic function for creating the annotation string for an unsupported argument.
function get_unsupported_comp_annotation(component::String)
    return "$(component): (OBS: Currently not supported, and hence empty)"
end


### Handles Independent Variables ###

# Checks if the reaction system have any independent variable. True for all valid reaction systems.
function has_iv(rn::ReactionSystem)
    return true
end

# Extract a string which declares the system's independent variable.
function get_iv_string(rn::ReactionSystem)
    return "@variables $(num_2_string(ModelingToolkit.get_iv(rn)))"
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
    return length(rn.sivs) != 0
end

# Extract a string which declares the system's spatial independent variables.
function get_sivs_string(rn::ReactionSystem)
    return "sivs = @variables$(wrapped_nums_2_string(rn.sivs))"
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
    return length(species(rn)) != 0
end

# Extract a string which declares the system's species.
function get_species_string(rn::ReactionSystem)
    return "sps = @species$(wrapped_nums_2_string(species(rn)))"
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
    return length(unknowns(rn)) > length(species(rn))
end

# Extract a string which declares the system's variables.
function get_variables_string(rn::ReactionSystem)
    variables = filter(!isspecies, unknowns(rn))
    return "vars = @variables$(wrapped_nums_2_string(variables))"
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
    return length(parameters(rn)) != 0
end

# Extract a string which declares the system's parameters.
function get_parameters_string(rn::ReactionSystem)
    return "ps = @parameters$(wrapped_nums_2_string(parameters(rn)))"
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
    uncall_dict = make_uncall_dict(unknowns(rn))
    rxs_string = "rxs = ["
    for rx in reactions(rn)
        rxs_string = rxs_string * "\n\t" * reaction_string(rx, uncall_dict) * ","
    end

    # Updates the string (including removing the last `,`) and returns in.
    return rxs_string[1:end-1] * "\n]"
end


# Creates a string that corresponds to the declaration of a single `Reaction`.
function reaction_string(rx::Reaction, uncall_dict)
    # Prepares the `Reaction` declaration components.
    rate = num_2_string(rx.rate; uncall_dict)
    substrates = isempty(rx.substrates) ? "nothing" : nums_2_uncalled_strin_vec(rx.substrates)    
    products = isempty(rx.products) ? "nothing" : nums_2_uncalled_strin_vec(rx.products)
    substoich = isempty(rx.substoich) ? "nothing" : nums_2_uncalled_strin_vec(rx.substoich)
    prodstoich = isempty(rx.prodstoich) ? "nothing" : nums_2_uncalled_strin_vec(rx.prodstoich)

    # Creates the full expression, including adding kwargs (`only_use_rate` and `metadata`).
    rx_string = "Reaction($rate, $(substrates), $(products), $(substoich), $(prodstoich)"
    if rx.only_use_rate
        rx_string = rx_string * "; only_use_rate = true"
    end
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
    return length(equations(rn)) > length(reactions(rn))
end

# Extract a string which declares the system's equations.
function get_equations_string(rn::ReactionSystem)
    get_unsupported_comp_string("equations")
end

# Creates an annotation for the system's equations.
function get_equations_annotation(rn::ReactionSystem)
    get_unsupported_comp_annotation("Equations")
end

# Combines the 3 equations-related functions in a constant tuple.
EQUATIONS_FS = (has_equations, get_equations_string, get_equations_annotation)