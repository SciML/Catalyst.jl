### Field Serialisation Support Functions ###

# Function which handles the addition of a single component to the file string.
function push_component(file_text::String, rn::ReactionSystem, annotate::Bool, comp_funcs::Tuple)
    has_component, get_comp_string, get_comp_annotation = comp_funcs
    has_component(rn) || (return (file_text, false))
    write_string = "\n" * get_comp_string(rn)
    annotate && (write_string = "\n\n# " * get_comp_annotation(rn) * write_string)
    return (file_text * write_string, true)
end

# Generic function for creating an string for an unsupported argument.
function get_unsupported_comp_string(component::String)
    @warn "Writing ReactionSystem models with $(component) is currently not supported. This field is not written to the file."
    return ""
end

# Generic function for creating the annotation string for an unsupported argument.
function get_unsupported_comp_annotation(component::String)
    return "$(component): (OBS: Currently not supported, and hence empty)"
end


### Symbolics String Conversions ###

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

# Converts a numeric expression to a string. `uncall_dict` is required to convert stuff like 
# `p1*X(t) + p2` to `p1*X + p2`.
function num_2_string(num; uncall_dict = make_uncall_dict(Symbolics.get_variables(num)))
    uncalled_num = substitute(num, uncall_dict)
    return repr(uncalled_num)
end

### Generic Expression Handling ###

# Function for converting a Symbolics of the form `X(t)` to the form `X`.
function uncall(var) 
    istree(var) ? Sym{Real}(Symbolics.getname(var)) : var
end

# For a list of symbolics, makes an "uncall dict". It permits the conversion of e.g. `X(t)` to `X`.
function make_uncall_dict(syms)
    return Dict([sym => uncall(Symbolics.unwrap(sym)) for sym in syms])
end
