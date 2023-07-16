const MT = ModelingToolkit
struct CompoundSpecies end
struct CompoundData end

Symbolics.option_to_metadata_type(::Val{:iscompound}) = CompoundSpecies
Symbolics.option_to_metadata_type(::Val{:components}) = CompoundData
Symbolics.option_to_metadata_type(::Val{:coefficients}) = CompoundData

macro compound(species_expr, arr_expr...)
    # Ensure the species name is a valid expression
    if !(species_expr isa Expr && species_expr.head == :call)
        error("Invalid species name in @compound macro")
    end

    # Parse the species name to extract the species name and argument
    species_name = species_expr.args[1]
    species_arg = species_expr.args[2]

    # Construct the expressions that define the species
    species_expr = Expr(:macrocall, Symbol("@species"), LineNumberNode(0),
                        Expr(:call, species_name, species_arg))

    # Construct the expression to set the iscompound metadata
    setmetadata_expr = :($(species_name) = ModelingToolkit.setmetadata($(species_name),
                                                                       Catalyst.CompoundSpecies,
                                                                       true))

    # Ensure the expressions are evaluated in the correct scope by escaping them
    escaped_species_expr = esc(species_expr)
    escaped_setmetadata_expr = esc(setmetadata_expr)

    # Construct the array from the remaining arguments
    arr = Expr(:vect, (arr_expr)...)

    # Construct the expression to set the components metadata
    setcomponents_expr = :($(species_name) = ModelingToolkit.setmetadata($(species_name),
                                                                         Catalyst.CompoundData,
                                                                         $arr))

    # Ensure the expression is evaluated in the correct scope by escaping it
    escaped_setcomponents_expr = esc(setcomponents_expr)

    # Return a block that contains the escaped expressions
    return Expr(:block, escaped_species_expr, escaped_setmetadata_expr,
                escaped_setcomponents_expr)
end

# Check if a species is a compound
iscompound(s::Num) = iscompound(MT.value(s))
function iscompound(s)
    MT.getmetadata(s, CompoundSpecies, false)
end

coefficients(s::Num) = coefficients(MT.value(s))
function coefficients(s)
    Z = MT.getmetadata(s, CompoundData)
    coefficients = [Symbolics.arguments(Symbolics.value(element))[1] for element in Z]
    return coefficients
end

components(s::Num) = components(MT.value(s))
function components(s)
    Z = MT.getmetadata(s, CompoundData)
    comp = Num[]

    for component in Z
        # extract the variables using the get_variables() function
        extracted = get_variables(component)

        # append the extracted variables to the new array
        append!(comp, extracted)
    end

    return comp
end

component_coefficients(s::Num) = component_coefficients(MT.value(s))
function component_coefficients(s)
    return [c => co for (c, co) in zip(components(s), coefficients(s))]
end
