const MT = ModelingToolkit
struct CompoundSpecies end
struct CompoundComponents end
struct CompoundCoefficients end

Symbolics.option_to_metadata_type(::Val{:iscompound}) = CompoundSpecies
Symbolics.option_to_metadata_type(::Val{:components}) = CompoundComponents
Symbolics.option_to_metadata_type(::Val{:coefficients}) = CompoundCoefficients

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
    coeffs = []
    species = []
    for expr in arr_expr
        if expr.head == :call && expr.args[1] == :*
            push!(coeffs, expr.args[2])
            push!(species, expr.args[3])
        else
            push!(coeffs, 1)
            push!(species, expr)
        end
    end
    coeffs_expr = Expr(:vect, coeffs...)
    species_expr = Expr(:vect, species...)

    # Ensure the expression is evaluated in the correct scope by escaping it
    escaped_setcomponentcoefficients_expr = esc(setcomponentcoefficients_expr)

    # Construct the expression to set the components metadata
    setcomponents_expr = :($(species_name) = ModelingToolkit.setmetadata($(species_name),
                                                                         Catalyst.CompoundComponents,
                                                                         $species_expr))

    # Ensure the expression is evaluated in the correct scope by escaping it
    escaped_setcomponents_expr = esc(setcomponents_expr)

    # Construct the expression to set the components metadata
    setcoefficients_expr = :($(species_name) = ModelingToolkit.setmetadata($(species_name),
                                                                         Catalyst.CompoundCoefficients,
                                                                         $coeffs_expr))

    escaped_setcoefficients_expr = esc(setcoefficients_expr)

    # Return a block that contains the escaped expressions
    return Expr(:block, escaped_species_expr, escaped_setmetadata_expr,
                escaped_setcomponents_expr, escaped_setcoefficients_expr, escaped_setcomponentcoefficients_expr)
end

# Check if a species is a compound
iscompound(s::Num) = iscompound(MT.value(s))
function iscompound(s)
    MT.getmetadata(s, CompoundSpecies, false)
end

coefficients(s::Num) = coefficients(MT.value(s))
function coefficients(s)
    MT.getmetadata(s, CompoundCoefficients)
end

components(s::Num) = components(MT.value(s))
function components(s)
    MT.getmetadata(s, CompoundComponents)
end

component_coefficients(s::Num) = component_coefficients(MT.value(s))
function component_coefficients(s)
    MT.getmetadata(s, CompoundComponentCoefficients)
end

component_coefficients(s::Num) = component_coefficients(MT.value(s))
function component_coefficients(s)
    return [c => co for (c, co) in zip(components(s), coefficients(s))]
end

# const MT = ModelingToolkit
# struct CompoundSpecies end
# struct CompoundComponents end
# struct CompoundCoefficients end
# struct CompoundComponentCoefficients end

# Symbolics.option_to_metadata_type(::Val{:iscompound}) = CompoundSpecies
# Symbolics.option_to_metadata_type(::Val{:components}) = CompoundComponents
# Symbolics.option_to_metadata_type(::Val{:coefficients}) = CompoundCoefficients
# Symbolics.option_to_metadata_type(::Val{:component_coefficients}) = CompoundComponentCoefficients

# macro compound(species_expr, arr_expr...)
#     # Ensure the species name is a valid expression
#     if !(species_expr isa Expr && species_expr.head == :call)
#         error("Invalid species name in @compound macro")
#     end

#     # Parse the species name to extract the species name and argument
#     species_name = species_expr.args[1]
#     species_arg = species_expr.args[2]

#     # Construct the expressions that define the species
#     species_expr = Expr(:macrocall, Symbol("@species"), LineNumberNode(0),
#                         Expr(:call, species_name, species_arg))

#     # Construct the expression to set the iscompound metadata
#     setmetadata_expr = :($(species_name) = ModelingToolkit.setmetadata($(species_name),
#                                                                        CompoundSpecies,
#                                                                        true))

#     # Ensure the expressions are evaluated in the correct scope by escaping them
#     escaped_species_expr = esc(species_expr)
#     escaped_setmetadata_expr = esc(setmetadata_expr)
#     escaped_arr_expr = esc(arr_expr)

#     # Parse each compound expression to extract the coefficients and components
#     coefficients_and_components = [parse_compound(expr) for expr in arr_expr]
    

    
#     # Extract coefficients and components from parsed result
#     coefficients = [cc[1] for cc in coefficients_and_components]
#     components = [cc[2] for cc in coefficients_and_components]
#     component_coefficients = [c => co for (c, co) in zip(components, coefficients)]


#     # Construct the expression to set the coefficients metadata
#     setcoefficients_expr = :($(species_name) = ModelingToolkit.setmetadata($(species_name),
#                                                                            CompoundCoefficients,
#                                                                            $coefficients))

#     # Ensure the expression is evaluated in the correct scope by escaping it
#     escaped_setcoefficients_expr = esc(setcoefficients_expr)

#     # Construct the expression to set the components metadata
#     setcomponents_expr = :($(species_name) = ModelingToolkit.setmetadata($(species_name),
#                                                                          CompoundComponents,
#                                                                          ($components)))

#     # Ensure the expression is evaluated in the correct scope by escaping it
#     escaped_setcomponents_expr = esc(setcomponents_expr)
    
#     # Construct the expression to set the component_coefficients metadata
#     setcomponent_coefficients_expr = :($(species_name) = ModelingToolkit.setmetadata($(species_name),
#                                                                            CompoundComponentCoefficients,
#                                                                            $component_coefficients)) 

#     # Ensure the expression is evaluated in the correct scope by escaping it
#     escaped_component_coefficients_expr = esc(setcomponent_coefficients_expr)

#     # Return a block that contains the escaped expressions
#     return Expr(:block, 
#                 escaped_species_expr,
#                 escaped_arr_expr,
#                 escaped_setmetadata_expr,
#                 escaped_setcoefficients_expr,
#                 escaped_setcomponents_expr,
#                 escaped_component_coefficients_expr
#                )
# end

# function parse_compound(compound_expr)
#     esc(compound_expr)
#     # Ensure compound_expr is a multiplication
#     Data = eval(compound_expr)
    
#     # extract the coefficient and component
#     coefficient = Symbolics.arguments(Symbolics.value(Data))[1]   # Integer
#     component = get_variables(Data)     # Species variable
    
#     return vcat(coefficient, component)
# end

# # Check if a species is a compound
# iscompound(s::Num) = iscompound(MT.value(s))
# function iscompound(s)
#     MT.getmetadata(s, CompoundSpecies, false)
# end

# coefficients(s::Num) = coefficients(MT.value(s))
# function coefficients(s)
#     MT.getmetadata(s, CompoundCoefficients)
# end

# components(s::Num) = components(MT.value(s))
# function components(s)
#     MT.getmetadata(s, CompoundComponents)
# end

# component_coefficients(s::Num) = component_coefficients(MT.value(s))
# function component_coefficients(s)
#     MT.getmetadata(s, CompoundComponentCoefficients)
# end


