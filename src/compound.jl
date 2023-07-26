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
        if isa(expr, Expr) && expr.head == :call && expr.args[1] == :*
            push!(coeffs, expr.args[2])
            push!(species, expr.args[3])
        else
            push!(coeffs, 1)
            push!(species, expr)
        end
    end

    coeffs_expr = Expr(:vect, coeffs...)
    species_expr = Expr(:vect, species...)


    # Construct the expression to set the components metadata
    setcomponents_expr = :($(species_name) = ModelingToolkit.setmetadata($(species_name),
                                                                         Catalyst.CompoundComponents,
                                                                         $species_expr))

    # Ensure the expression is evaluated in the correct scope by escaping it
    escaped_setcomponents_expr = esc(setcomponents_expr)

    # Construct the expression to set the coefficients metadata
    setcoefficients_expr = :($(species_name) = ModelingToolkit.setmetadata($(species_name),
                                                                           Catalyst.CompoundCoefficients,
                                                                           $coeffs_expr))

    escaped_setcoefficients_expr = esc(setcoefficients_expr)

    # Return a block that contains the escaped expressions
    return Expr(:block, escaped_species_expr, escaped_setmetadata_expr,
                escaped_setcomponents_expr, escaped_setcoefficients_expr)
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
    return [c => co for (c, co) in zip(components(s), coefficients(s))]
end

# function esc_dollars!(ex)
#     if ex isa Expr
#         if ex.head == :call && ex.args[1] == :* && ex.args[3] isa Expr && ex.args[3].head == :$
#             coeff = ex.args[2]
#             value = eval(ex.args[3].args[1])
#             return coeff * value
#         else
#             for i in 1:length(ex.args)
#                 ex.args[i] = esc_dollars!(ex.args[i])
#             end
#         end
#     end
#     ex
# end

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
#                                                                        Catalyst.CompoundSpecies,
#                                                                        true))

#     # Ensure the expressions are evaluated in the correct scope by escaping them
#     escaped_species_expr = esc(species_expr)
#     escaped_setmetadata_expr = esc(setmetadata_expr)
    
#     # Construct the array from the remaining arguments
#     coeffs = []
#     species = []

#     for expr in arr_expr
#         if isa(expr, Expr) && expr.head == :call && expr.args[1] == :*
#             coeff = expr.args[2]
#             specie = expr.args[3]
#             # check if the coefficient is a variable and get its value
#             if isa(coeff, Expr) && coeff.head == :$
#                 coeff = getfield(Main, coeff.args[1])
#             end
#             push!(coeffs, coeff)
#             push!(species, specie)
#         else
#             push!(coeffs, 1)
#             push!(species, expr)
#         end
#     end

#     coeffs_expr = Expr(:vect, coeffs...)
#     species_expr = Expr(:vect, species...)

#     # Construct the expression to set the components metadata
#     setcomponents_expr = :($(species_name) = ModelingToolkit.setmetadata($(species_name),
#                                                                          Catalyst.CompoundComponents,
#                                                                          $species_expr))

#     # Ensure the expression is evaluated in the correct scope by escaping it
#     escaped_setcomponents_expr = esc(setcomponents_expr)

#     # Construct the expression to set the coefficients metadata
#     setcoefficients_expr = :($(species_name) = ModelingToolkit.setmetadata($(species_name),
#                                                                            Catalyst.CompoundCoefficients,
#                                                                            $coeffs_expr))

#     escaped_setcoefficients_expr = esc(setcoefficients_expr)

#     # Return a block that contains the escaped expressions
#     return Expr(:block, escaped_species_expr, escaped_setmetadata_expr,
#                 escaped_setcomponents_expr, escaped_setcoefficients_expr)
# end


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
#                                                                        Catalyst.CompoundSpecies,
#                                                                        true))

#     # Ensure the expressions are evaluated in the correct scope by escaping them
#     escaped_species_expr = esc(species_expr)
#     escaped_setmetadata_expr = esc(setmetadata_expr)

#     # Construct the array from the remaining arguments
#     arr = Expr(:vect, (arr_expr)...)

#     coeffs = []
#     species = []

#     # Add process_expr! function within the macro
#     function process_expr!(ex)
#         if ex isa Expr && ex.head == :call && ex.args[1] == :*
#             # Extract the coefficient and the species
#             coeff = ex.args[2]
#             species = ex.args[3]

#             # Check if the species is an interpolated variable
#             if species isa Expr && species.head == :$
#                 # Replace the interpolated variable with its value
#                 species = getfield(Main, species.args[1])

#                 # Multiply the coefficient with the species
#                 return Expr(:call, :*, coeff, species)
#             end
#         end
#         return ex
#     end

#     for expr in arr_expr
#         # Call process_expr! on the expression to process it
#         processed_expr = process_expr!(expr)
        
#         if isa(processed_expr, Expr) && processed_expr.head == :call && processed_expr.args[1] == :*
#             push!(coeffs, processed_expr.args[2])
#             push!(species, processed_expr.args[3])
#         else
#             push!(coeffs, 1)
#             push!(species, processed_expr)
#         end
#     end

#     coeffs_expr = Expr(:vect, coeffs...)
#     species_expr = Expr(:vect, species...)

#     # Construct the expression to set the components metadata
#     setcomponents_expr = :($(species_name) = ModelingToolkit.setmetadata($(species_name),
#                                                                          Catalyst.CompoundComponents,
#                                                                          $species_expr))

#     # Ensure the expression is evaluated in the correct scope by escaping it
#     escaped_setcomponents_expr = esc(setcomponents_expr)

#     # Construct the expression to set the coefficients metadata
#     setcoefficients_expr = :($(species_name) = ModelingToolkit.setmetadata($(species_name),
#                                                                            Catalyst.CompoundCoefficients,
#                                                                            $coeffs_expr))

#     escaped_setcoefficients_expr = esc(setcoefficients_expr)

#     # Return a block that contains the escaped expressions
#     return Expr(:block, escaped_species_expr, escaped_setmetadata_expr,
#                 escaped_setcomponents_expr, escaped_setcoefficients_expr)
# end
