struct CompoundSpecies end
Symbolics.option_to_metadata_type(::Val{:iscompound}) = CompoundSpecies

macro compound(species_expr, arr_expr...)
    # Ensure the species name is a valid expression
    if !(species_expr isa Expr && species_expr.head == :call)
        error("Invalid species name in @compound macro")
    end

    # Parse the species name to extract the species name and argument
    species_name = species_expr.args[1]
    species_arg = species_expr.args[2]

    # Construct the expressions that define the species
    species_expr = Expr(:macrocall, Symbol("@species"), LineNumberNode(0), Expr(:call, species_name, species_arg))

    # Construct the expression to set the iscompound metadata
    setmetadata_expr = :($(species_name) = ModelingToolkit.setmetadata($(species_name), CompoundSpecies, true))

    # Ensure the expressions are evaluated in the correct scope by escaping them
    escaped_species_expr = esc(species_expr)
    escaped_setmetadata_expr = esc(setmetadata_expr)

    # Construct the array from the remaining arguments
    arr = Expr(:vect, map(esc, arr_expr)...)

    # Return a block that contains the escaped expressions
    return Expr(:block, escaped_species_expr, escaped_setmetadata_expr, :(arr = $arr))
end

iscompound(s::Num) = iscompound(ModelingToolkit.value(s))
function iscompound(s)
    getmetadata(s, CompoundSpecies, false)
end


# @variables t
# @parameters k
# @species C(t) H(t) O(t)
# @compound C6H12O2(t) 6C 12H 2O

# iscompound(C6H12O2)
# typeof(C6H12O2)