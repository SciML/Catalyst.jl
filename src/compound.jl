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

### Balancing Code

function create_matrix(reaction::Catalyst.Reaction)
    compounds = [reaction.substrates; reaction.products]
    atoms = []
    n_atoms = 0
    A = zeros(Int, 0, length(compounds))

    for (j, compound) in enumerate(compounds)
        if iscompound(compound)
            pairs = component_coefficients(compound)
            if pairs == Nothing 
                continue
            end
        else 
            pairs = [(compound, 1)]
        end

        for pair in pairs
            atom, coeff = pair
            i = findfirst(x -> isequal(x, atom), atoms)
            if i === nothing  
                push!(atoms, atom)
                n_atoms += 1
                A = [A; zeros(Int, 1, length(compounds))]
                i = n_atoms
            end
            coeff = any(map(p -> isequal(p, compounds[j]), reaction.products)) ? -coeff : coeff
            A[i, j] = coeff
        end
    end
        # Append a row with last element as 1 
        new_row = zeros(Int, 1, size(A, 2))
        new_row[end] = 1
        A = vcat(A, new_row)
    return A
end

function get_stoich(reaction::Reaction)
    # Create the matrix A using create_matrix function.
    A = create_matrix(reaction)
    
    # Create the vector b. The last element is 1 and others are 0.
    B = zeros(Int64, size(A,1))
    B[end] = 1

    # Concatenate B to A to form an augmented matrix   
    AB = [A B]

    # Apply the Bareiss algorithm
    ModelingToolkit.bareiss!(AB)

    # Extract the transformed A and B
    A_transformed = AB[:, 1:end-1]
    B_transformed = AB[:, end]

    # Convert A_transformed to rational numbers
    A_transformed_rational = Rational.(A_transformed)

    # Convert B_transformed to rational numbers
    B_transformed_rational = Rational.(B_transformed)

    # Perform the division
    X = A_transformed_rational \ B_transformed_rational

    # Get the denominators of the rational numbers in X
    denominators = denominator.(X)

    # Compute the LCM of the denominators
    lcm_value = reduce(lcm, denominators)

    # Multiply each element in X by the LCM of the denominators
    X_multiplied = X .* lcm_value

    # Convert the rational numbers to integers
    X_integers = numerator.(X_multiplied)

    return X_integers
end

function balance(reaction::Reaction)
    # Calculate the stoichiometric coefficients for the balanced reaction.
    stoichiometries = get_stoich(reaction)

    # Divide the stoichiometry vector into substrate and product stoichiometries.
    substoich = stoichiometries[1:length(reaction.substrates)]
    prodstoich = stoichiometries[(length(reaction.substrates)) + 1:end]

    # Create a new reaction with the balanced stoichiometries
    balanced_reaction = Reaction(reaction.rate, reaction.substrates, reaction.products, substoich, prodstoich)

    # Return the balanced reaction
    return balanced_reaction
end

# function get_stoich(reaction::Reaction) ## Gaussian Elimination
#     # Create the matrix A using create_matrix function.
#     A = create_matrix(reaction)
    
#     # Create the vector b. The last element is 1 and others are 0.
#     B = zeros(Int64, size(A,1))
#     B[end] = 1

#     # Convert A and B to sparse matrices
#     A = sparse(A)

#     # Apply the LU factorization (Gaussian elimination)
#     lu_f = lu(A)

#     # Convert B to a dense vector
#     B_dense = Vector(B)

#     # Solve for X using the LU factorization
#     X = lu_f \ B_dense

#     # Get the smallest positive value in X
#     smallest_value = minimum(X[X .> 0])  

#     # Normalize X to be integers
#     X_normalized = round.(Int64, X / smallest_value)

#     return X_normalized
# end

# function get_stoich(reaction::Reaction) ##Bareiss Algorithm (Does not avoid floating points)
#     # Create the matrix A using create_matrix function.
#     A = create_matrix(reaction)
    
#     # Create the vector b. The last element is 1 and others are 0.
#     B = zeros(Int64, size(A,1))
#     B[end] = 1

#     # Concatenate B to A to form an augmented matrix   
#     AB = [A B]

#     # Apply the Bareiss algorithm
#     ModelingToolkit.bareiss!(AB)

#     # Extract the transformed A and B
#     A_transformed = AB[:, 1:end-1]
#     B_transformed = AB[:, end]

#     # Solve for X using back substitution
#     X = A_transformed \ B_transformed

#     # Get the smallest positive value in X
#     smallest_value = minimum(X[X .> 0])  

#     # Normalize X to be integers
#     X_normalized = round.(Int64, X / smallest_value)

#     return X_normalized
# end

# @variables t
# @parameters k
# @species C(t) H(t) O(t) 
# @compound O2(t) 2O
# @compound CO2(t) 1C 2O
# @compound H2O(t) 2H 1O
# @compound C6H12O6(t) 6C 12H 6O
# rx = Reaction(k,[CO2,H2O],[C6H12O6,O2])

# using LinearAlgebra
# using SparseArrays
# using SuiteSparse.UMFPACK

# gcd_value = gcd(X...)
# X_scaled = X ./ gcd_value

# X = [1,1,0.2,1]

