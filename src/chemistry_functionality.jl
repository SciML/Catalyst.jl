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
const COMPOUND_OF_COMPOUND_ERROR = ErrorException("Reaction balancing does not currently work for reactions involving compounds of compounds.")

# note this does not correctly handle compounds of compounds currently
function create_matrix(reaction::Catalyst.Reaction)
    @unpack substrates, products = reaction
    unique_atoms = [] # Array to store unique atoms
    n_atoms = 0
    ncompounds = length(substrates) + length(products)
    A = zeros(Int, 0, ncompounds)

    coeffsign = 1
    jbase = 0
    for compounds in (substrates, products)
        for (j2, compound) in enumerate(compounds)
            j = jbase + j2

            if iscompound(compound)
                atoms = components(compound)
                any(iscompound, atoms) && throw(COMPOUND_OF_COMPOUND_ERROR)
                coeffs = coefficients(compound)
                (atoms == nothing || coeffs == nothing) && continue
            else
                # If not a compound, assume coefficient of 1
                atoms = [compound]
                coeffs = [1]
            end

            for (atom,coeff) in zip(atoms, coeffs)
                # Extract atom and coefficient from the pair
                i = findfirst(x -> isequal(x, atom), unique_atoms)
                if i === nothing
                    # Add the atom to the atoms array if it's not already present
                    push!(unique_atoms, atom)
                    n_atoms += 1
                    A = [A; zeros(Int, 1, ncompounds)]
                    i = n_atoms
                end

                # Adjust coefficient based on whether the compound is a product or substrate
                coeff *= coeffsign

                A[i, j] = coeff
            end
        end

        # update for iterating through products
        coeffsign = -1
        jbase = length(substrates)
    end

    return A
end

function get_balanced_stoich(reaction::Reaction)
    # Create the reaction matrix A that is m atoms by n compounds
    A = create_matrix(reaction)

    # get an integer nullspace basis
    X = ModelingToolkit.nullspace(A)
    nullity = size(X, 2)

    stoichvecs = Vector{Vector{Int64}}()
    for (j, col) in pairs(eachcol(X))
        signs = unique(col)

        # If there is only one basis vector and the signs are not all the same this means
        # we have a solution that would require moving at least one substrate to be a
        # product (or vice-versa). We therefore do not return anything in this case.
        # If there are multiple basis vectors we don't currently determine if we can
        # construct a linear combination giving a solution, so we just return them.
        if (nullity > 1) || (all(>=(0), signs) || all(<=(0), signs))
            coefs = abs.(col)
            common_divisor = reduce(gcd, coefs)
            coefs .= div.(coefs, common_divisor)
            push!(stoichvecs, coefs)
        end
    end

    return stoichvecs
end

function balance_reaction(reaction::Reaction)
    # Calculate the stoichiometric coefficients for the balanced reaction.
    stoichiometries = get_balanced_stoich(reaction)

    balancedrxs = Vector{Reaction}(undef, length(stoichiometries))

    # Iterate over each stoichiometry vector and create a reaction
    for (i,stoich) in enumerate(stoichiometries)
        # Divide the stoichiometry vector into substrate and product stoichiometries.
        substoich = stoich[1:length(reaction.substrates)]
        prodstoich = stoich[(length(reaction.substrates) + 1):end]

        # Create a new reaction with the balanced stoichiometries
        balancedrx = Reaction(reaction.rate, reaction.substrates,
                              reaction.products, substoich, prodstoich)

        # Add the reaction to the vector of all reactions
        balancedrxs[i] = balancedrx
    end

    isempty(balancedrxs) && (@warn "Unable to balance reaction.")
    (length(balancedrxs) > 1) && (@warn "Infinite balanced reactions from ($reaction) are possible, returning a basis for them. Note that we do not check if they preserve the set of substrates and products from the original reaction.")
    return balancedrxs
end
