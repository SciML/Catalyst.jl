### Declares Compound Related Metadata ###
struct CompoundSpecies end
struct CompoundComponents end
struct CompoundCoefficients end

Symbolics.option_to_metadata_type(::Val{:iscompound}) = CompoundSpecies
Symbolics.option_to_metadata_type(::Val{:components}) = CompoundComponents
Symbolics.option_to_metadata_type(::Val{:coefficients}) = CompoundCoefficients

### Create @compound Macro(s) ###

"""
    @compound

Macro that creates a compound species, which is composed of smaller component species.

Example:
```julia
@variables t
@species C(t) O(t)
@compound CO2(t) ~ C + 2O
```

Notes: 
- The component species must be defined before using the `@compound` macro.
"""
macro compound(expr)
    make_compound(MacroTools.striplines(expr))
end

# Declares compound error messages:
const COMPOUND_CREATION_ERROR_BASE = "Malformed input to @compound. Should use form like e.g. \"@compound CO2 ~ C + 2O\"."
const COMPOUND_CREATION_ERROR_BAD_SEPARATOR = "Malformed input to @compound. Left-hand side (the compound) and the right-hand side (the components) should be separated by a \"~\" (e.g. \"@compound CO2 ~ C + 2O\"). If the left hand side contains metadata and/or default values, this should be enclosed by \"()\" (e.g. \"@compound (CO2 = 1.0, [output=true]) ~ C + 2O\")."
const COMPOUND_CREATION_ERROR_DEPENDENT_VAR_REQUIRED = "When the components (collectively) have more than 1 independent variable, independent variables have to be specified for the compound (e.g. `@compound CO2(t,s) ~ C + 2O`). This is required since the @compound macro cannot infer the correct order of the independent variables."

# Function managing the @compound macro.
function make_compound(expr)
    # Error checks.
    (expr isa Expr) || error(COMPOUND_CREATION_ERROR_BASE)
    ((expr.head == :call) && (expr.args[1] == :~) && (length(expr.args) == 3)) || error(COMPOUND_CREATION_ERROR_BAD_SEPARATOR)

    # Loops through all components, add the component and the coefficients to the corresponding vectors
    # Cannot extract directly using e.g. "getfield.(composition, :reactant)" because then 
    # we get something like :([:C, :O]), rather than :([C, O]).
    composition = Catalyst.recursive_find_reactants!(expr.args[3], 1, Vector{ReactantStruct}(undef, 0))
    components = :([])                                      # Becomes something like :([C, O]).                                         
    coefficients = :([])                                    # Becomes something like :([1, 2]). 
    for comp in composition
        push!(components.args, comp.reactant)
        push!(coefficients.args, comp.stoichiometry)
    end

    # Extracts:
    # - The compound species name (species_name, e.g. `:CO2`).
    # - Any ivs attached to it (ivs, e.g. `[]` or `[t,x]`).
    # - The expression which creates the compound (species_expr, e.g. `CO2 = 1.0, [metadata=true]`).
    species_expr = expr.args[2]
    species_name, ivs, _, _ = find_varinfo_in_declaration(expr.args[2])

    # If no ivs were given, inserts `(..)` (e.g. turning `CO` to `CO(..)`).
    isempty(ivs) && (species_expr = insert_independent_variable(species_expr, :(..)))      

    # Expression which when evaluated gives a vector with all the ivs of the components.
    ivs_get_expr = :(unique(reduce(vcat,[arguments(ModelingToolkit.unwrap(comp)) for comp in $components])))    
    
    # Creates the found expressions that will create the compound species.
    # The `Expr(:escape, :(...))` is required so that the expressions are evaluated in 
    # the scope the users use the macro in (to e.g. detect already exiting species).     
    # Creates something like (where `compound_ivs` and `component_ivs` evaluates to all the compound's and components' ivs):
    #   `@species CO2(..)`
    #   `isempty([])` && length(component_ivs) && error("When ...)
    #   `CO2 = CO2(component_ivs..)` 
    #   `issetequal(compound_ivs, component_ivs) || error("The ...)` 
    #   `CO2 = ModelingToolkit.setmetadata(CO2, Catalyst.CompoundSpecies, true)`
    #   `CO2 = ModelingToolkit.setmetadata(CO2, Catalyst.CompoundSpecies, [C, O])`
    #   `CO2 = ModelingToolkit.setmetadata(CO2, Catalyst.CompoundSpecies, [1, 2])`
    species_declaration_expr = Expr(:escape, :(@species $species_expr)) 
    multiple_ivs_error_check_expr = Expr(:escape, :($(isempty(ivs)) && (length($ivs_get_expr) > 1) && error($COMPOUND_CREATION_ERROR_DEPENDENT_VAR_REQUIRED)))      
    iv_designation_expr = Expr(:escape, :($(isempty(ivs)) && ($species_name = $(species_name)($(ivs_get_expr)...)))) 
    iv_check_expr = Expr(:escape, :(issetequal(arguments(ModelingToolkit.unwrap($species_name)), $ivs_get_expr) || error("The independent variable(S) provided to the compound ($(arguments(ModelingToolkit.unwrap($species_name)))), and those of its components ($($ivs_get_expr)))), are not identical.")))
    compound_designation_expr = Expr(:escape, :($species_name = ModelingToolkit.setmetadata($species_name, Catalyst.CompoundSpecies, true)))
    components_designation_expr = Expr(:escape, :($species_name = ModelingToolkit.setmetadata($species_name, Catalyst.CompoundComponents, $components)))
    coefficients_designation_expr = Expr(:escape, :($species_name = ModelingToolkit.setmetadata($species_name, Catalyst.CompoundCoefficients, $coefficients))) 

    # Currently, non-t independent variables are not supported for compounds. If there are any like these, we throw an error:
    non_t_iv_error_check_expr = Expr(:escape, :(issetequal(unique(reduce(vcat, arguments.(ModelingToolkit.unwrap.($components)))), [t]) || error("Currently, compounds depending on components that are not \"t\" are not supported.")))

println(non_t_iv_error_check_expr)

    # Returns the rephrased expression.
    return quote
        $non_t_iv_error_check_expr
        $species_declaration_expr
        $multiple_ivs_error_check_expr
        $iv_designation_expr
        $iv_check_expr
        $compound_designation_expr
        $components_designation_expr
        $coefficients_designation_expr
    end
end

"""
    @compounds

Macro that creates several compound species, which each is composed of smaller component species. Uses the same syntax as `@compound`, but with one compound species one each line.

Example:
```julia
@variables t
@species C(t) H(t) O(t) 
@compounds
    CH4(t) = C + 4H
    O2(t) = 2O
    CO2(t) = C + 2O
    H2O(t) = 2H + O
end
```

Notes: 
- The component species must be defined before using the `@compound` macro.
"""
macro compounds(expr)
    make_compounds(MacroTools.striplines(expr))
end

# Function managing the @compound macro.
function make_compounds(expr)
    # Creates an empty block containing the output call.
    compound_declarations = Expr(:block)

    # For each compound in `expr`, creates the set of 7 compound creation lines (using `make_compound`).
    # Next, loops through all 7*[Number of compounds] lines and add them to compound_declarations.
    compound_calls = [Catalyst.make_compound(line) for line in expr.args] 
    for compound_call in compound_calls, line in MacroTools.striplines(compound_call).args
        push!(compound_declarations.args, line)
    end

    # The output of the macros should be a vector with the compounds (same as for e.g. "@species begin ... end", also require for things to work in the DSL).
    # Creates an output vector, and loops through all compounds, adding them to it.
    push!(compound_declarations.args, :($(Expr(:escape, :([])))))
    compound_syms = :([])
    for arg in expr.args
        push!(compound_syms.args, find_varinfo_in_declaration(arg.args[2])[1])
    end
    push!(compound_declarations.args, :($(Expr(:escape, :($(compound_syms))))))

    # The output needs to be converted to Vector{Num} (from  Vector{SymbolicUtils.BasicSymbolic{Real}}) to be consistent with e.g. @variables.
    compound_declarations.args[end] = :([ModelingToolkit.wrap(cmp) for cmp in $(compound_declarations.args[end])])

    # Returns output that.
    return compound_declarations
end

## Compound Getters ###

"""
    iscompound(s)

Returns `true` if the input is a compound species (else false).
"""
iscompound(s::Num) = iscompound(MT.value(s))
function iscompound(s)
    MT.getmetadata(s, CompoundSpecies, false)
end

"""
    components(s)

Returns a vector with a list of all the components of a compound species (created using e.g. the @compound macro).
"""
components(s::Num) = components(MT.value(s))
function components(s)
    MT.getmetadata(s, CompoundComponents)
end

"""
    coefficients(s)

Returns a vector with a list of all the stoichiometric coefficients of the components of a compound species (created using e.g. the @compound macro).
"""
coefficients(s::Num) = coefficients(MT.value(s))
function coefficients(s)
    MT.getmetadata(s, CompoundCoefficients)
end

"""
    component_coefficients(s)

Returns a Vector{Pari{Symbol,Int64}}, listing a compounds species (created using e.g. the @compound macro) all the coefficients and their stoichiometric coefficients.
"""
component_coefficients(s::Num) = component_coefficients(MT.value(s))
function component_coefficients(s)
    return [c => co for (c, co) in zip(components(s), coefficients(s))]
end


### Reaction Balancing Functionality ###

### Reaction Balancing Functionality ###

# Reaction balancing error.
const COMPOUND_OF_COMPOUND_ERROR = ErrorException("Reaction balancing does not currently work for reactions involving compounds of compounds.")

# Note this does not correctly handle compounds of compounds currently.
# Internal function used by "balance_reaction" (via "get_balanced_stoich").
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

# Internal function used by "balance_reaction".
function get_balanced_stoich(reaction::Reaction)
    # Create the reaction matrix A that is m atoms by n compounds
    A = create_matrix(reaction)

    # get an integer nullspace basis
    X = ModelingToolkit.nullspace(A)
    nullity = size(X, 2)

    stoichvecs = Vector{Vector{Int64}}()
    for (j, col) in enumerate(eachcol(X))
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

"""
    balance_reaction(reaction::Reaction)

Returns a vector of all possible stoichiometrically balanced `Reaction` objects for the given `Reaction`.

Example:
```julia
@variables t
@species Si(t) Cl(t) H(t) O(t)
@compound SiCl4 ~ Si + 4Cl
@compound H2O ~ 2H + O
@compound H4SiO4 ~ 4H + Si + 4O
@compound HCl ~ H + Cl
rx = @reaction 1.0, SiCl4 + H2O --> H4SiO4 HCl
balance_reaction(rx) # Exactly one solution.
```

```julia
@variables t
@species C(t) H(t) O(t)
@compound CO ~ C + O
@compound CO2 ~ C + 2O
@compound H2 ~ 2H
@compound CH4 ~ C + 4H
@compound H2O ~ 2H + O
rx = @reaction 1.0, CO + CO2 + H2--> CH4 H2O
balance_reaction(rx) # Multiple solutions.
```

```julia
@variables t
@species Fe(t) S(t) O(t) H(t) N(t)
@compound FeS2 ~ Fe + 2S
@compound HNO3 ~ H + N + 3O
@compound Fe2S3O12 ~ 2Fe + 3S + 12O
@compound NO ~ N + O
@compound H2SO4 ~ 2H + S + 4O
rx = @reaction 1.0, FeS2 + HNO3 --> Fe2S3O12 NO + H2SO4
brxs = balance_reaction(rx) # No solution.
```

Notes:
- Balancing reactions that contain compounds of compounds is currently not supported.
- A reaction may not always yield a single solution; it could have an infinite number of solutions or none at all. When there are multiple solutions, a vector of all possible `Reaction` objects is returned. However, substrates and products may be interchanged as we currently do not solve for a linear combination that maintains the set of substrates and products. 
- If the reaction cannot be balanced, an empty `Reaction` vector is returned.
"""
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
