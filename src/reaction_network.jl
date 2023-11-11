### Temporary deprecation warning - Eventually to be removed. ###
deprication_message = """
@reaction_network notation where parameters are declared after "end", e.g. like:

```julia
@reaction_network begin
    p, 0 --> X
    d, X --> 0
end p d
```

has been deprecated in favor of a notation where the parameters are inferred, e.g:

```julia
@reaction_network begin
    p, 0 --> X
    d, X --> 0
end
```

Parameters and species can be explicitly indicated using the @parameters and @species
macros, e.g:

```julia
@reaction_network begin
    @parameters p d
    @species X(t)
    p, 0 --> X
    d, X --> 0
end
```
"""
macro reaction_network(name::Symbol, ex::Expr, parameters...)
    error(deprication_message)
end
macro reaction_network(name::Expr, ex::Expr, parameters...)
    error(deprication_message)
end
macro reaction_network(ex::Expr, parameters...)
    error(deprication_message)
end

"""
Macro that inputs an expression corresponding to a reaction network and outputs
a `ReactionNetwork` that can be used as input to generation of ODE, SDE, and
Jump problems.

Most arrows accepted (both right, left, and bi-drectional arrows). Note that
while --> is a correct arrow, neither <-- nor <--> works. Using non-filled
arrows (⇐, ⟽, ⇒, ⟾, ⇔, ⟺) will disable mass kinetics and let you cutomize
reaction rates yourself. Use 0 or ∅ for degradation/creation to/from nothing.

Example systems:

    ### Basic Usage ###
    rn = @reaction_network begin           # Creates a ReactionSystem.
        2.0, X + Y --> XY                  # This will have reaction rate corresponding to 2.0*[X][Y]
        2.0, XY ← X + Y                    # Identical to 2.0, X + Y --> XY
    end

    ### Manipulating Reaction Rates ###
    rn = @reaction_network begin
        2.0, X + Y ⟾ XY                   # Ignores mass kinetics. This will have reaction rate corresponding to 2.0.
        2.0X, X + Y --> XY                 # Reaction rate needs not be constant. This will have reaction rate corresponding to 2.0*[X]*[X]*[Y].
        XY+log(X)^2, X + Y --> XY          # Reaction rate accepts quite complicated expressions.
        hill(XY,2,2,2), X + Y --> XY       # Reaction inis activated by XY according to a hill function. hill(x,v,K,N).
        mm(XY,2,2), X + Y --> XY           # Reaction inis activated by XY according to a michaelis menten function. mm(x,v,K).
    end

    ### Multiple Reactions on a Single Line ###
    rn = @reaction_network begin
        (2.0,1.0), X + Y ↔ XY              # Identical to reactions (2.0, X + Y --> XY) and (1.0, XY --> X + Y).
        2.0, (X,Y) --> 0                   # This corresponds to both X and Y degrading at rate 2.0.
        (2.0, 1.0), (X,Y) --> 0            # This corresponds to X and Y degrading at rates 2.0 and 1.0, respectively.
        2.0, (X1,Y1) --> (X2,Y2)           # X1 and Y1 becomes X2 and Y2, respectively, at rate 2.0.
    end

    ### Adding Parameters ###
    kB = 2.0; kD = 1.0
    p = [kB, kD]
    p = []
    rn = @reaction_network begin
        (kB, kD), X + Y ↔ XY               # Lets you define parameters outside on network. Parameters can be changed without recalling the network.
    end

    ### Defining New Functions ###
    my_hill_repression(x, v, k, n) = v*k^n/(k^n+x^n)

    # may be necessary to
    # @register_symbolic my_hill_repression(x, v, k, n)
    # see https://docs.sciml.ai/ModelingToolkit/stable/basics/Validation/#User-Defined-Registered-Functions-and-Types

    r = @reaction_network MyReactionType begin
        my_hill_repression(x, v_x, k_x, n_x), 0 --> x
    end

    ### Simulating Reaction Networks ###
    probODE = ODEProblem(rn, args...; kwargs...)        # Using multiple dispatch the reaction network can be used as input to create ODE, SDE and Jump problems.
    probSDE = SDEProblem(rn, args...; kwargs...)
    probJump = JumpProblem(prob,aggregator::Direct,rn)
"""

### Declares various options and constants. ###

# Declare various arrow types symbols used for the empty set (also 0).
const empty_set = Set{Symbol}([:∅])
const fwd_arrows = Set{Symbol}([:>, :(=>), :→, :↣, :↦, :⇾, :⟶, :⟼, :⥟, :⥟, :⇀, :⇁, :⇒, :⟾])
const bwd_arrows = Set{Symbol}([:<, :(<=), :←, :↢, :↤, :⇽, :⟵, :⟻, :⥚, :⥞, :↼, :↽, :⇐, :⟽,
                                   Symbol("<--")])
const double_arrows = Set{Symbol}([:↔, :⟷, :⇄, :⇆, :⇌, :⇋, :⇔, :⟺, Symbol("<-->")])
const pure_rate_arrows = Set{Symbol}([:(=>), :(<=), :⇐, :⟽, :⇒, :⟾, :⇔, :⟺])

const CONSERVED_CONSTANT_SYMBOL = :Γ

# Declares symbols which may neither be used as parameters not varriables.
const forbidden_symbols_skip = Set([:ℯ, :pi, :π, :t, :∅])
const forbidden_symbols_error = union(Set([:im, :nothing, CONSERVED_CONSTANT_SYMBOL]),
                                      forbidden_symbols_skip)
const forbidden_variables_error = let
    fvars = copy(forbidden_symbols_error)
    delete!(fvars, :t)
    fvars
end

# Declares the keys used for various options.
const option_keys = (:species, :parameters, :variables, :ivs, :compounds)

### The @species macro, basically a copy of the @variables macro. ###
macro species(ex...)
    vars = Symbolics._parse_vars(:variables, Real, ex)

    # vector of symbols that get defined
    lastarg = vars.args[end]

    # start adding metadata statements where the vector of symbols was previously declared
    idx = length(vars.args)
    resize!(vars.args, idx + length(lastarg.args) + 1)
    for sym in lastarg.args
        vars.args[idx] = :($sym = ModelingToolkit.wrap(setmetadata(ModelingToolkit.value($sym),
                                                                   Catalyst.VariableSpecies,
                                                                   true)))
        idx += 1
    end

    # check nothing was declared isconstantspecies
    ex = quote
        all(!Catalyst.isconstant ∘ ModelingToolkit.value, $lastarg) ||
            throw(ArgumentError("isconstantspecies metadata can only be used with parameters."))
    end
    vars.args[idx] = ex
    idx += 1

    # put back the vector of the new species symbols
    vars.args[idx] = lastarg

    esc(vars)
end

### The main macro, takes reaction network notation and returns a ReactionSystem. ###
"""
    @reaction_network

Generates a [`ReactionSystem`](@ref dsl_description) that encodes a chemical reaction
network.

See [The Reaction DSL](@ref dsl_description) documentation for details on
parameters to the macro.

Examples:
```julia
# a basic SIR model, with name SIR
sir_model = @reaction_network SIR begin
    c1, s + i --> 2i
    c2, i --> r
end

# a basic SIR model, with random generated name
sir_model = @reaction_network begin
    c1, s + i --> 2i
    c2, i --> r
end

# an empty network with name empty
emptyrn = @reaction_network empty

# an empty network with random generated name
emptyrn = @reaction_network
```
"""
macro reaction_network(name::Symbol, ex::Expr)
    make_reaction_system(MacroTools.striplines(ex); name = :($(QuoteNode(name))))
end

# allows @reaction_network $name begin ... to interpolate variables storing a name
macro reaction_network(name::Expr, ex::Expr)
    make_reaction_system(MacroTools.striplines(ex); name = :($(esc(name.args[1]))))
end

macro reaction_network(ex::Expr)
    ex = MacroTools.striplines(ex)

    # no name but equations: @reaction_network begin ... end ...
    if ex.head == :block
        make_reaction_system(ex)
    else  # empty but has interpolated name: @reaction_network $name
        networkname = :($(esc(ex.args[1])))
        return Expr(:block, :(@parameters t),
                    :(ReactionSystem(Reaction[], t, [], []; name = $networkname)))
    end
end

#Returns a empty network (with, or without, a declared name)
# @reaction_network name
macro reaction_network(name::Symbol = gensym(:ReactionSystem))
    return Expr(:block, :(@parameters t),
                :(ReactionSystem(Reaction[], t, [], []; name = $(QuoteNode(name)))))
end

### Macros used for manipulating, and successively builing up, reaction systems. ###
@doc raw"""
    @reaction

Generates a single [`Reaction`](@ref) object.

Examples:
```julia
rx = @reaction k*v, A + B --> C + D

# is equivalent to
@parameters k v
@variables t
@species A(t) B(t) C(t) D(t)
rx == Reaction(k*v, [A,B], [C,D])
```
Here `k` and `v` will be parameters and `A`, `B`, `C` and `D` will be variables.
Interpolation of existing parameters/variables also works
```julia
@parameters k b
@variables t
@species A(t)
ex = k*A^2 + t
rx = @reaction b*$ex*$A, $A --> C
```

Notes:
- Any symbols arising in the rate expression that aren't interpolated are treated as
  parameters. In the reaction part (`α*A + B --> C + D`), coefficients are treated as
  parameters, e.g. `α`, and rightmost symbols as species, e.g. `A,B,C,D`.
- Works with any *single* arrow types supported by [`@reaction_network`](@ref).
- Interpolation of Julia variables into the macro works similar to the `@reaction_network`
  macro. See [The Reaction DSL](@ref dsl_description) tutorial for more details.
"""
macro reaction(ex)
    make_reaction(ex)
end

"""
    @add_reactions

Adds the reactions declared to a preexisting [`ReactionSystem`](@ref). Note, mutates the
original network.

Notes:
- To instead generate a new network by combining two existing networks use
  `ModelingToolkit.extend`.

Example:
```julia
rn = @reaction_network begin
    @parameters G
    π, 2*A --> B
    end

# add this reaction into rn
@add_reactions rn begin
    k*A, C --> D
end
```
"""
macro add_reactions(rn::Symbol, ex::Expr)
    :(merge!($(esc(rn)), $(make_reaction_system(MacroTools.striplines(ex)))))
end

### Internal DSL structures for representing reactants and reactions. ###

#Structure containing information about one reactant in one reaction.
struct ReactantStruct
    reactant::Union{Symbol, Expr}
    stoichiometry::ExprValues
end
#Structure containing information about one Reaction. Contain all its substrates and products as well as its rate. Contains a specialized constructor.
struct ReactionStruct
    substrates::Vector{ReactantStruct}
    products::Vector{ReactantStruct}
    rate::ExprValues
    only_use_rate::Bool

    function ReactionStruct(sub_line::ExprValues, prod_line::ExprValues, rate::ExprValues,
                            only_use_rate::Bool)
        sub = recursive_find_reactants!(sub_line, 1, Vector{ReactantStruct}(undef, 0))
        prod = recursive_find_reactants!(prod_line, 1, Vector{ReactantStruct}(undef, 0))
        new(sub, prod, rate, only_use_rate)
    end
end

### Functions rephrasing the macro input as a ReactionSystem structure. ###

function forbidden_variable_check(v)
    !isempty(intersect(forbidden_variables_error, v)) &&
        error("The following symbol(s) are used as variables: " *
              ((map(s -> "'" * string(s) * "', ",
                    intersect(forbidden_variables_error, v))...)) *
              "this is not permited.")
end

function forbidden_symbol_check(v)
    !isempty(intersect(forbidden_symbols_error, v)) &&
        error("The following symbol(s) are used as species or parameters: " *
              ((map(s -> "'" * string(s) * "', ",
                    intersect(forbidden_symbols_error, v))...)) *
              "this is not permited.")
    nothing
end

function unique_symbol_check(syms)
    allunique(syms) ||
        error("Reaction network independent variables, parameters, species, and variables must all have distinct names, but a duplicate has been detected. ")
    nothing
end

# takes a ModelingToolkit declaration macro like @parameters and returns an expression
# that calls the macro and then scalarizes all the symbols created into a vector of Nums
function scalarize_macro(nonempty, ex, name)
    namesym = gensym(name)
    if nonempty
        symvec = gensym()
        ex = quote
            $symvec = $ex
            $namesym = reduce(vcat, Symbolics.scalarize($symvec))
        end
    else
        ex = :($namesym = Num[])
    end
    ex, namesym
end

# Function for creating a ReactionSystem structure (used by the @reaction_network macro).
function make_reaction_system(ex::Expr; name = :(gensym(:ReactionSystem)))

    # Handle interpolation of variables
    ex = esc_dollars!(ex)

    # Read lines with reactions and options.
    reaction_lines = Expr[x for x in ex.args if x.head == :tuple]
    option_lines = Expr[x for x in ex.args if x.head == :macrocall]

    # Get macro options.
    options = Dict(map(arg -> Symbol(String(arg.args[1])[2:end]) => arg,
                       option_lines))

    # Reads compounds options.
    compound_expr, compound_species = read_compound_options(options)

    # Parses reactions, species, and parameters.
    reactions = get_reactions(reaction_lines)
    species_declared = [extract_syms(options, :species); compound_species]
    parameters_declared = extract_syms(options, :parameters)
    variables = extract_syms(options, :variables)

    # handle independent variables
    if haskey(options, :ivs)
        ivs = Tuple(extract_syms(options, :ivs))
        ivexpr = copy(options[:ivs])
        ivexpr.args[1] = Symbol("@", "variables")
    else
        ivs = (DEFAULT_IV_SYM,)
        ivexpr = :(@variables $(DEFAULT_IV_SYM))
    end
    tiv = ivs[1]
    sivs = (length(ivs) > 1) ? Expr(:vect, ivs[2:end]...) : nothing

    declared_syms = Set(Iterators.flatten((parameters_declared, species_declared,
                                           variables)))
    species_extracted, parameters_extracted = extract_species_and_parameters!(reactions,
                                                                              declared_syms)
    species = vcat(species_declared, species_extracted)
    parameters = vcat(parameters_declared, parameters_extracted)

    # Checks for input errors.
    (sum(length.([reaction_lines, option_lines])) != length(ex.args)) &&
        error("@reaction_network input contain $(length(ex.args) - sum(length.([reaction_lines,option_lines]))) malformed lines.")
    any(!in(opt_in, option_keys) for opt_in in keys(options)) &&
        error("The following unsupported options were used: $(filter(opt_in->!in(opt_in,option_keys), keys(options)))")
    forbidden_symbol_check(union(species, parameters))
    forbidden_variable_check(variables)
    unique_symbol_check(union(species, parameters, variables, ivs))

    # Creates expressions corresponding to actual code from the internal DSL representation.
    sexprs = get_sexpr(species_extracted, options; iv_symbols = ivs)
    vexprs = haskey(options, :variables) ? options[:variables] : :()
    pexprs = get_pexpr(parameters_extracted, options)
    ps, pssym = scalarize_macro(!isempty(parameters), pexprs, "ps")
    vars, varssym = scalarize_macro(!isempty(variables), vexprs, "vars")
    sps, spssym = scalarize_macro(!isempty(species), sexprs, "specs")
    comps, compssym = scalarize_macro(!isempty(compound_species), compound_expr, "comps")
    rxexprs = :(Catalyst.CatalystEqType[])
    for reaction in reactions
        push!(rxexprs.args, get_rxexprs(reaction))
    end

    # Returns the rephrased expression.
    quote
        $ps
        $ivexpr
        $vars
        $sps
        $comps

        Catalyst.make_ReactionSystem_internal($rxexprs, $tiv, union($spssym, $varssym, $compssym),
                                              $pssym; name = $name,
                                              spatial_ivs = $sivs)
    end
end

# Function for creating a Reaction structure (used by the @reaction macro).
function make_reaction(ex::Expr)

    # Handle interpolation of variables
    ex = esc_dollars!(ex)

    # Parses reactions, species, and parameters.
    reaction = get_reaction(ex)
    species, parameters = extract_species_and_parameters!([reaction], [])

    # Checks for input errors.
    forbidden_symbol_check(union(species, parameters))

    # Creates expressions corresponding to actual code from the internal DSL representation.
    sexprs = get_sexpr(species, Dict{Symbol, Expr}())
    pexprs = get_pexpr(parameters, Dict{Symbol, Expr}())
    rxexpr = get_rxexprs(reaction)
    iv = :(@variables $(DEFAULT_IV_SYM))

    # Returns the rephrased expression.
    quote
        $pexprs
        $iv
        $sexprs
        $rxexpr
    end
end

### Auxiliary function called by the main make_reaction_system and make_reaction functions. ###

# Function that handles variable interpolation.
function esc_dollars!(ex)
    if ex isa Expr
        if ex.head == :$
            return esc(:($(ex.args[1])))
        else
            for i in 1:length(ex.args)
                ex.args[i] = esc_dollars!(ex.args[i])
            end
        end
    end
    ex
end

# When compound species are declared using the "@compound begin ... end" option, get a list of the compound species, and also the expression that crates them.
function read_compound_options(opts)
    # If the compound option is used retrive a list of compound species, and the option that creeates them
    if haskey(opts, :compounds)
        compound_expr = opts[:compounds]
        compound_species = [arg.args[1].args[1] for arg in compound_expr.args[3].args] # Loops through where in the "@compound begin ... end" the compound species names are.
    else
        compound_expr = :()
        compound_species = Union{Symbol, Expr}[]
    end
    return compound_expr, compound_species
end

# When the user have used the @species (or @parameters) options, extract species (or
# parameters) from its input.
function extract_syms(opts, vartype::Symbol)
    if haskey(opts, vartype)
        ex = opts[vartype]
        vars = Symbolics._parse_vars(vartype, Real, ex.args[3:end])
        syms = Vector{Union{Symbol, Expr}}(vars.args[end].args)
    else
        syms = Union{Symbol, Expr}[]
    end
    return syms
end

# Function looping through all reactions, to find undeclared symbols (species or
# parameters), and assign them to the right category.
function extract_species_and_parameters!(reactions, excluded_syms)
    species = OrderedSet{Union{Symbol, Expr}}()
    for reaction in reactions
        for reactant in Iterators.flatten((reaction.substrates, reaction.products))
            add_syms_from_expr!(species, reactant.reactant, excluded_syms)
        end
    end

    foreach(s -> push!(excluded_syms, s), species)
    parameters = OrderedSet{Union{Symbol, Expr}}()
    for reaction in reactions
        add_syms_from_expr!(parameters, reaction.rate, excluded_syms)
        for reactant in Iterators.flatten((reaction.substrates, reaction.products))
            add_syms_from_expr!(parameters, reactant.stoichiometry, excluded_syms)
        end
    end

    collect(species), collect(parameters)
end
# Function called by extract_species_and_parameters!, recursively loops through an
# expression and find symbols (adding them to the push_symbols vector).
function add_syms_from_expr!(push_symbols::AbstractSet, rateex::ExprValues, excluded_syms)
    if rateex isa Symbol
        if !(rateex in forbidden_symbols_skip) && !(rateex in excluded_syms)
            push!(push_symbols, rateex)
        end
    elseif rateex isa Expr
        # note, this (correctly) skips $(...) expressions
        for i in 2:length(rateex.args)
            add_syms_from_expr!(push_symbols, rateex.args[i], excluded_syms)
        end
    end
    nothing
end

# Given the species that were extracted from the reactions, and the options dictionary, creates the @species ... expression for the macro output.
function get_sexpr(species_extracted, options, key = :species;
                   iv_symbols = (DEFAULT_IV_SYM,))
    if haskey(options, key)
        sexprs = options[key]
    elseif isempty(species_extracted)
        sexprs = :()
    else
        sexprs = Expr(:macrocall, Symbol("@", key), LineNumberNode(0))
    end
    foreach(s -> (s isa Symbol) && push!(sexprs.args, Expr(:call, s, iv_symbols...)),
            species_extracted)
    sexprs
end

# Given the parameters that were extracted from the reactions, and the options dictionary, creates the @parameters ... expression for the macro output.
function get_pexpr(parameters_extracted, options)
    pexprs = (haskey(options, :parameters) ? options[:parameters] :
              (isempty(parameters_extracted) ? :() : :(@parameters)))
    foreach(p -> push!(pexprs.args, p), parameters_extracted)
    pexprs
end

# Creates the reaction vector declaration statement.
function get_rxexprs(rxstruct)
    subs_init = isempty(rxstruct.substrates) ? nothing : :([])
    subs_stoich_init = deepcopy(subs_init)
    prod_init = isempty(rxstruct.products) ? nothing : :([])
    prod_stoich_init = deepcopy(prod_init)
    reaction_func = :(Reaction($(recursive_expand_functions!(rxstruct.rate)), $subs_init,
                               $prod_init, $subs_stoich_init, $prod_stoich_init,
                               only_use_rate = $(rxstruct.only_use_rate)))
    for sub in rxstruct.substrates
        push!(reaction_func.args[3].args, sub.reactant)
        push!(reaction_func.args[5].args, sub.stoichiometry)
    end
    for prod in rxstruct.products
        push!(reaction_func.args[4].args, prod.reactant)
        push!(reaction_func.args[6].args, prod.stoichiometry)
    end
    reaction_func
end

### Functions for extracting the reactions from a DSL expression, and putting them ReactionStruct vector. ###

# Reads a line and creates the corresponding ReactionStruct.
function get_reaction(line)
    (rate, r_line) = line.args
    (r_line.head == :-->) && (r_line = Expr(:call, :→, r_line.args[1], r_line.args[2]))

    arrow = r_line.args[1]
    in(arrow, double_arrows) && error("Double arrows not allowed for single reactions.")

    only_use_rate = in(arrow, pure_rate_arrows)
    if in(arrow, fwd_arrows)
        rs = create_ReactionStruct(r_line.args[2], r_line.args[3], rate, only_use_rate)
    elseif in(arrow, bwd_arrows)
        rs = create_ReactionStruct(r_line.args[3], r_line.args[2], rate, only_use_rate)
    else
        throw("Malformed reaction, invalid arrow type used in: $(MacroTools.striplines(line))")
    end

    rs
end
# Generates a vector containing a number of reaction structures, each containing the information about one reaction.
function get_reactions(exprs::Vector{Expr}, reactions = Vector{ReactionStruct}(undef, 0))
    for line in exprs
        (rate, r_line) = line.args
        (r_line.head == :-->) && (r_line = Expr(:call, :→, r_line.args[1], r_line.args[2]))

        arrow = r_line.args[1]
        only_use_rate = in(arrow, pure_rate_arrows)
        if in(arrow, double_arrows)
            (typeof(rate) == Expr && rate.head == :tuple) ||
                error("Error: Must provide a tuple of reaction rates when declaring a bi-directional reaction.")
            push_reactions!(reactions, r_line.args[2], r_line.args[3], rate.args[1],
                            only_use_rate)
            push_reactions!(reactions, r_line.args[3], r_line.args[2], rate.args[2],
                            only_use_rate)
        elseif in(arrow, fwd_arrows)
            push_reactions!(reactions, r_line.args[2], r_line.args[3], rate, only_use_rate)
        elseif in(arrow, bwd_arrows)
            push_reactions!(reactions, r_line.args[3], r_line.args[2], rate, only_use_rate)
        else
            throw("Malformed reaction, invalid arrow type used in: $(MacroTools.striplines(line))")
        end
    end
    reactions
end

# Creates a ReactionStruct from the information in a single line.
function create_ReactionStruct(sub_line::ExprValues, prod_line::ExprValues,
                               rate::ExprValues, only_use_rate::Bool)
    all(==(1), (tup_leng(sub_line), tup_leng(prod_line), tup_leng(rate))) ||
        error("Malformed reaction, line appears to be defining multiple reactions incorrectly: rate=$rate, subs=$sub_line, prods=$prod_line.")
    ReactionStruct(get_tup_arg(sub_line, 1), get_tup_arg(prod_line, 1),
                   get_tup_arg(rate, 1), only_use_rate)
end

#Takes a reaction line and creates reactions from it and pushes those to the reaction array. Used to create multiple reactions from, for instance, 1.0, (X,Y) --> 0.
function push_reactions!(reactions::Vector{ReactionStruct}, sub_line::ExprValues,
                         prod_line::ExprValues, rate::ExprValues, only_use_rate::Bool)
    lengs = (tup_leng(sub_line), tup_leng(prod_line), tup_leng(rate))
    for i in 1:maximum(lengs)
        (count(lengs .== 1) + count(lengs .== maximum(lengs)) < 3) &&
            (throw("Malformed reaction, rate=$rate, subs=$sub_line, prods=$prod_line."))
        push!(reactions,
              ReactionStruct(get_tup_arg(sub_line, i), get_tup_arg(prod_line, i),
                             get_tup_arg(rate, i), only_use_rate))
    end
end

function processmult(op, mult, stoich)
    if (mult isa Number) && (stoich isa Number)
        op(mult, stoich)
    else
        :($op($mult, $stoich))
    end
end

# Recursive function that loops through the reaction line and finds the reactants and their
# stoichiometry. Recursion makes it able to handle weird cases like 2(X+Y+3(Z+XY)).
function recursive_find_reactants!(ex::ExprValues, mult::ExprValues,
                                   reactants::Vector{ReactantStruct})
    if typeof(ex) != Expr || (ex.head == :escape) || (ex.head == :ref)
        (ex == 0 || in(ex, empty_set)) && (return reactants)
        if any(ex == reactant.reactant for reactant in reactants)
            idx = findall(x -> x == ex, getfield.(reactants, :reactant))[1]
            reactants[idx] = ReactantStruct(ex,
                                            processmult(+, mult,
                                                        reactants[idx].stoichiometry))
        else
            push!(reactants, ReactantStruct(ex, mult))
        end
    elseif ex.args[1] == :*
        if length(ex.args) == 3
            recursive_find_reactants!(ex.args[3], processmult(*, mult, ex.args[2]),
                                      reactants)
        else
            newmult = processmult(*, mult, Expr(:call, ex.args[1:(end - 1)]...))
            recursive_find_reactants!(ex.args[end], newmult, reactants)
        end
    elseif ex.args[1] == :+
        for i in 2:length(ex.args)
            recursive_find_reactants!(ex.args[i], mult, reactants)
        end
    else
        throw("Malformed reaction, bad operator: $(ex.args[1]) found in stochiometry expression $ex.")
    end
    reactants
end

### Functionality for expanding function call to custom and specific functions ###

#Recursively traverses an expression and replaces special function call like "hill(...)" with the actual corresponding expression.
function recursive_expand_functions!(expr::ExprValues)
    (typeof(expr) != Expr) && (return expr)
    foreach(i -> expr.args[i] = recursive_expand_functions!(expr.args[i]),
            1:length(expr.args))
    if expr.head == :call
        !isdefined(Catalyst, expr.args[1]) && (expr.args[1] = esc(expr.args[1]))
    end
    expr
end

# ### Old functions (for deleting).

# function get_rx_species_deletethis(rxs, ps)
#     pset = Set(ps)
#     species_set = Set{Symbol}()
#     for rx in rxs
#         find_species_in_rate!(species_set, rx.rate, pset)
#         for sub in rx.substrates
#             find_species_in_rate!(species_set, sub.stoichiometry, pset)
#         end
#         for prod in rx.products
#             find_species_in_rate!(species_set, prod.stoichiometry, pset)
#         end
#     end
#     collect(species_set)
# end

# function find_species_in_rate!_deletethis(sset, rateex::ExprValues, ps)
#     if rateex isa Symbol
#         if !(rateex in forbidden_symbols) && !(rateex in ps)
#             push!(sset, rateex)
#         end
#     elseif rateex isa Expr
#         # note, this (correctly) skips $(...) expressions
#         for i in 2:length(rateex.args)
#             find_species_in_rate!(sset, rateex.args[i], ps)
#         end
#     end
#     nothing
# end

# function get_reactants_deletethis(reaction::ReactionStruct,
#                                   reactants = Vector{Union{Symbol, Expr}}())
#     for reactant in Iterators.flatten((reaction.substrates, reaction.products))
#         !in(reactant.reactant, reactants) && push!(reactants, reactant.reactant)
#     end
#     return reactants
# end

# # Extract the reactants from the set of reactions.
# function get_reactants_deletethis(reactions::Vector{ReactionStruct},
#                                   reactants = Vector{Union{Symbol, Expr}}())
#     for reaction in reactions
#         get_reactants(reaction, reactants)
#     end
#     return reactants
# end

# # Gets the species/parameter symbols from the reactions (when the user has omitted the designation of these).
# function extract_species(reactions::Vector{ReactionStruct}, parameters::Vector,
#                          species = Vector{Union{Symbol, Expr}}())
#     for reaction in reactions,
#         reactant in Iterators.flatten((reaction.substrates, reaction.products))

#         find_parameters_in_expr!(species, reaction.rate, parameters)
#         !in(reactant.reactant, species) && !in(reactant.reactant, parameters) &&
#             push!(species, reactant.reactant)
#     end
#     return species
# end
# function extract_parameters(reactions::Vector{ReactionStruct},
#                             species::Vector{Union{Symbol, Expr}},
#                             parameters = Vector{Symbol}())
#     for rx in reactions
#         find_parameters_in_expr!(parameters, rx.rate, species)
#         for sub in rx.substrates
#             find_parameters_in_expr!(parameters, sub.stoichiometry, species)
#         end
#         for prod in rx.products
#             find_parameters_in_expr!(parameters, prod.stoichiometry, species)
#         end
#     end
#     return parameters
# end

# # Goes through an expression, and returns the paramters in it.
# function find_parameters_in_expr!(parameters, rateex::ExprValues,
#                                   species::Vector)
#     if rateex isa Symbol
#         if !(rateex in forbidden_symbols) && !(rateex in species)
#             push!(parameters, rateex)
#         end
#     elseif rateex isa Expr
#         # note, this (correctly) skips $(...) expressions
#         for i in 2:length(rateex.args)
#             find_parameters_in_expr!(parameters, rateex.args[i], species)
#         end
#     end
#     nothing
# end

# # Loops through the users species and parameter inputs, and checks if any have default values.
# function make_default_args(options)
#     defaults = :(Dict([]))
#     haskey(options, :species) && for arg in options[:species]
#         (arg isa Symbol) && continue
#         (arg.head != :(=)) && continue
#         push!(defaults.args[2].args, :($(arg.args[1]) => $(arg.args[2])))
#     end
#     haskey(options, :parameters) && for arg in options[:parameters]
#         (arg isa Symbol) && continue
#         (arg.head != :(=)) && continue
#         push!(defaults.args[2].args, :($(arg.args[1]) => $(arg.args[2])))
#     end
#     return defaults
# end
