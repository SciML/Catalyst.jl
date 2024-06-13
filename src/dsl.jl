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

### Constant Declarations ###

# Declare various arrow types symbols used for the empty set (also 0).
const empty_set = Set{Symbol}([:∅])
const fwd_arrows = Set{Symbol}([:>, :(=>), :→, :↣, :↦, :⇾, :⟶, :⟼, :⥟, :⥟, :⇀, :⇁, :⇒, :⟾])
const bwd_arrows = Set{Symbol}([:<, :(<=), :←, :↢, :↤, :⇽, :⟵, :⟻, :⥚, :⥞, :↼, :↽, :⇐, :⟽,
    Symbol("<--")])
const double_arrows = Set{Symbol}([:↔, :⟷, :⇄, :⇆, :⇌, :⇋, :⇔, :⟺, Symbol("<-->")])
const pure_rate_arrows = Set{Symbol}([:(=>), :(<=), :⇐, :⟽, :⇒, :⟾, :⇔, :⟺])

# Declares the keys used for various options.
const option_keys = (:species, :parameters, :variables, :ivs, :compounds, :observables,
    :default_noise_scaling, :differentials, :equations,
    :continuous_events, :discrete_events, :combinatoric_ratelaws)

### `@species` Macro ###

# The @species macro, basically a copy of the @variables macro.
macro species(ex...)
    vars = Symbolics._parse_vars(:variables, Real, ex)

    # vector of symbols that get defined
    lastarg = vars.args[end]

    # start adding metadata statements where the vector of symbols was previously declared
    idx = length(vars.args)
    resize!(vars.args, idx + length(lastarg.args) + 1)
    for sym in lastarg.args
        vars.args[idx] = :($sym = ModelingToolkit.wrap(setmetadata(
            ModelingToolkit.value($sym), Catalyst.VariableSpecies, true)))
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

### `@reaction_network` and `@network_component` Macros ###

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

ReactionSystems generated through `@reaction_network` are complete.
"""
macro reaction_network(name::Symbol, ex::Expr)
    :(complete($(make_reaction_system(
        MacroTools.striplines(ex); name = :($(QuoteNode(name)))))))
end

# Allows @reaction_network $name begin ... to interpolate variables storing a name.
macro reaction_network(name::Expr, ex::Expr)
    :(complete($(make_reaction_system(
        MacroTools.striplines(ex); name = :($(esc(name.args[1])))))))
end

macro reaction_network(ex::Expr)
    ex = MacroTools.striplines(ex)

    # no name but equations: @reaction_network begin ... end ...
    if ex.head == :block
        :(complete($(make_reaction_system(ex))))
    else  # empty but has interpolated name: @reaction_network $name
        networkname = :($(esc(ex.args[1])))
        return Expr(:block, :(@parameters t),
            :(complete(ReactionSystem(Reaction[], t, [], []; name = $networkname))))
    end
end

# Returns a empty network (with, or without, a declared name).
macro reaction_network(name::Symbol = gensym(:ReactionSystem))
    return Expr(:block, :(@parameters t),
        :(complete(ReactionSystem(Reaction[], t, [], []; name = $(QuoteNode(name))))))
end

# Ideally, it would have been possible to combine the @reaction_network and @network_component macros.
# However, this issue: https://github.com/JuliaLang/julia/issues/37691 causes problem with interpolations
# if we make the @reaction_network macro call the @network_component macro.
"""
    @network_component

As @reaction_network, but the output system is not complete.
"""
macro network_component(name::Symbol, ex::Expr)
    make_reaction_system(MacroTools.striplines(ex); name = :($(QuoteNode(name))))
end

# Allows @network_component $name begin ... to interpolate variables storing a name.
macro network_component(name::Expr, ex::Expr)
    make_reaction_system(MacroTools.striplines(ex); name = :($(esc(name.args[1]))))
end

macro network_component(ex::Expr)
    ex = MacroTools.striplines(ex)

    # no name but equations: @network_component begin ... end ...
    if ex.head == :block
        make_reaction_system(ex)
    else  # empty but has interpolated name: @network_component $name
        networkname = :($(esc(ex.args[1])))
        return Expr(:block, :(@parameters t),
            :(ReactionSystem(Reaction[], t, [], []; name = $networkname)))
    end
end

# Returns a empty network (with, or without, a declared name).
macro network_component(name::Symbol = gensym(:ReactionSystem))
    return Expr(:block, :(@parameters t),
        :(ReactionSystem(Reaction[], t, [], []; name = $(QuoteNode(name)))))
end

### Internal DSL Structures ###

# Structure containing information about one reactant in one reaction.
struct ReactantStruct
    reactant::Union{Symbol, Expr}
    stoichiometry::ExprValues
end

# Structure containing information about one Reaction. Contain all its substrates and products as well as its rate. Contains a specialized constructor.
struct ReactionStruct
    substrates::Vector{ReactantStruct}
    products::Vector{ReactantStruct}
    rate::ExprValues
    metadata::Expr

    function ReactionStruct(sub_line::ExprValues, prod_line::ExprValues, rate::ExprValues,
            metadata_line::ExprValues)
        sub = recursive_find_reactants!(sub_line, 1, Vector{ReactantStruct}(undef, 0))
        prod = recursive_find_reactants!(prod_line, 1, Vector{ReactantStruct}(undef, 0))
        metadata = extract_metadata(metadata_line)
        new(sub, prod, rate, metadata)
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

function processmult(op, mult, stoich)
    if (mult isa Number) && (stoich isa Number)
        op(mult, stoich)
    else
        :($op($mult, $stoich))
    end
end

# Finds the metadata from a metadata expresion (`[key=val, ...]`) and returns as a Vector{Pair{Symbol,ExprValues}}.
function extract_metadata(metadata_line::Expr)
    metadata = :([])
    for arg in metadata_line.args
        (arg.head != :(=)) &&
            error("Malformatted metadata line: $metadata_line. Each entry in the vector should contain a `=`.")
        (arg.args[1] isa Symbol) ||
            error("Malformatted metadata entry: $arg. Entries left-hand-side should be a single symbol.")
        push!(metadata.args, :($(QuoteNode(arg.args[1])) => $(arg.args[2])))
    end
    return metadata
end

### DSL Internal Master Function ###

# Function for creating a ReactionSystem structure (used by the @reaction_network macro).
function make_reaction_system(ex::Expr; name = :(gensym(:ReactionSystem)))

    # Handle interpolation of variables
    ex = esc_dollars!(ex)

    # Read lines with reactions and options.
    reaction_lines = Expr[x for x in ex.args if x.head == :tuple]
    option_lines = Expr[x for x in ex.args if x.head == :macrocall]

    # Get macro options.
    if length(unique(arg.args[1] for arg in option_lines)) < length(option_lines)
        error("Some options where given multiple times.")
    end
    options = Dict(map(arg -> Symbol(String(arg.args[1])[2:end]) => arg,
        option_lines))

    # Reads options.
    default_reaction_metadata = :([])
    check_default_noise_scaling!(default_reaction_metadata, options)
    compound_expr, compound_species = read_compound_options(options)
    continuous_events_expr = read_events_option(options, :continuous_events)
    discrete_events_expr = read_events_option(options, :discrete_events)

    # Parses reactions, species, and parameters.
    reactions = get_reactions(reaction_lines)
    species_declared = [extract_syms(options, :species); compound_species]
    parameters_declared = extract_syms(options, :parameters)
    variables_declared = extract_syms(options, :variables)

    # Reads more options.
    vars_extracted, add_default_diff, equations = read_equations_options(
        options, variables_declared)
    variables = vcat(variables_declared, vars_extracted)

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
    all_ivs = (isnothing(sivs) ? [tiv] : [tiv; sivs.args])

    if haskey(options, :combinatoric_ratelaws)
        combinatoric_ratelaws = options[:combinatoric_ratelaws].args[end]
    else
        combinatoric_ratelaws = true
    end

    # Reads more options.
    observed_vars, observed_eqs, obs_syms = read_observed_options(
        options, [species_declared; variables], all_ivs)

    declared_syms = Set(Iterators.flatten((parameters_declared, species_declared,
        variables)))
    species_extracted, parameters_extracted = extract_species_and_parameters!(reactions,
        declared_syms)
    species = vcat(species_declared, species_extracted)
    parameters = vcat(parameters_declared, parameters_extracted)

    # Create differential expression.
    diffexpr = create_differential_expr(
        options, add_default_diff, [species; parameters; variables], tiv)

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
    vexprs = get_sexpr(vars_extracted, options, :variables; iv_symbols = ivs)
    pexprs = get_pexpr(parameters_extracted, options)
    ps, pssym = scalarize_macro(!isempty(parameters), pexprs, "ps")
    vars, varssym = scalarize_macro(!isempty(variables), vexprs, "vars")
    sps, spssym = scalarize_macro(!isempty(species), sexprs, "specs")
    comps, compssym = scalarize_macro(!isempty(compound_species), compound_expr, "comps")
    rxexprs = :(Catalyst.CatalystEqType[])
    for reaction in reactions
        push!(rxexprs.args, get_rxexprs(reaction))
    end
    for equation in equations
        push!(rxexprs.args, equation)
    end

    quote
        $ps
        $ivexpr
        $vars
        $sps
        $observed_vars
        $comps
        $diffexpr

        Catalyst.remake_ReactionSystem_internal(
            Catalyst.make_ReactionSystem_internal(
                $rxexprs, $tiv, setdiff(union($spssym, $varssym, $compssym), $obs_syms),
                $pssym; name = $name, spatial_ivs = $sivs, observed = $observed_eqs,
                continuous_events = $continuous_events_expr,
                discrete_events = $discrete_events_expr,
                combinatoric_ratelaws = $combinatoric_ratelaws);
            default_reaction_metadata = $default_reaction_metadata
        )
    end
end

### DSL Reaction Reading Functions ###

# Generates a vector containing a number of reaction structures, each containing the information about one reaction.
function get_reactions(exprs::Vector{Expr}, reactions = Vector{ReactionStruct}(undef, 0))
    for line in exprs
        # Reads core reaction information.
        arrow, rate, reaction, metadata = read_reaction_line(line)

        # Checks the type of arrow used, and creates the corresponding reaction(s). Returns them in an array.
        if in(arrow, double_arrows)
            if typeof(rate) != Expr || rate.head != :tuple
                error("Error: Must provide a tuple of reaction rates when declaring a bi-directional reaction.")
            end
            push_reactions!(reactions, reaction.args[2], reaction.args[3],
                rate.args[1], metadata.args[1], arrow)
            push_reactions!(reactions, reaction.args[3], reaction.args[2],
                rate.args[2], metadata.args[2], arrow)
        elseif in(arrow, fwd_arrows)
            push_reactions!(reactions, reaction.args[2], reaction.args[3],
                rate, metadata, arrow)
        elseif in(arrow, bwd_arrows)
            push_reactions!(reactions, reaction.args[3], reaction.args[2],
                rate, metadata, arrow)
        else
            throw("Malformed reaction, invalid arrow type used in: $(MacroTools.striplines(line))")
        end
    end
    return reactions
end

# Extracts the rate, reaction, and metadata fields (the last one optional) from a reaction line.
function read_reaction_line(line::Expr)
    # Handles rate, reaction, and arrow.
    # Special routine required for  the`-->` case, which creates different expression from all other cases.
    rate = line.args[1]
    reaction = line.args[2]
    if reaction.head == :-->
        reaction = Expr(:call, :→, reaction.args[1], reaction.args[2])
    end
    arrow = reaction.args[1]

    # Handles metadata. If not provided, empty metadata is created.
    if length(line.args) == 2
        metadata = in(arrow, double_arrows) ? :(([], [])) : :([])
    elseif length(line.args) == 3
        metadata = line.args[3]
    else
        error("The following reaction line: \"$line\" was malformed. It should have a form \"rate, reaction\" or a form \"rate, reaction, metadata\". It has neither.")
    end

    return arrow, rate, reaction, metadata
end

# Takes a reaction line and creates reaction(s) from it and pushes those to the reaction array.
# Used to create multiple reactions from, for instance, `k, (X,Y) --> 0`.
function push_reactions!(reactions::Vector{ReactionStruct}, sub_line::ExprValues,
        prod_line::ExprValues, rate::ExprValues, metadata::ExprValues, arrow::Symbol)
    # The rates, substrates, products, and metadata may be in a tupple form (e.g. `k, (X,Y) --> 0`).
    # This finds the length of these tuples (or 1 if not in tuple forms). Errors if lengs inconsistent.
    lengs = (tup_leng(sub_line), tup_leng(prod_line), tup_leng(rate), tup_leng(metadata))
    if any(!(leng == 1 || leng == maximum(lengs)) for leng in lengs)
        throw("Malformed reaction, rate=$rate, subs=$sub_line, prods=$prod_line, metadata=$metadata.")
    end

    # Loops through each reaction encoded by the reaction composites. Adds the reaction to the reactions vector.
    for i in 1:maximum(lengs)
        # If the `only_use_rate` metadata was not provided, this has to be infered from the arrow used.
        metadata_i = get_tup_arg(metadata, i)
        if all(arg -> arg.args[1] != :only_use_rate, metadata_i.args)
            push!(metadata_i.args, :(only_use_rate = $(in(arrow, pure_rate_arrows))))
        end

        # Checks that metadata fields are unqiue.
        if !allunique(arg.args[1] for arg in metadata_i.args)
            error("Some reaction metadata fields where repeated: $(metadata_entries)")
        end

        push!(reactions,
            ReactionStruct(get_tup_arg(sub_line, i),
                get_tup_arg(prod_line, i), get_tup_arg(rate, i), metadata_i))
    end
end

### DSL Species and Parameters Extraction ###

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

### DSL Output Expression Builders ###

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
        metadata = $(rxstruct.metadata)))
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

### DSL Option Handling ###

# Checks if the `@default_noise_scaling` option is used. If so, read its input and adds it as a
# default metadata value to the `default_reaction_metadata` vector.
function check_default_noise_scaling!(default_reaction_metadata, options)
    if haskey(options, :default_noise_scaling)
        if (length(options[:default_noise_scaling].args) != 3) # Becasue of how expressions are, 3 is used.
            error("@default_noise_scaling should only have a single input, this appears not to be the case: \"$(options[:default_noise_scaling])\"")
        end
        push!(default_reaction_metadata.args,
            :(:noise_scaling => $(options[:default_noise_scaling].args[3])))
    end
end

# When compound species are declared using the "@compound begin ... end" option, get a list of the compound species, and also the expression that crates them.
function read_compound_options(opts)
    # If the compound option is used retrive a list of compound species (need to be added to the reaction system's species), and the option that creates them (used to declare them as compounds at the end).
    if haskey(opts, :compounds)
        compound_expr = opts[:compounds]
        # Find compound species names, and append the independent variable.
        compound_species = [find_varinfo_in_declaration(arg.args[2])[1]
                            for arg in compound_expr.args[3].args]
    else  # If option is not used, return empty vectors and expressions.
        compound_expr = :()
        compound_species = Union{Symbol, Expr}[]
    end
    return compound_expr, compound_species
end

# Read the events (continious or discrete) provided as options to the DSL. Returns an expression which evalutes to these.
function read_events_option(options, event_type::Symbol)
    # Prepares the events, if required to, converts them to block form.
    if event_type ∉ [:continuous_events, :discrete_events]
        error("Trying to read an unsupported event type.")
    end
    events_input = haskey(options, event_type) ? options[event_type].args[3] :
                   MacroTools.striplines(:(begin end))
    events_input = option_block_form(events_input)

    # Goes throgh the events, checks for errors, and adds them to the output vector.
    events_expr = :([])
    for arg in events_input.args
        # Formatting error checks.
        # NOTE: Maybe we should move these deeper into the system (rather than the DSL), throwing errors more generally?
        if (arg isa Expr) && (arg.head != :call) || (arg.args[1] != :(=>)) ||
           (length(arg.args) != 3)
            error("Events should be on form `condition => affect`, separated by a `=>`. This appears not to be the case for: $(arg).")
        end
        if (arg isa Expr) && (arg.args[2] isa Expr) && (arg.args[2].head != :vect) &&
           (event_type == :continuous_events)
            error("The condition part of continious events (the left-hand side) must be a vector. This is not the case for: $(arg).")
        end
        if (arg isa Expr) && (arg.args[3] isa Expr) && (arg.args[3].head != :vect)
            error("The affect part of all events (the righ-hand side) must be a vector. This is not the case for: $(arg).")
        end

        # Adds the correctly formatted event to the event creation expression.
        push!(events_expr.args, arg)
    end

    return events_expr
end

# Reads the variables options. Outputs:
# `vars_extracted`: A vector with extracted variables (lhs in pure differential equations only).
# `dtexpr`: If a differentialequation is defined, the default derrivative (D ~ Differential(t)) must be defined.
# `equations`: a vector with the equations provided.
function read_equations_options(options, variables_declared)
    # Prepares the equations. First, extracts equations from provided option (converting to block form if requried).
    # Next, uses MTK's `parse_equations!` function to split input into a vector with the equations.
    eqs_input = haskey(options, :equations) ? options[:equations].args[3] : :(begin end)
    eqs_input = option_block_form(eqs_input)
    equations = Expr[]
    ModelingToolkit.parse_equations!(Expr(:block), equations,
        Dict{Symbol, Any}(), eqs_input)

    # Loops through all equations, checks for lhs of the form `D(X) ~ ...`.
    # When this is the case, the variable X and differential D are extracted (for automatic declaration).
    # Also performs simple error checks.
    vars_extracted = Vector{Symbol}()
    add_default_diff = false
    for eq in equations
        if (eq.head != :call) || (eq.args[1] != :~)
            error("Malformed equation: \"$eq\". Equation's left hand and right hand sides should be separated by a \"~\".")
        end

        # Checks if the equation have the format D(X) ~ ... (where X is a symbol). This means that the
        # default differential has been used. X is added as a declared variable to the system, and
        # we make a note that a differential D = Differential(iv) should be made as well.
        lhs = eq.args[2]
        # if lhs: is an expression. Is a function call. The function's name is D. Calls a single symbol.
        if (lhs isa Expr) && (lhs.head == :call) && (lhs.args[1] == :D) &&
           (lhs.args[2] isa Symbol)
            diff_var = lhs.args[2]
            if in(diff_var, forbidden_symbols_error)
                error("A forbidden symbol ($(diff_var)) was used as an variable in this differential equation: $eq")
            end
            add_default_diff = true
            in(diff_var, variables_declared) || push!(vars_extracted, diff_var)
        end
    end

    return vars_extracted, add_default_diff, equations
end

# Creates an expression declaring differentials. Here, `tiv` is the time independent variables,
# which is used by the default differential (if it is used).
function create_differential_expr(options, add_default_diff, used_syms, tiv)
    # Creates the differential expression.
    # If differentials was provided as options, this is used as the initial expression.
    # If the default differential (D(...)) was used in equations, this is added to the expression.
    diffexpr = (haskey(options, :differentials) ? options[:differentials].args[3] :
                MacroTools.striplines(:(begin end)))
    diffexpr = option_block_form(diffexpr)

    # Goes through all differentials, checking that they are correctly formatted and their symbol is not used elsewhere.
    for dexpr in diffexpr.args
        (dexpr.head != :(=)) &&
            error("Differential declaration must have form like D = Differential(t), instead \"$(dexpr)\" was given.")
        (dexpr.args[1] isa Symbol) ||
            error("Differential left-hand side must be a single symbol, instead \"$(dexpr.args[1])\" was given.")
        in(dexpr.args[1], used_syms) &&
            error("Differential name ($(dexpr.args[1])) is also a species, variable, or parameter. This is ambigious and not allowed.")
        in(dexpr.args[1], forbidden_symbols_error) &&
            error("A forbidden symbol ($(dexpr.args[1])) was used as a differential name.")
    end

    # If the default differential D has been used, but not pre-declared using the @differenitals
    # options, add this declaration to the list of declared differentials.
    if add_default_diff && !any(diff_dec.args[1] == :D for diff_dec in diffexpr.args)
        push!(diffexpr.args, :(D = Differential($(tiv))))
    end

    return diffexpr
end

# Reads the observables options. Outputs an expression ofr creating the obervable variables, and a vector of observable equations.
function read_observed_options(options, species_n_vars_declared, ivs_sorted)
    if haskey(options, :observables)
        # Gets list of observable equations and prepares variable declaration expression.
        # (`options[:observables]` inlucdes `@observables`, `.args[3]` removes this part)
        observed_eqs = make_observed_eqs(options[:observables].args[3])
        observed_vars = Expr(:block, :(@variables))
        obs_syms = :([])

        for (idx, obs_eq) in enumerate(observed_eqs.args)
            # Extract the observable, checks errors, and continues the loop if the observable has been declared.
            obs_name, ivs, defaults, metadata = find_varinfo_in_declaration(obs_eq.args[2])
            isempty(ivs) ||
                error("An observable ($obs_name) was given independent variable(s). These should not be given, as they are inferred automatically.")
            isnothing(defaults) ||
                error("An observable ($obs_name) was given a default value. This is forbidden.")
            (obs_name in forbidden_symbols_error) &&
                error("A forbidden symbol ($(obs_eq.args[2])) was used as an observable name.")

            # Error checks.
            if (obs_name in species_n_vars_declared) && is_escaped_expr(obs_eq.args[2])
                error("An interpoalted observable have been used, which has also been explicitly delcared within the system using eitehr @species or @variables. This is not permited.")
            end
            if ((obs_name in species_n_vars_declared) || is_escaped_expr(obs_eq.args[2])) &&
               !isnothing(metadata)
                error("Metadata was provided to observable $obs_name in the `@observables` macro. However, the obervable was also declared separately (using either @species or @variables). When this is done, metadata should instead be provided within the original @species or @variable declaration.")
            end

            # This bits adds the observables to the @variables vector which is given as output.
            # For Observables that have already been declared using @species/@variables,
            # or are interpolated, this parts should not be carried out.
            if !((obs_name in species_n_vars_declared) || is_escaped_expr(obs_eq.args[2]))
                # Appends (..) to the observable (which is later replaced with the extracted ivs).
                # Adds the observable to the first line of the output expression (starting with `@variables`).
                obs_expr = insert_independent_variable(obs_eq.args[2], :(..))
                push!(observed_vars.args[1].args, obs_expr)

                # Adds a line to the `observed_vars` expression, setting the ivs for this observable.
                # Cannot extract directly using e.g. "getfield.(dependants_structs, :reactant)" because
                # then we get something like :([:X1, :X2]), rather than :([X1, X2]).
                dep_var_expr = :(filter(!MT.isparameter,
                    Symbolics.get_variables($(obs_eq.args[3]))))
                ivs_get_expr = :(unique(reduce(
                    vcat, [arguments(MT.unwrap(dep))
                           for dep in $dep_var_expr])))
                ivs_get_expr_sorted = :(sort($(ivs_get_expr);
                    by = iv -> findfirst(MT.getname(iv) == ivs for ivs in $ivs_sorted)))
                push!(observed_vars.args,
                    :($obs_name = $(obs_name)($(ivs_get_expr_sorted)...)))
            end

            # In case metadata was given, this must be cleared from `observed_eqs`.
            # For interpolated observables (I.e. $X ~ ...) this should and cannot be done.
            is_escaped_expr(obs_eq.args[2]) || (observed_eqs.args[idx].args[2] = obs_name)

            # Adds the observable to the list of observable names.
            # This is required for filtering away so these are not added to the ReactionSystem's species list.
            # Again, avoid this check if we have interpoalted the variable.
            is_escaped_expr(obs_eq.args[2]) || push!(obs_syms.args, obs_name)
        end

        # If nothing was added to `observed_vars`, it has to be modified not to throw an error.
        (length(observed_vars.args) == 1) && (observed_vars = :())
    else
        # If option is not used, return empty expression and vector.
        observed_vars = :()
        observed_eqs = :([])
        obs_syms = :([])
    end
    return observed_vars, observed_eqs, obs_syms
end

# From the input to the @observables options, creates a vector containing one equation for each observable.
# Checks separate cases for "@obervables O ~ ..." and "@obervables begin ... end". Other cases errors.
function make_observed_eqs(observables_expr)
    if observables_expr.head == :call
        return :([$(observables_expr)])
    elseif observables_expr.head == :block
        observed_eqs = :([])
        for arg in observables_expr.args
            push!(observed_eqs.args, arg)
        end
        return observed_eqs
    else
        error("Malformed observables option usage: $(observables_expr).")
    end
end

### `@reaction` Macro & its Internals ###

@doc raw"""
@reaction

Generates a single [`Reaction`](@ref) object.

Examples:
```julia
rx = @reaction k*v, A + B --> C + D

# is equivalent to
t = default_t()
@parameters k v
@species A(t) B(t) C(t) D(t)
rx == Reaction(k*v, [A,B], [C,D])
```
Here `k` and `v` will be parameters and `A`, `B`, `C` and `D` will be variables.
Interpolation of existing parameters/variables also works
```julia
t = default_t()
@parameters k b
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

# Reads a single line and creates the corresponding ReactionStruct.
function get_reaction(line)
    reaction = get_reactions([line])
    if (length(reaction) != 1)
        error("Malformed reaction. @reaction macro only creates a single reaction. E.g. double arrows, such as `<-->` are not supported.")
    end
    return reaction[1]
end

### Generic Expression Manipulation ###

# Recursively traverses an expression and replaces special function call like "hill(...)" with the actual corresponding expression.
function recursive_expand_functions!(expr::ExprValues)
    (typeof(expr) != Expr) && (return expr)
    foreach(i -> expr.args[i] = recursive_expand_functions!(expr.args[i]),
        1:length(expr.args))
    if expr.head == :call
        !isdefined(Catalyst, expr.args[1]) && (expr.args[1] = esc(expr.args[1]))
    end
    expr
end

# Returns the length of a expression tuple, or 1 if it is not an expression tuple (probably a Symbol/Numerical).
function tup_leng(ex::ExprValues)
    (typeof(ex) == Expr && ex.head == :tuple) && (return length(ex.args))
    return 1
end

# Gets the ith element in a expression tuple, or returns the input itself if it is not an expression tuple
# (probably a  Symbol/Numerical).
function get_tup_arg(ex::ExprValues, i::Int)
    (tup_leng(ex) == 1) && (return ex)
    return ex.args[i]
end
