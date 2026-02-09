### Constants Declarations ###

# Declare various arrow types symbols used for the empty set (also 0).
const empty_set = Set{Symbol}([:∅, :Ø])
const fwd_arrows = Set{Symbol}([:>, :(=>), :→, :↣, :↦, :⇾, :⟶, :⟼, :⥟, :⥟, :⇀, :⇁, :⇒, :⟾])
const bwd_arrows = Set{Symbol}([:<, :(<=), :←, :↢, :↤, :⇽, :⟵, :⟻, :⥚, :⥞, :↼, :↽, :⇐, :⟽, Symbol("<--")])
const double_arrows = Set{Symbol}([:↔, :⟷, :⇄, :⇆, :⇌, :⇋, :⇔, :⟺, Symbol("<-->")])
const pure_rate_arrows = Set{Symbol}([:(=>), :(<=), :⇐, :⟽, :⇒, :⟾, :⇔, :⟺])

# Declares the keys used for various options.
const option_keys = (
    :species, :parameters, :variables, :discretes, :ivs, :compounds, :observables,
    :default_noise_scaling, :differentials, :equations, :continuous_events, :discrete_events,
    :brownians, :poissonians, :combinatoric_ratelaws, :require_declaration,
)

### `@species` Macro ###

# The @species macro, basically a copy of the @variables macro.
macro species(ex...)
    return Symbolics.parse_vars(:variables, Real, ex, tospecies)
end

### `@reaction_network` and `@network_component` Macros ###

"""
    @reaction_network

Macro for generating chemical reaction network models (Catalyst `ReactionSystem`s). See the
([DSL introduction](https://docs.sciml.ai/Catalyst/stable/model_creation/dsl_basics/)
and [advantage usage](https://docs.sciml.ai/Catalyst/stable/model_creation/dsl_advanced/)) sections of
the Catalyst documentation for more details on the domain-specific language (DSL) that the
macro implements. The macro's output (a `ReactionSystem` structure) is central to Catalyst
and its functionality. How to e.g. simulate these is described in the [Catalyst documentation](https://docs.sciml.ai/Catalyst/stable/).

Returns:
- A Catalyst `ReactionSystem`, i.e. a symbolic model for the reaction network. The returned
system is marked `complete`. To obtain a `ReactionSystem` that is not marked complete, for
example to then use in compositional modelling, see the otherwise equivalent `@network_component`
macro.

Examples:
Here we create a basic SIR model. It contains two reactions (infection and recovery):
```julia
sir_model = @reaction_network begin
    c1, S + I --> 2I
    c2, I --> R
end
```

Next, we create a self-activation loop. Here, a single component (`X`) activates its own production
with a Michaelis-Menten function:
```julia
sa_loop = @reaction_network begin
    mm(X,v,K), 0 --> X
    d, X --> 0
end
```
This model also contains production and degradation reactions, where `0` denotes that there are
either no substrates or no products in a reaction.

Options:
In addition to reactions, the macro also supports "option" inputs (permitting e.g. the addition
of observables). Each option is designated by a tag starting with a `@` followed by its input.
A list of options can be found [here](https://docs.sciml.ai/Catalyst/stable/api/#api_dsl_options).
"""
macro reaction_network(name::Symbol, network_expr::Expr)
    return make_rs_expr(QuoteNode(name), network_expr)
end

# The case where the name contains an interpolation.
macro reaction_network(name::Expr, network_expr::Expr)
    return make_rs_expr(esc(name.args[1]), network_expr)
end

# The case where nothing, or only a name, is provided.
macro reaction_network(name::Symbol = gensym(:ReactionSystem))
    return make_rs_expr(QuoteNode(name))
end

# Handles two disjoint cases.
macro reaction_network(expr::Expr)
    # Case 1: The input is a name with interpolation.
    (expr.head != :block) && return make_rs_expr(esc(expr.args[1]))
    # Case 2: The input is a reaction network (and no name is provided).
    return make_rs_expr(:(gensym(:ReactionSystem)), expr)
end

"""
    @network_component

Equivalent to `@reaction_network` except the generated `ReactionSystem` is not marked as
complete.
"""
macro network_component(name::Symbol, network_expr::Expr)
    return make_rs_expr(QuoteNode(name), network_expr; complete = false)
end

# The case where the name contains an interpolation.
macro network_component(name::Expr, network_expr::Expr)
    return make_rs_expr(esc(name.args[1]), network_expr; complete = false)
end

# The case where nothing, or only a name, is provided.
macro network_component(name::Symbol = gensym(:ReactionSystem))
    return make_rs_expr(QuoteNode(name); complete = false)
end

# Handles two disjoint cases.
macro network_component(expr::Expr)
    # Case 1: The input is a name with interpolation.
    (expr.head != :block) && return make_rs_expr(esc(expr.args[1]); complete = false)
    # Case 2: The input is a reaction network (and no name is provided).
    return make_rs_expr(:(gensym(:ReactionSystem)), expr; complete = false)
end

### DSL Macros Helper Functions ###

# For when no reaction network has been designated. Generates an empty network.
function make_rs_expr(name; complete = true)
    rs_expr = :(ReactionSystem(Reaction[], t, [], []; name = $name))
    complete && (rs_expr = :(complete($rs_expr)))
    return Expr(:block, :(t = default_t()), rs_expr)
end

# When both a name and a network expression are generated, dispatch these to the internal
# `make_reaction_system` function.
function make_rs_expr(name, network_expr; complete = true)
    rs_expr = make_reaction_system(striplines(network_expr), name)
    complete && (rs_expr = :(complete($rs_expr)))
    return rs_expr
end

### Internal DSL Structures ###

# Internal structure containing information about one reactant in one reaction.
struct DSLReactant
    reactant::Union{Symbol, Expr}
    stoichiometry::ExprValues
end

# Internal structure containing information about one Reaction. Contains all its substrates and
# products as well as its rate and potential metadata. Uses a specialized constructor.
struct DSLReaction
    substrates::Vector{DSLReactant}
    products::Vector{DSLReactant}
    rate::ExprValues
    metadata::Expr
    rxexpr::Expr

    function DSLReaction(
            sub_line::ExprValues, prod_line::ExprValues,
            rate::ExprValues, metadata_line::ExprValues, rx_line::Expr
        )
        subs = recursive_find_reactants!(sub_line, 1, Vector{DSLReactant}(undef, 0))
        prods = recursive_find_reactants!(prod_line, 1, Vector{DSLReactant}(undef, 0))
        metadata = extract_metadata(metadata_line)
        return new(subs, prods, rate, metadata, rx_line)
    end
end

# Recursive function that loops through the reaction line and finds the reactants and their
# stoichiometry. Recursion makes it able to handle weird cases like 2(X + Y + 3(Z + XY)). The
# reactants are stored in the `reactants` vector. As the expression tree is parsed, the
# stoichiometry is updated and new reactants are added.
function recursive_find_reactants!(
        ex::ExprValues, mult::ExprValues,
        reactants::Vector{DSLReactant}
    )
    # We have reached the end of the expression tree and can finalise and return the reactants.
    if (typeof(ex) != Expr) || (ex.head == :escape) || (ex.head == :ref)
        # The final bit of the expression is not a relevant reactant, no additions are required.
        (ex == 0 || in(ex, empty_set)) && (return reactants)

        # If the expression corresponds to a reactant on our list, increase its multiplicity.
        idx = findfirst(r.reactant == ex for r in reactants)
        if !isnothing(idx)
            newmult = processmult(+, mult, reactants[idx].stoichiometry)
            reactants[idx] = DSLReactant(ex, newmult)

            # If the expression corresponds to a new reactant, add it to the list.
        else
            push!(reactants, DSLReactant(ex, mult))
        end

        # If we have encountered a multiplication (i.e. a stoichiometry and a set of reactants).
    elseif ex.args[1] == :*
        # The normal case (e.g. 3*X or 3*(X+Y)). Update the current multiplicity and continue.
        if length(ex.args) == 3
            newmult = processmult(*, mult, ex.args[2])
            recursive_find_reactants!(ex.args[3], newmult, reactants)
            # More complicated cases (e.g. 2*3*X). Yes, `ex.args[1:(end - 1)]` should start at 1 (not 2).
        else
            newmult = processmult(*, mult, Expr(:call, ex.args[1:(end - 1)]...))
            recursive_find_reactants!(ex.args[end], newmult, reactants)
        end
        # If we have encountered a sum of different reactants, apply recursion on each.
    elseif ex.args[1] == :+
        for i in 2:length(ex.args)
            recursive_find_reactants!(ex.args[i], mult, reactants)
        end
    else
        throw("Malformed reaction, bad operator: $(ex.args[1]) found in stoichiometry expression $ex.")
    end
    return reactants
end

# Helper function for updating the multiplicity throughout recursion (handles e.g. parametric
# stoichiometries). The `op` argument is an operation (e.g. `*`, but could also e.g. be `+`).
function processmult(op, mult, stoich)
    return if (mult isa Number) && (stoich isa Number)
        op(mult, stoich)
    else
        :($op($mult, $stoich))
    end
end

# Finds the metadata from a metadata expression (`[key=val, ...]`) and returns as a
# `Vector{Pair{Symbol,ExprValues}}``.
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

### Specialised @require_declaration Option Error ###
struct UndeclaredSymbolicError <: Exception
    msg::String
end

function Base.showerror(io::IO, err::UndeclaredSymbolicError)
    print(io, "UndeclaredSymbolicError: ")
    return print(io, err.msg)
end

### DSL Internal Master Function ###

# Function for creating a ReactionSystem structure (used by the @reaction_network macro).
function make_reaction_system(ex::Expr, name)
    # Handle interpolation of variables in the input.
    ex = esc_dollars!(ex)

    # Extract the lines with reactions, the lines with options, and the options. Check for input errors.
    reaction_lines = Expr[x for x in ex.args if x.head == :tuple]
    option_lines = Expr[x for x in ex.args if x.head == :macrocall]
    allunique(arg.args[1] for arg in option_lines) ||
        error("Some options where given multiple times.")
    numlines = length(reaction_lines) + length(option_lines)
    (numlines != length(ex.args)) &&
        error("@reaction_network input contain $(length(ex.args) - numlines) malformed lines.")
    options = Dict(Symbol(String(arg.args[1])[2:end]) => arg for arg in option_lines)
    any(!in(option_keys), keys(options)) &&
        error("The following unsupported options were used: $(filter(opt_in -> !in(opt_in, option_keys), keys(options)))")

    # Read options that explicitly declare some symbol (e.g. `@species`). Compiles a list of
    # all declared symbols and checks that there have been no double-declarations.
    sps_declared = extract_syms(options, :species)
    ps_declared = extract_syms(options, :parameters)
    vs_declared = extract_syms(options, :variables)
    discs_declared = extract_syms(options, :discretes)
    tiv, sivs, ivs, ivsexpr = read_ivs_option(options)
    cmpexpr_init, cmps_declared = read_compounds_option(options)
    diffsexpr, diffs_declared = read_differentials_option(options)
    brownsexpr_init, browns_declared = read_brownians_option(options)
    poissexpr_init, poiss_declared = read_poissonians_option(options)
    syms_declared = collect(
        Iterators.flatten(
            (
                cmps_declared, sps_declared, ps_declared,
                vs_declared, discs_declared, ivs, diffs_declared, browns_declared, poiss_declared,
            )
        )
    )
    if !allunique(syms_declared)
        nonunique_syms = [s for s in syms_declared if count(x -> x == s, syms_declared) > 1]
        error("The following symbols $(unique(nonunique_syms)) have explicitly been declared as multiple types of components (e.g. occur in at least two of the `@species`, `@parameters`, `@variables`, `@ivs`, `@compounds`, `@differentials`). This is not allowed.")
    end

    # Reads the reactions and equation. From these, infer species, variables, and parameters.
    requiredec = haskey(options, :require_declaration)
    reactions = get_reactions(reaction_lines)
    sps_inferred, ps_pre_inferred, stoich_ps = extract_sps_and_ps(reactions, syms_declared; requiredec)
    vs_inferred, diffs_inferred, equations = read_equations_option!(
        diffsexpr, options,
        union(syms_declared, sps_inferred), tiv; requiredec
    )
    ps_inferred = setdiff(ps_pre_inferred, vs_inferred, diffs_inferred)
    syms_inferred = union(sps_inferred, ps_inferred, vs_inferred, diffs_inferred)
    all_syms = union(syms_declared, syms_inferred)
    validate_poissonian_rate_syms(options, all_syms)
    obsexpr, obs_eqs, obs_syms = read_observables_option(
        options, ivs,
        union(sps_declared, vs_declared), all_syms; requiredec
    )

    # Read options not related to the declaration or inference of symbols.
    discs_inferred = Vector{Symbol}()
    continuous_events_expr = read_events_option!(options, discs_inferred, ps_inferred, discs_declared, :continuous_events)
    discrete_events_expr = read_events_option!(options, discs_inferred, ps_inferred, discs_declared, :discrete_events)
    default_reaction_metadata = read_default_noise_scaling_option(options)
    combinatoric_ratelaws = read_combinatoric_ratelaws_option(options)

    # Creates expressions corresponding to actual code from the internal DSL representation.
    psexpr_init = get_psexpr(ps_inferred, stoich_ps, options)
    spsexpr_init = get_usexpr(sps_inferred, options; ivs)
    vsexpr_init = get_usexpr(vs_inferred, options, :variables; ivs)
    discsexpr_init = get_usexpr(discs_inferred, options, :discretes; ivs)
    psexpr, psvar = assign_var_to_symvar_declaration(psexpr_init, "ps", scalarize = false)
    spsexpr, spsvar = assign_var_to_symvar_declaration(spsexpr_init, "specs")
    vsexpr, vsvar = assign_var_to_symvar_declaration(vsexpr_init, "vars")
    discsexpr, discsvar = assign_var_to_symvar_declaration(discsexpr_init, "discs")
    cmpsexpr, cmpsvar = assign_var_to_symvar_declaration(cmpexpr_init, "comps")
    brownsexpr, brownsvar = assign_var_to_symvar_declaration(brownsexpr_init, "brownians", scalarize = false)
    poissexpr, poissvar = assign_var_to_symvar_declaration(poissexpr_init, "poissonians", scalarize = false)
    rxsexprs = get_rxexprs(reactions, equations, all_syms)

    # Assemblies the full expression that declares all required symbolic variables, and
    # then the output `ReactionSystem`.
    return MacroTools.flatten(
        striplines(
            quote
                # Inserts the expressions which generate the `ReactionSystem` input.
                $ivsexpr
                $psexpr
                $vsexpr
                $spsexpr
                $discsexpr
                $obsexpr
                $cmpsexpr
                $diffsexpr
                $brownsexpr
                $poissexpr

                # Stores each kwarg in a variable. Not necessary, but useful when debugging generated code.
                name = $name
                spatial_ivs = $sivs
                rx_eq_vec = $rxsexprs
                us = setdiff(union($spsvar, $vsvar, $cmpsvar), $obs_syms)
                ps = union($psvar, $discsvar)
                _observed = $obs_eqs
                _continuous_events = $continuous_events_expr
                _discrete_events = $discrete_events_expr
                _combinatoric_ratelaws = $combinatoric_ratelaws
                _default_reaction_metadata = $default_reaction_metadata

                remake_ReactionSystem_internal(
                    make_ReactionSystem_internal(
                        rx_eq_vec, $tiv, us, ps, $brownsvar; poissonians = $poissvar,
                        name, spatial_ivs, observed = _observed, continuous_events = _continuous_events,
                        discrete_events = _discrete_events, combinatoric_ratelaws = _combinatoric_ratelaws
                    );
                    default_reaction_metadata = _default_reaction_metadata
                )
            end
        )
    )
end

### DSL Reaction Reading Functions ###

# Generates a vector of reaction structures, each containing information about one reaction.
function get_reactions(exprs::Vector{Expr})
    # Declares an array to which we add all found reactions.
    reactions = Vector{DSLReaction}(undef, 0)

    # Loops through each line of reactions. Extracts and adds each line's reactions to `reactions`.
    for line in exprs
        # Reads core reaction information.
        arrow, rate, reaction, metadata = read_reaction_line(line)

        # Currently, reaction bundling where rates (but neither substrates nor products) are
        # bundled, is disabled. See discussion in https://github.com/SciML/Catalyst.jl/issues/1219.
        if !in(arrow, double_arrows) && Meta.isexpr(rate, :tuple) &&
                !Meta.isexpr(reaction.args[2], :tuple) && !Meta.isexpr(reaction.args[3], :tuple)
            error("Bundling of reactions with multiple rates but singular substrates and product sets is disallowed. This error is potentially due to a bidirectional (`<-->`) reaction being incorrectly typed as `-->`.")
        end

        # Checks which type of line is used, and calls `push_reactions!` on the processed line.
        if in(arrow, double_arrows)
            (typeof(rate) != Expr || rate.head != :tuple) &&
                error("Error: Must provide a tuple of reaction rates when declaring a bi-directional reaction.")
            (typeof(metadata) != Expr || metadata.head != :tuple) &&
                error("Error: Must provide a tuple of reaction metadata when declaring a bi-directional reaction.")

            push_reactions!(
                reactions, reaction.args[2], reaction.args[3],
                rate.args[1], metadata.args[1], arrow, line
            )
            push_reactions!(
                reactions, reaction.args[3], reaction.args[2],
                rate.args[2], metadata.args[2], arrow, line
            )
        elseif in(arrow, fwd_arrows)
            push_reactions!(
                reactions, reaction.args[2], reaction.args[3],
                rate, metadata, arrow, line
            )
        elseif in(arrow, bwd_arrows)
            push_reactions!(
                reactions, reaction.args[3], reaction.args[2],
                rate, metadata, arrow, line
            )
        else
            throw("Malformed reaction, invalid arrow type used in: $(striplines(line))")
        end
    end
    return reactions
end

# Extract the rate, reaction, and metadata fields (the last one optional) from a reaction line.
function read_reaction_line(line::Expr)
    # Handles rate, reaction, and arrow. A special routine is required for  the`-->` case,
    # which creates an expression different from what the other arrows create.
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

# Takes a reaction line and creates reaction(s) from it and pushes those to the reaction vector.
# Used to create multiple reactions from bundled reactions (like `k, (X,Y) --> 0`).
function push_reactions!(
        reactions::Vector{DSLReaction}, subs::ExprValues,
        prods::ExprValues, rate::ExprValues, metadata::ExprValues, arrow::Symbol, line::Expr
    )
    # The rates, substrates, products, and metadata may be in a tuple form (e.g. `k, (X,Y) --> 0`).
    # This finds these tuples' lengths (or 1 for non-tuple forms). Inconsistent lengths yield error.
    lengs = (tup_leng(subs), tup_leng(prods), tup_leng(rate), tup_leng(metadata))
    maxlen = maximum(lengs)
    any(!(leng == 1 || leng == maxlen) for leng in lengs) &&
        error("Malformed reaction, rate: $rate, subs: $subs, prods: $prods, metadata: $metadata.")

    # Loops through each reaction encoded by the reaction's different components.
    # Creates a `DSLReaction` representation and adds it to `reactions`.
    for i in 1:maxlen
        # If the `only_use_rate` metadata was not provided, this must be inferred from the arrow.
        metadata_i = get_tup_arg(metadata, i)
        if all(arg.args[1] != :only_use_rate for arg in metadata_i.args)
            push!(metadata_i.args, :(only_use_rate = $(in(arrow, pure_rate_arrows))))
        end

        # Checks that metadata fields are unique.
        if !allunique(arg.args[1] for arg in metadata_i.args)
            error("Some reaction metadata fields where repeated: $(metadata_entries)")
        end

        # Extracts substrates, products, and rates for the i'th reaction.
        subs_i, prods_i, rate_i = get_tup_arg.((subs, prods, rate), i)
        push!(reactions, DSLReaction(subs_i, prods_i, rate_i, metadata_i, line))
    end
    return
end

### DSL Species and Parameters Extraction ###

# When the user have used the @species (or @parameters) options, extract species (or
# parameters) from its input.
function extract_syms(opts, vartype::Symbol)
    # If the corresponding option has been used, use `Symbolics._parse_vars` to find all
    # variable within it (returning them in a vector).
    return if haskey(opts, vartype)
        ex = opts[vartype]
        vars = Symbolics._parse_vars(vartype, Real, ex.args[3:end])
        Vector{Union{Symbol, Expr}}(vars.args[end].args)
    else
        Union{Symbol, Expr}[]
    end
end

# Function looping through all reactions, to find undeclared symbols (species or
# parameters) and assign them to the right category.
# `stoich_ps` records parameters used in stoichiometries (if these are not declare separately,
# these are infered to be integers)..
function extract_sps_and_ps(reactions, excluded_syms; requiredec = false)
    # Loops through all reactants and extract undeclared ones as species.
    species = OrderedSet{Union{Symbol, Expr}}()
    for reaction in reactions
        for reactant in Iterators.flatten((reaction.substrates, reaction.products))
            add_syms_from_expr!(species, reactant.reactant, excluded_syms)
        end
        (!isempty(species) && requiredec) &&
            throw(UndeclaredSymbolicError("Unrecognized reactant $(join(species, ", ")) detected in reaction expression: \"$(string(reaction.rxexpr))\". Since the flag @require_declaration is declared, all species must be explicitly declared with the @species option."))
    end
    excluded_syms = union(excluded_syms, species)

    # Loops through all rates and stoichiometries, extracting used symbols as parameters.
    parameters = OrderedSet{Union{Symbol, Expr}}()
    stoich_ps = OrderedSet{Union{Symbol, Expr}}()
    for reaction in reactions
        add_syms_from_expr!(parameters, reaction.rate, excluded_syms)
        (!isempty(parameters) && requiredec) &&
            throw(UndeclaredSymbolicError("Unrecognized symbol $(join(parameters, ", ")) detected in rate expression: $(reaction.rate) for the following reaction expression: \"$(string(reaction.rxexpr))\". Since the flag @require_declaration is declared, all parameters must be explicitly declared with the @parameters option."))
        for reactant in Iterators.flatten((reaction.substrates, reaction.products))
            add_syms_from_expr!(parameters, reactant.stoichiometry, excluded_syms, stoich_ps)
            (!isempty(parameters) && requiredec) &&
                throw(UndeclaredSymbolicError("Unrecognized symbol $(join(parameters, ", ")) detected in the stoichiometry for reactant $(reactant.reactant) in the following reaction expression: \"$(string(reaction.rxexpr))\". Since the flag @require_declaration is declared, all parameters must be explicitly declared with the @parameters option."))
        end
    end

    return collect(species), collect(parameters), collect(stoich_ps)
end

# Function called by `extract_sps_and_ps`, recursively loops through an expression and find
# symbols (adding them to the push_symbols vector). Returns `nothing` to ensure type stability.
function add_syms_from_expr!(push_symbols::AbstractSet, expr::ExprValues, excluded_syms, push_symbols2 = nothing)
    # If we have encountered a Symbol in the recursion, we can try extracting it.
    if expr isa Symbol
        if !(expr in forbidden_symbols_skip) && !(expr in excluded_syms)
            push!(push_symbols, expr)
            isnothing(push_symbols2) || push!(push_symbols2, expr)
        end
    elseif expr isa Expr
        # note, this (correctly) skips $(...) expressions
        for i in 2:length(expr.args)
            add_syms_from_expr!(push_symbols, expr.args[i], excluded_syms, push_symbols2)
        end
    end
    return nothing
end

### DSL Output Expression Builders ###

# Given the parameters that were extracted from the reactions, and the options dictionary,
# creates the `@parameters ...` expression for the macro output.
function get_psexpr(parameters_extracted, stoich_ps, options)
    pexprs = if haskey(options, :parameters)
        options[:parameters]
    elseif isempty(parameters_extracted)
        :()
    else
        :(@parameters)
    end
    arg_vec = ((length(pexprs.args) > 2) && Meta.isexpr(pexprs.args[3], :block)) ?
        (pexprs.args[3].args) : (pexprs.args)
    foreach(p -> push!(arg_vec, p), setdiff(parameters_extracted, stoich_ps))
    foreach(p -> push!(arg_vec, :($p::Int64)), stoich_ps)
    return pexprs
end

# Given the extracted species (or variables) and the option dictionary, create the
# `@species ...` (or `@variables ..`) expression which would declare these.
# If `key = :variables`, does this for variables (and not species).
function get_usexpr(us_extracted, options, key = :species; ivs = (DEFAULT_IV_SYM,))
    usexpr = if haskey(options, key)
        options[key]
    elseif isempty(us_extracted)
        :()
    else
        Expr(:macrocall, Symbol("@", key), LineNumberNode(0))
    end
    arg_vec = ((length(usexpr.args) > 2) && Meta.isexpr(usexpr.args[3], :block)) ?
        (usexpr.args[3].args) : (usexpr.args)
    for u in us_extracted
        u isa Symbol && push!(arg_vec, Expr(:call, u, ivs...))
    end
    return usexpr
end

# From the system reactions (as `DSLReaction`s) and equations (as expressions),
# creates the expression that evaluates to the reaction (+ equations) vector.
function get_rxexprs(reactions, equations, all_syms)
    rxsexprs = :(Catalyst.CatalystEqType[])
    foreach(rx -> push!(rxsexprs.args, get_rxexpr(rx)), reactions)
    foreach(eq -> push!(rxsexprs.args, escape_equation!(eq, all_syms)), equations)
    return rxsexprs
end

# From a `DSLReaction` struct, creates the expression which evaluates to the creation
# of the corresponding reaction.
function get_rxexpr(rx::DSLReaction)
    # Initiates the `Reaction` expression.
    rate = recursive_escape_functions!(rx.rate)
    subs_init = isempty(rx.substrates) ? nothing : :([])
    subs_stoich_init = deepcopy(subs_init)
    prod_init = isempty(rx.products) ? nothing : :([])
    prod_stoich_init = deepcopy(prod_init)
    rx_constructor = :(
        Reaction(
            $rate, $subs_init, $prod_init, $subs_stoich_init,
            $prod_stoich_init; metadata = $(rx.metadata)
        )
    )

    # Loops through all products and substrates, and adds them (and their stoichiometries)
    # to the `Reaction` expression.
    for sub in rx.substrates
        push!(rx_constructor.args[4].args, sub.reactant)
        push!(rx_constructor.args[6].args, sub.stoichiometry)
    end
    for prod in rx.products
        push!(rx_constructor.args[5].args, prod.reactant)
        push!(rx_constructor.args[7].args, prod.stoichiometry)
    end
    return rx_constructor
end

# Takes a ModelingToolkit declaration macro (like @parameters ...) and return and expression:
# That calls the macro and then scalarizes all the symbols created into a vector of Nums.
# stores the created symbolic variables in a variable (whose name is generated from `name`).
# It will also return the name used for the variable that stores the symbolic variables.
# If requested, performs scalarization.
function assign_var_to_symvar_declaration(expr_init, name; scalarize = true)
    # Generates a random variable name which (in generated code) will store the produced
    # symbolic variables (e.g. `var"##ps#384"`).
    namesym = gensym(name)

    # If the input expression is non-empty, wrap it with additional information.
    if expr_init != :(())
        if scalarize
            symvec = gensym()
            expr = quote
                $symvec = $expr_init
                $namesym = reduce(vcat, Symbolics.scalarize($symvec))
            end
        else
            expr = quote
                $namesym = $expr_init
            end
        end
    else
        expr = :($namesym = Num[])
    end
    return expr, namesym
end

# Recursively escape functions within equations of an equation written using user-defined functions.
# Does not escape special function calls like "hill(...)" and differential operators. Does
# also not escape stuff corresponding to e.g. species or parameters (required for good error
# for when e.g. a species is used as a differential, or for time delays in the future).
function escape_equation!(eqexpr::Expr, all_syms)
    eqexpr.args[2] = recursive_escape_functions!(eqexpr.args[2], all_syms)
    eqexpr.args[3] = recursive_escape_functions!(eqexpr.args[3], all_syms)
    return eqexpr
end

### DSL Option Handling ###

# Finds the time independent variable, and any potential spatial independent variables.
# Returns these (individually and combined), as well as an expression for declaring them.
function read_ivs_option(options)
    # Creates the independent variables expressions (depends on whether the `ivs` option was used).
    if haskey(options, :ivs)
        ivs = Tuple(extract_syms(options, :ivs))
        ivsexpr = copy(options[:ivs])
        ivsexpr.args[1] = Symbol("@", "independent_variables")
    else
        ivs = (DEFAULT_IV_SYM,)
        ivsexpr = :($(DEFAULT_IV_SYM) = default_t())
    end

    # Extracts the independent variables symbols (time and spatial), and returns the output.
    tiv = ivs[1]
    sivs = (length(ivs) > 1) ? Expr(:vect, ivs[2:end]...) : nothing
    return tiv, sivs, ivs, ivsexpr
end

# When compound species are declared using the "@compound begin ... end" option, get a list
# of the compound species, and also the expression that creates them.
function read_compounds_option(options)
    # If the compound option is used, retrieve a list of compound species and  the option line
    # that creates them (used to declare them as compounds at the end). Due to some expression
    # handling, in the case of a single compound we must change to the `@compound` macro.
    if haskey(options, :compounds)
        cmpexpr_init = options[:compounds]
        cmpexpr_init.args[3] = option_block_form(get_block_option(cmpexpr_init))
        cmps_declared = [
            find_varinfo_in_declaration(arg.args[2])[1]
                for arg in cmpexpr_init.args[3].args
        ]
        (length(cmps_declared) == 1) && (cmpexpr_init.args[1] = Symbol("@compound"))
    else  # If option is not used, return empty vectors and expressions.
        cmpexpr_init = :()
        cmps_declared = Union{Symbol, Expr}[]
    end
    return cmpexpr_init, cmps_declared
end

# Creates an expression declaring differentials. Here, `tiv` is the time independent variables,
# which is used by the default differential (if it is used).
function read_differentials_option(options)
    # Creates the differential expression.
    # If differentials were provided as options, this is used as the initial expression.
    # If the default differential (D(...)) was used in equations, this is added to the expression.
    diffsexpr = (
        haskey(options, :differentials) ?
            get_block_option(options[:differentials]) : striplines(:(begin end))
    )
    diffsexpr = option_block_form(diffsexpr)

    # Goes through all differentials, checking that they are correctly formatted. Adds their
    # symbol to the list of declared differential symbols.
    diffs_declared = Union{Symbol, Expr}[]
    for dexpr in diffsexpr.args
        (dexpr.head != :(=)) &&
            error("Differential declaration must have form like D = Differential(t), instead \"$(dexpr)\" was given.")
        (dexpr.args[1] isa Symbol) ||
            error("Differential left-hand side must be a single symbol, instead \"$(dexpr.args[1])\" was given.")
        in(dexpr.args[1], forbidden_symbols_error) &&
            error("A forbidden symbol ($(dexpr.args[1])) was used as a differential name.")
        push!(diffs_declared, dexpr.args[1])
    end

    return diffsexpr, diffs_declared
end

# Creates the initial expression for declaring brownians. Also extracts any symbols
# declared as brownians by the `@brownian` option.
function read_brownians_option(options)
    browns_declared = extract_syms(options, :brownians)
    brownsexpr_init = haskey(options, :brownians) ? options[:brownians] : :()
    return brownsexpr_init, browns_declared
end

# Creates the initial expression for declaring poissonians. Also extracts any symbols
# declared as poissonians by the `@poissonians` option. Unlike brownians (plain symbols),
# poissonians use call-expression syntax (e.g. `dN(λ)`), so custom name extraction is needed.
function read_poissonians_option(options)
    poiss_declared = extract_poissonian_names(options)
    poissexpr_init = haskey(options, :poissonians) ? options[:poissonians] : :()
    return poissexpr_init, poiss_declared
end

# Extract symbol names from @poissonians call expressions (e.g., dN(λ) → :dN).
function extract_poissonian_names(options)
    !haskey(options, :poissonians) && return Union{Symbol, Expr}[]
    ex = options[:poissonians]
    names = Union{Symbol, Expr}[]
    for arg in ex.args[3:end]  # Skip macrocall head and LineNumberNode
        arg isa LineNumberNode && continue
        if arg isa Expr && arg.head == :block
            for inner in arg.args
                inner isa LineNumberNode && continue
                (inner isa Expr && inner.head == :call) && push!(names, inner.args[1])
            end
        elseif arg isa Expr && arg.head == :call
            push!(names, arg.args[1])
        end
    end
    return names
end

# Extract bare symbols from poissonian rate expressions, skipping escaped (interpolated) subtrees.
# Used to validate that all symbols in rates are pre-declared.
function extract_poissonian_rate_syms(options)
    !haskey(options, :poissonians) && return Symbol[]
    rate_exprs = Any[]
    ex = options[:poissonians]
    for arg in ex.args[3:end]
        arg isa LineNumberNode && continue
        if arg isa Expr && arg.head == :block
            for inner in arg.args
                inner isa LineNumberNode && continue
                if inner isa Expr && inner.head == :call && length(inner.args) >= 2
                    push!(rate_exprs, inner.args[2])
                end
            end
        elseif arg isa Expr && arg.head == :call && length(arg.args) >= 2
            push!(rate_exprs, arg.args[2])
        end
    end
    syms = Symbol[]
    for rexpr in rate_exprs
        _collect_symbols!(syms, rexpr)
    end
    return syms
end

# Recursively collect bare Symbol names from an expression, skipping escaped nodes.
function _collect_symbols!(syms::Vector{Symbol}, ex)
    return if ex isa Symbol
        push!(syms, ex)
    elseif ex isa Expr
        is_escaped_expr(ex) && return
        # For function/operator calls, skip the call head (e.g. `*`, `sin`) and
        # only inspect argument expressions for user-declared symbols.
        args = (ex.head == :call) ? ex.args[2:end] : ex.args
        for arg in args
            _collect_symbols!(syms, arg)
        end
    end
end

# Validates that all non-interpolated symbols in @poissonians rate expressions are declared
# or inferred. Throws an `UndeclaredSymbolicError` for any unrecognized symbols.
function validate_poissonian_rate_syms(options, all_syms)
    poiss_rate_syms = extract_poissonian_rate_syms(options)
    undeclared = setdiff(poiss_rate_syms, all_syms)
    return if !isempty(undeclared)
        throw(
            UndeclaredSymbolicError(
                "Unrecognized symbol(s) $(join(undeclared, ", ")) in a `@poissonians` rate " *
                    "expression. Symbols in poissonian rates must be pre-declared (e.g. via " *
                    "`@parameters`, `@species`, or `@variables`) or interpolated."
            )
        )
    end
end

# Reads the variables options. Outputs a list of the variables inferred from the equations,
# as well as the equation vector. If the default differential was used, update the `diffsexpr`
# expression so that this declares this as well.
function read_equations_option!(
        diffsexpr, options, syms_unavailable, tiv; requiredec = false
    )
    # Prepares the equations. First, extract equations from the provided option (converting to block form if required).
    # Next, uses `parse_equations!` function to split input into a vector with the equations.
    eqs_input = haskey(options, :equations) ? get_block_option(options[:equations]) :
        MacroTools.striplines(:(begin end))
    eqs_input = option_block_form(eqs_input)
    equations = eqs_input.args

    # Loops through all equations, checks for lhs of the form `D(X) ~ ...`.
    # When this is the case, the variable X and differential D are extracted (for automatic declaration).
    # Also performs simple error checks.
    vs_inferred = OrderedSet{Union{Symbol, Expr}}()
    add_default_diff = false
    for eq in equations
        if (eq.head != :call) || (eq.args[1] != :~)
            error("Malformed equation: \"$eq\". Equation's left hand and right-hand sides should be separated by a \"~\".")
        end

        # If the default differential (`D`) is used, record that it should be declared later on.
        if (:D ∉ syms_unavailable) && find_D_call(eq)
            requiredec && throw(
                UndeclaredSymbolicError(
                    "Unrecognized symbol D was used as a differential in an equation: \"$eq\". Since the @require_declaration flag is set, all differentials in equations must be explicitly declared using the @differentials option."
                )
            )
            add_default_diff = true
            push!(syms_unavailable, :D)
        end

        # Any undeclared symbolic variables encountered should be extracted as variables.
        add_syms_from_expr!(vs_inferred, eq, syms_unavailable)
        (!isempty(vs_inferred) && requiredec) && throw(
            UndeclaredSymbolicError(
                "Unrecognized symbol $(join(vs_inferred, ", ")) detected in equation expression: \"$(string(eq))\". Since the flag @require_declaration is declared, all symbolic variables must be explicitly declared with the @species, @variables, and @parameters options."
            )
        )
    end

    # If `D` differential is used, add it to the differential expression and inferred differentials list.
    diffs_inferred = Union{Symbol, Expr}[]
    if add_default_diff && !any(diff_dec.args[1] == :D for diff_dec in diffsexpr.args)
        diffs_inferred = [:D]
        push!(diffsexpr.args, :(D = Differential($(tiv))))
    end

    return collect(vs_inferred), diffs_inferred, equations
end

# Searches an expression `expr` and returns true if it has any subexpression `D(...)` (where `...` can be anything).
# Used to determine whether the default differential D has been used in any equation provided to `@equations`.
function find_D_call(expr)
    return if Meta.isexpr(expr, :call) && expr.args[1] == :D
        true
    elseif expr isa Expr
        any(find_D_call, expr.args)
    else
        false
    end
end

# Reads the observables options. Outputs an expression for creating the observable variables,
# a vector containing the observable equations, and a list of all observable symbols (this
# list contains both those declared separately or inferred from the `@observables` option` input`).
function read_observables_option(
        options, all_ivs, us_declared, all_syms; requiredec = false
    )
    syms_unavailable = setdiff(all_syms, us_declared)
    if haskey(options, :observables)
        # Gets list of observable equations and prepares variable declaration expression.
        # (`options[:observables]` includes `@observables`, `.args[3]` removes this part)
        obs_eqs = make_obs_eqs(get_block_option(options[:observables]))
        obsexpr = Expr(:block, :(@variables))
        obs_syms = :([])

        for (idx, obs_eq) in enumerate(obs_eqs.args)
            # Extract the observable, checks for errors.
            obs_name, ivs, _, defaults,
                metadata = find_varinfo_in_declaration(obs_eq.args[2])

            # Error checks.
            (requiredec && !in(obs_name, us_declared)) &&
                throw(UndeclaredSymbolicError("An undeclared symbol ($obs_name) was used as an observable in the following observable equation: \"$obs_eq\". Since the flag @require_declaration is set, all observables must be declared with either the @species or @variables options."))
            isempty(ivs) ||
                error("An observable ($obs_name) was given independent variable(s). These should not be given, as they are inferred automatically.")
            isnothing(defaults) ||
                error("An observable ($obs_name) was given a default value. This is forbidden.")
            in(obs_name, forbidden_symbols_error) &&
                error("A forbidden symbol ($(obs_eq.args[2])) was used as an observable name.")
            in(obs_name, syms_unavailable) &&
                error("An observable ($obs_name) uses a name that already have been already been declared or inferred as another model property.")
            (obs_name in us_declared) && is_escaped_expr(obs_eq.args[2]) &&
                error("An interpolated observable have been used, which has also been explicitly declared within the system using either @species or @variables. This is not permitted.")
            ((obs_name in us_declared) || is_escaped_expr(obs_eq.args[2])) &&
                !isnothing(metadata) &&
                error("Metadata was provided to observable $obs_name in the `@observables` macro. However, the observable was also declared separately (using either @species or @variables). When this is done, metadata should instead be provided within the original @species or @variable declaration.")

            # This bit adds the observables to the @variables vector which is given as output.
            # For Observables that have already been declared using @species/@variables,
            # or are interpolated, this parts should not be carried out.
            if !((obs_name in us_declared) || is_escaped_expr(obs_eq.args[2]))
                # Creates an expression which extracts the ivs of the species & variables the
                # observable depends on, and splats them out in the correct order.
                dep_var_expr = :(
                    filter(
                        !MT.isparameter,
                        Symbolics.get_variables($(obs_eq.args[3]))
                    )
                )
                ivs_get_expr = :(
                    unique(
                        reduce(
                            vcat, [
                                sorted_arguments(unwrap(dep))
                                    for dep in $dep_var_expr
                            ]
                        )
                    )
                )
                ivs_get_expr_sorted = :(
                    sort(
                        $(ivs_get_expr);
                        by = iv -> findfirst(MT.getname(iv) == ivs for ivs in $all_ivs)
                    )
                )

                obs_expr = insert_independent_variable(obs_eq.args[2], :($ivs_get_expr_sorted...))
                push!(obsexpr.args[1].args, obs_expr)
            end

            # In case metadata was given, this must be cleared from `obs_eqs`.
            # For interpolated observables (I.e. $X ~ ...) this should and cannot be done.
            is_escaped_expr(obs_eq.args[2]) || (obs_eqs.args[idx].args[2] = obs_name)

            # Adds the observable to the list of observable names.
            # This is required for filtering away so these are not added to the ReactionSystem's species list.
            # Again, avoid this check if we have interpolated the variable.
            is_escaped_expr(obs_eq.args[2]) || push!(obs_syms.args, obs_name)
        end

        # If nothing was added to `obsexpr`, it has to be modified not to throw an error.
        (striplines(obsexpr) == striplines(Expr(:block, :(@variables)))) &&
            (obsexpr = :())
    else
        # If observables option is not used, return empty expression and vector.
        obsexpr = :()
        obs_eqs = :([])
        obs_syms = :([])
    end

    return obsexpr, obs_eqs, obs_syms
end

# From the input to the @observables options, create a vector containing one equation for
# each observable. `option_block_form` handles if single line declaration of `@observables`,
# i.e. without a `begin ... end` block, was used.
function make_obs_eqs(observables_expr)
    observables_expr = option_block_form(observables_expr)
    obs_eqs = :([])
    foreach(arg -> push!(obs_eqs.args, arg), observables_expr.args)
    return obs_eqs
end

# Helper function to detect if an expression already contains a `Pre()` call at parse time.
# Used to avoid double-wrapping when users explicitly write `Pre()` in event affects.
function expr_contains_pre(expr)
    if expr isa Expr
        if expr.head == :call && expr.args[1] == :Pre
            return true
        end
        return any(expr_contains_pre, expr.args)
    end
    return false
end

# Read the events (continuous or discrete) provided as options to the DSL. Returns an expression which evaluates to these.
# Infered parameters that are updated byu the event should be declared using e.g. `@discretes p(t)`.
# `read_events_option!` moves these from `ps_inferred` to `discs_inferred`
function read_events_option!(options, discs_inferred::Vector, ps_inferred::Vector, discs_declared::Vector, event_type::Symbol)
    # Prepares the events, if required to, converts them to block form.
    if event_type ∉ [:continuous_events, :discrete_events]
        error("Trying to read an unsupported event type.")
    end
    events_input = haskey(options, event_type) ? get_block_option(options[event_type]) :
        striplines(:(begin end))
    events_input = option_block_form(events_input)

    # Goes through the events, checks for errors, and adds them to the output vector.
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
            error("The condition part of continuous events (the left-hand side) must be a vector. This is not the case for: $(arg).")
        end
        if (arg isa Expr) && (arg.args[3] isa Expr) && (arg.args[3].head != :vect)
            error("The affect part of all events (the right-hand side) must be a vector. This is not the case for: $(arg).")
        end

        # Goes through all affects, checking formatting, recording discrete parameters, and
        # adding `Pre(...)` statements where necessary.
        disc_ps = :([])
        affects = :([])
        for affect in arg.args[3].args
            Meta.isexpr(affect, :call) ||
                error("Event affects must be assignments (e.g. `X ~ X + 1`). This is not the case for: $(affect).")
            (affect.args[2] isa Symbol) ||
                error("The Catalyst DSL currently only supports assignment events where the LHS is a single symbol. This is not the case for: $(affect).")

            # If the event updates an inferred parameter, this should be moved to an inferred discrete.
            if affect.args[2] in ps_inferred
                push!(discs_inferred, affect.args[2])
                deleteat!(ps_inferred, findfirst(==(affect.args[2]), ps_inferred))
            end

            # If the event updates an inferred parameter or declared discrete, it should be in `discrete_parameters`.
            (affect.args[2] in [discs_inferred; discs_declared]) && push!(disc_ps.args, affect.args[2])

            # Creates the affect RHS (adds `Pre` if it doesn't already contain Pre).
            rhs = affect.args[3]
            if !(rhs isa Number) && !expr_contains_pre(rhs)
                rhs = :(Pre($(rhs)))
            end
            push!(affects.args, :($(affect.args[2]) ~ $rhs))
        end

        # Adds the correctly formatted event to the event creation expression.
        event_func = (
            event_type == :continuous_events ? :(MT.SymbolicContinuousCallback) :
                :(MT.SymbolicDiscreteCallback)
        )
        event = :($event_func($(arg.args[2]) => $affects; discrete_parameters = $disc_ps))
        push!(events_expr.args, event)
    end

    return events_expr
end

# Returns the `default_reaction_metadata` output. Technically Catalyst's code could have been made
# more generic to account for other default reaction metadata. Practically, this will likely
# be the only relevant reaction metadata to have a default value via the DSL. If another becomes
# relevant, the code can be rewritten to take this into account.
# Checks if the `@default_noise_scaling` option is used. If so, use it as the default value of
# the `default_noise_scaling` reaction metadata, otherwise, returns an empty vector.
function read_default_noise_scaling_option(options)
    if haskey(options, :default_noise_scaling)
        (length(options[:default_noise_scaling].args) != 3) &&
            error("@default_noise_scaling should only have a single expression as its input, this appears not to be the case: \"$(options[:default_noise_scaling])\"")
        return :([:noise_scaling => $(options[:default_noise_scaling].args[3])])
    end
    return :([])
end

# Reads the combinatorial ratelaw options, which determines if a combinatorial rate law should
# be used or not. If not provided, use the default (true).
function read_combinatoric_ratelaws_option(options)
    return haskey(options, :combinatoric_ratelaws) ?
        get_block_option(options[:combinatoric_ratelaws]) : true
end

### `@reaction` Macro & its Internals ###

"""
    @reaction

Macro for generating a single [`Reaction`](@ref) object using a similar syntax as the `@reaction_network`
macro (but permitting only a single reaction). A more detailed introduction to the syntax can
be found in the description of `@reaction_network`.

The `@reaction` macro is followed by a single line consisting of three parts:
- A rate (at which the reaction occurs).
- Any number of substrates (which are consumed by the reaction).
- Any number of products (which are produced by the reaction).

The output is a reaction (just like created using the `Reaction` constructor).

Examples:
Here we create a simple binding reaction and store it in the variable rx:
```julia
rx = @reaction k, X + Y --> XY
```
The macro will automatically deduce `X`, `Y`, and `XY` to be species (as these occur as reactants)
and `k` as a parameter (as it does not occur as a reactant).

The `@reaction` macro provides a more concise notation to the `Reaction` constructor. I.e. here
we create the same reaction using both approaches, and also confirm that they are identical.
```julia
# Creates a reaction using the `@reaction` macro.
rx = @reaction k*v, A + B --> C + D

# Creates a reaction using the `Reaction` constructor.
t = default_t()
@parameters k v
@species A(t) B(t) C(t) D(t)
rx2 = Reaction(k*v, [A, B], [C, D])

# Confirms that the two approaches yield identical results:
rx1 == rx2
```
Interpolation of already declared symbolic variables into `@reaction` is possible:
```julia
t = default_t()
@parameters k b
@species A(t)
ex = k*A^2 + t
rx = @reaction b*\$ex*\$A, \$A --> C
```

Notes:
- `@reaction` does not support bi-directional type reactions (using `<-->`) or reaction bundling
(e.g. `d, (X,Y) --> 0`).
- Interpolation of Julia variables into the macro works similarly to the `@reaction_network`
macro. See [The Reaction DSL](@ref dsl_description) tutorial for more details.
"""
macro reaction(ex)
    return make_reaction(ex)
end

# Function for creating a Reaction structure (used by the @reaction macro).
function make_reaction(ex::Expr)

    # Handle interpolation of variables in the input.
    ex = esc_dollars!(ex)

    # Parses reactions. Extracts species and parameters within it.
    reaction = get_reaction(ex)
    species, parameters, stoich_ps = extract_sps_and_ps([reaction], [])

    # Checks for input errors. Needed here but not in `@reaction_network` as `ReactionSystem` performs this check but `Reaction` doesn't.
    forbidden_symbol_check(union(species, parameters))

    # Creates expressions corresponding to code for declaring the parameters, species, and reaction.
    spexprs = get_usexpr(species, Dict{Symbol, Expr}())
    pexprs = get_psexpr(parameters, stoich_ps, Dict{Symbol, Expr}())
    rxexpr = get_rxexpr(reaction)
    iv = :($(DEFAULT_IV_SYM) = default_t())

    # Returns a rephrased expression which generates the `Reaction`.
    return quote
        $pexprs
        $iv
        $spexprs
        $rxexpr
    end
end

# Reads a single line and creates the corresponding DSLReaction.
function get_reaction(line)
    reaction = get_reactions([line])
    (length(reaction) != 1) &&
        error("Malformed reaction. @reaction macro only creates a single reaction. E.g. double arrows, such as `<-->` are not supported.")
    return only(reaction)
end

### Generic Expression Manipulation ###

# Recursively traverses an expression and escapes all the user-defined functions.
# Special function calls like "hill(...)" are not expanded.
function recursive_escape_functions!(expr::ExprValues, syms_skip = [])
    (typeof(expr) != Expr) && (return expr)
    foreach(
        i -> expr.args[i] = recursive_escape_functions!(expr.args[i], syms_skip),
        1:length(expr.args)
    )
    if (expr.head == :call) && (expr.args[1] isa Symbol) &&!isdefined(Catalyst, expr.args[1]) &&
            expr.args[1] ∉ syms_skip
        expr.args[1] = esc(expr.args[1])
    end
    return expr
end

# Returns the length of an expression tuple, or 1 if it is not an expression tuple (probably a Symbol/Numerical).
function tup_leng(ex::ExprValues)
    (typeof(ex) == Expr && ex.head == :tuple) && (return length(ex.args))
    return 1
end

# Gets the ith element in an expression tuple, or returns the input itself if it is not an expression tuple
# (probably a  Symbol/Numerical). This is used to handle bundled reactions (like `d, (X,Y) --> 0`).
function get_tup_arg(ex::ExprValues, i::Int)
    (tup_leng(ex) == 1) && (return ex)
    return ex.args[i]
end
