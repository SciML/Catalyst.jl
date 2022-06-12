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
    end kB, kD

    ### Defining New Functions ###
    my_hill_repression(x, v, k, n) = v*k^n/(k^n+x^n)

    # may be necessary to
    # @register_symbolic my_hill_repression(x, v, k, n)
    # see https://mtk.sciml.ai/stable/tutorials/symbolic_functions/#Registering-Functions-1

    r = @reaction_network MyReactionType begin
        my_hill_repression(x, v_x, k_x, n_x), 0 --> x
    end v_x k_x n_x

    ### Simulating Reaction Networks ###
    probODE = ODEProblem(rn, args...; kwargs...)        # Using multiple dispatch the reaction network can be used as input to create ODE, SDE and Jump problems.
    probSDE = SDEProblem(rn, args...; kwargs...)
    probJump = JumpProblem(prob,aggregator::Direct,rn)
"""

# Declare various arrow types symbols used for the empty set (also 0).
const empty_set = Set{Symbol}([:∅])
const fwd_arrows = Set{Symbol}([:>, :→, :↣, :↦, :⇾, :⟶, :⟼, :⥟, :⥟, :⇀, :⇁, :⇒, :⟾])
const bwd_arrows = Set{Symbol}([:<, :←, :↢, :↤, :⇽, :⟵, :⟻, :⥚, :⥞, :↼, :↽, :⇐, :⟽, Symbol("<--")])
const double_arrows = Set{Symbol}([:↔, :⟷, :⇄, :⇆, :⇌, :⇋, :⇔, :⟺, Symbol("<-->")])
const pure_rate_arrows = Set{Symbol}([:⇐, :⟽, :⇒, :⟾, :⇔, :⟺])

# Declares symbols which may neither be used as parameters not varriables.
forbidden_symbols = [:t, :π, :pi, :ℯ, :im, :nothing, :∅]


### The main macro, takes reaction network notation and returns a ReactionSystem. ###
"""
    @reaction_network

Generates a [`ReactionSystem`](@ref) that encodes a chemical reaction
network.

See [The Reaction DSL](@ref) documentation for details on
parameters to the macro.

Examples:
```julia
# a basic SIR model, with name SIR
sir_model = @reaction_network SIR begin
    c1, s + i --> 2i
    c2, i --> r
end c1 c2

# a basic SIR model, with random generated name
sir_model = @reaction_network begin
    c1, s + i --> 2i
    c2, i --> r
end c1 c2

# an empty network with name empty
emptyrn = @reaction_network empty

# an empty network with random generated name
emptyrn = @reaction_network
```
"""
macro reaction_network(name::Symbol, ex::Expr, parameters...)
    make_reaction_system(MacroTools.striplines(ex), parameters; name=:($(QuoteNode(name))))
end

# allows @reaction_network $name begin ... to interpolate variables storing a name
macro reaction_network(name::Expr, ex::Expr, parameters...)
    make_reaction_system(MacroTools.striplines(ex), parameters; name=:($(esc(name.args[1]))))
end

macro reaction_network(ex::Expr, parameters...)
    ex = MacroTools.striplines(ex)

    # no name but equations: @reaction_network begin ... end ...
    if ex.head == :block
        make_reaction_system(ex, parameters)
    else  # empty but has interpolated name: @reaction_network $name
        networkname = :($(esc(ex.args[1])))
        return Expr(:block,:(@parameters t),
                    :(ReactionSystem(Reaction[],
                                    t,
                                    [],
                                    [];
                                    name=$networkname)))
    end
end

#Returns a empty network (with, or without, a declared name)
# @reaction_network name
macro reaction_network(name::Symbol=gensym(:ReactionSystem))
    return Expr(:block,:(@parameters t),
                :(ReactionSystem(Reaction[],
                                 t,
                                 [],
                                 [];
                                 name=$(QuoteNode(name)))))
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
@variables t A(t) B(t) C(t) D(t)
rx == Reaction(k*v, [A,B], [C,D])
```
Here `k` and `v` will be parameters and `A`, `B`, `C` and `D` will be variables.
Interpolation of existing parameters/variables also works
```julia
@parameters k b
@variables t A(t)
ex = k*A^2 + t
rx = @reaction b*$ex*$A, $A --> C
```

Notes:
- Any symbols arising in the rate expression that aren't interpolated are treated as
  parameters. In the reaction part (`α*A + B --> C + D`), coefficients are treated as
  parameters, e.g. `α`, and rightmost symbols as species, e.g. `A,B,C,D`.
- Works with any *single* arrow types supported by [`@reaction_network`](@ref).
- Interpolation of Julia variables into the macro works similar to the `@reaction_network`
  macro. See [The Reaction DSL](@ref) tutorial for more details.
"""
macro reaction(ex)
    make_reaction(ex)
end

"""
    @add_reactions

Adds the reactions declared to a preexisting [`ReactionSystem`](@ref). All
parameters used in the added reactions need to be declared after the
reactions.

See the [Catalyst.jl for Reaction Network Modeling](@ref) documentation for details on
parameters to the macro.
"""
macro add_reactions(rn::Symbol, ex::Expr, parameters...)
    :(merge!($(esc(rn)),$(make_reaction_system(MacroTools.striplines(ex), parameters))))
end



### Sturctures containing information about a reactant, and a reaction, respectively.

#Structure containing information about one reactant in one reaction.
struct ReactantStruct
    reactant::Union{Symbol,Expr}
    stoichiometry::ExprValues
end
#Structure containing information about one Reaction. Contain all its substrates and products as well as its rate. Contains a specialized constructor.
struct ReactionStruct
    substrates::Vector{ReactantStruct}
    products::Vector{ReactantStruct}
    rate::ExprValues
    only_use_rate::Bool

    function ReactionStruct(sub_line::ExprValues, prod_line::ExprValues, rate::ExprValues, only_use_rate::Bool)
        sub = recursive_find_reactants!(sub_line,1,Vector{ReactantStruct}(undef,0))
        prod = recursive_find_reactants!(prod_line,1,Vector{ReactantStruct}(undef,0))
        new(sub, prod, rate, only_use_rate)
    end
end


### Functions that process the input and rephrase it as a reaction system ###
function esc_dollars!(ex)
    if ex isa Expr
        if ex.head == :$
            return esc(:($(ex.args[1])))
        else
            for i = 1:length(ex.args)
                ex.args[i] = esc_dollars!(ex.args[i])
            end
        end
    end
    ex
end

function get_rxexprs(rxstruct)
    subs_init = isempty(rxstruct.substrates) ? nothing : :([]); subs_stoich_init = deepcopy(subs_init)
    prod_init = isempty(rxstruct.products) ? nothing : :([]); prod_stoich_init = deepcopy(prod_init)
    reaction_func = :(Reaction($(recursive_expand_functions!(rxstruct.rate)), $subs_init, $prod_init, $subs_stoich_init, $prod_stoich_init, only_use_rate=$(rxstruct.only_use_rate)))
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

function get_pexprs(psyms)
    pexprs = isempty(psyms) ? :() : :(@parameters)
    if !isempty(psyms)
        foreach(psym-> push!(pexprs.args, psym), psyms)
    end
    pexprs
end

function get_sexprs(ssyms)
    sexprs = :(@variables t)
    foreach(s -> (s isa Symbol) && push!(sexprs.args, Expr(:call,s,:t)), ssyms)
    sexprs
end

# Takes the reactions, and rephrases it as a "ReactionSystem" call, as designated by the ModelingToolkit IR.
function make_reaction_system(ex::Expr, parameters; name=:(gensym(:ReactionSystem)))

    # handle interpolation of variables
    ex = esc_dollars!(ex)

    # parse DSL lines
    reactions = get_reactions(ex)
    reactants = get_reactants(reactions)
    allspecies = union(reactants, get_rx_species(reactions,parameters))
    !isempty(intersect(forbidden_symbols,union(allspecies,parameters))) &&
        error("The following symbol(s) are used as species or parameters: "*((map(s -> "'"*string(s)*"', ",intersect(forbidden_symbols,union(species,parameters)))...))*"this is not permited.")

    pexprs = get_pexprs(parameters)    # parameters
    sexprs = get_sexprs(allspecies)    # species

    # ReactionSystem
    rxexprs = :($(make_ReactionSystem_internal)([],t,nothing,[]; name=$(name)))
    foreach(parameter -> push!(rxexprs.args[6].args,parameter), parameters)
    for reaction in reactions
        push!(rxexprs.args[3].args, get_rxexprs(reaction))
    end

    quote
        $pexprs
        $sexprs
        $rxexprs
    end
end

function make_reaction(ex::Expr)

    # handle interpolation of variables
    ex = esc_dollars!(ex)

    # parse DSL lines
    reaction   = get_reaction(ex)
    allspecies = get_reactants(reaction)   # species defined by stoich
    parameters = get_rx_species([reaction],Symbol[])  # anything in a rate is a parameter
    !isempty(intersect(forbidden_symbols,union(allspecies,parameters))) &&
        error("The following symbol(s) are used as species or parameters: "*((map(s -> "'"*string(s)*"', ",intersect(forbidden_symbols,union(species,parameters)))...))*"this is not permited.")

    pexprs = get_pexprs(parameters)    # parameters
    sexprs = get_sexprs(allspecies)    # species
    rxexpr = get_rxexprs(reaction)     # reaction

    quote
        $pexprs
        $sexprs
        $rxexpr
    end
end

function get_rx_species(rxs, ps)
    pset = Set(ps)
    species_set = Set{Symbol}()
    for rx in rxs
        find_species_in_rate!(species_set, rx.rate, pset)
        for sub in rx.substrates
            find_species_in_rate!(species_set, sub.stoichiometry, pset)
        end
        for prod in rx.products
            find_species_in_rate!(species_set, prod.stoichiometry, pset)
        end
    end
    collect(species_set)
end

function find_species_in_rate!(sset, rateex::ExprValues, ps)
    if rateex isa Symbol
        if !(rateex in forbidden_symbols) && !(rateex in ps)
            push!(sset, rateex)
        end
    elseif rateex isa Expr
        # note, this (correctly) skips $(...) expressions
        for i = 2:length(rateex.args)
            find_species_in_rate!(sset, rateex.args[i], ps)
        end
    end
    nothing
end

function get_reaction(line)
    (line.head != :tuple) && error("Malformed reaction line: $(MacroTools.striplines(line))")
    (rate,r_line) = line.args
    (r_line.head  == :-->) && (r_line = Expr(:call,:→,r_line.args[1],r_line.args[2]))

    arrow = r_line.args[1]
    in(arrow,double_arrows) && error("Double arrows not allowed for single reactions.")

    only_use_rate = in(arrow,pure_rate_arrows)
    if in(arrow,fwd_arrows)
        rs = create_ReactionStruct(r_line.args[2], r_line.args[3], rate, only_use_rate)
    elseif in(arrow,bwd_arrows)
        rs = create_ReactionStruct(r_line.args[3], r_line.args[2], rate, only_use_rate)
    else
        throw("Malformed reaction, invalid arrow type used in: $(MacroTools.striplines(line))")
    end

    rs
end

#Generates a vector containing a number of reaction structures, each containing the information about one reaction.
function get_reactions(ex::Expr, reactions = Vector{ReactionStruct}(undef,0))
    for line in ex.args
        (line.head != :tuple) && error("Malformed reaction line: $(MacroTools.striplines(line))")
        (rate,r_line) = line.args
        (r_line.head  == :-->) && (r_line = Expr(:call,:→,r_line.args[1],r_line.args[2]))

        arrow = r_line.args[1]
        only_use_rate = in(arrow,pure_rate_arrows)
        if in(arrow,double_arrows)
            (typeof(rate) == Expr && rate.head == :tuple) || error("Error: Must provide a tuple of reaction rates when declaring a bi-directional reaction.")
            push_reactions!(reactions, r_line.args[2], r_line.args[3], rate.args[1], only_use_rate)
            push_reactions!(reactions, r_line.args[3], r_line.args[2], rate.args[2], only_use_rate)
        elseif in(arrow,fwd_arrows)
            push_reactions!(reactions, r_line.args[2], r_line.args[3], rate, only_use_rate)
        elseif in(arrow,bwd_arrows)
            push_reactions!(reactions, r_line.args[3], r_line.args[2], rate, only_use_rate)
        else
            throw("Malformed reaction, invalid arrow type used in: $(MacroTools.striplines(line))")
        end
    end
    return reactions
end

function create_ReactionStruct(sub_line::ExprValues, prod_line::ExprValues, rate::ExprValues, only_use_rate::Bool)
    all(==(1), (tup_leng(sub_line), tup_leng(prod_line), tup_leng(rate))) ||
        error("Malformed reaction, line appears to be defining multiple reactions incorrectly: rate=$rate, subs=$sub_line, prods=$prod_line.")
    ReactionStruct(get_tup_arg(sub_line,1), get_tup_arg(prod_line,1), get_tup_arg(rate,1), only_use_rate)
end

#Takes a reaction line and creates reactions from it and pushes those to the reaction array. Used to create multiple reactions from, for instance, 1.0, (X,Y) --> 0.
function push_reactions!(reactions::Vector{ReactionStruct}, sub_line::ExprValues, prod_line::ExprValues, rate::ExprValues, only_use_rate::Bool)
    lengs = (tup_leng(sub_line), tup_leng(prod_line), tup_leng(rate))
    for i = 1:maximum(lengs)
        (count(lengs.==1) + count(lengs.==maximum(lengs)) < 3) && (throw("Malformed reaction, rate=$rate, subs=$sub_line, prods=$prod_line."))
        push!(reactions, ReactionStruct(get_tup_arg(sub_line,i), get_tup_arg(prod_line,i), get_tup_arg(rate,i), only_use_rate))
    end
end

function processmult(op, mult, stoich)
    if (mult isa Number) && (stoich isa Number)
        op(mult, stoich)
    else
        :($op($mult,$stoich))
    end
end

#Recursive function that loops through the reaction line and finds the reactants and their stoichiometry. Recursion makes it able to handle weird cases like 2(X+Y+3(Z+XY)).
function recursive_find_reactants!(ex::ExprValues, mult::ExprValues, reactants::Vector{ReactantStruct})
    if typeof(ex)!=Expr || (ex.head == :escape)
        (ex == 0 || in(ex,empty_set)) && (return reactants)
        if any(ex==reactant.reactant for reactant in reactants)
            idx = findall(x -> x==ex, getfield.(reactants,:reactant))[1]
            reactants[idx] = ReactantStruct(ex,processmult(+,mult,reactants[idx].stoichiometry))
        else
            push!(reactants, ReactantStruct(ex,mult))
        end
    elseif ex.args[1] == :*
        if length(ex.args) == 3
            recursive_find_reactants!(ex.args[3],processmult(*,mult,ex.args[2]),reactants)
        else
            newmult = processmult(*, mult, Expr(:call,ex.args[1:end-1]...))
            recursive_find_reactants!(ex.args[end],newmult,reactants)
        end
    elseif ex.args[1] == :+
        for i = 2:length(ex.args)
            recursive_find_reactants!(ex.args[i],mult,reactants)
        end
    else
        throw("Malformed reaction, bad operator: $(ex.args[1]) found in stochiometry expression $ex.")
    end
    return reactants
end


function get_reactants(reaction::ReactionStruct, reactants=Vector{Union{Symbol,Expr}}())
    for reactant in Iterators.flatten((reaction.substrates,reaction.products))
        !in(reactant.reactant,reactants) && push!(reactants,reactant.reactant)
    end
    return reactants
end

# Extract the reactants from the set of reactions.
function get_reactants(reactions::Vector{ReactionStruct}, reactants=Vector{Union{Symbol,Expr}}())
    for reaction in reactions
        get_reactants(reaction, reactants)
    end
    return reactants
end

### Functionality for expanding function call to custom and specific functions ###

#Recursively traverses an expression and replaces special function call like "hill(...)" with the actual corresponding expression.
function recursive_expand_functions!(expr::ExprValues)
    (typeof(expr)!=Expr) && (return expr)
    foreach(i -> expr.args[i] = recursive_expand_functions!(expr.args[i]), 1:length(expr.args))
    if expr.head == :call
        !isdefined(Catalyst,expr.args[1]) && (expr.args[1] = esc(expr.args[1]))
    end
    return expr
end
