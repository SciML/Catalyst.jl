"""
Macro that inputs an expression corresponding to a reaction network and output a Reaction Network Structure that can be used as input to generation of SDE and ODE and Jump problems.
Most arrows accepted (both right, left and bi drectional arrows).
Note that while --> is a correct arrow, neither <-- nor <--> works.
Using non-filled arrows (⇐, ⟽, ⇒, ⟾, ⇔, ⟺) will disable mass kinetics and lets you cutomize reaction rates yourself.
Use 0 or ∅ for degradation/creation to/from nothing.
Example systems:
    ### Basic Usage ###
    rn = @reaction_network rType begin #Creates a reaction network of type rType.
        2.0, X + Y --> XY                  #This will have reaction rate corresponding to 2.0*[X][Y]
        2.0, XY ← X + Y                    #Identical to 2.0, X + Y --> XY
    end

    ### Manipulating Reaction Rates ###
    rn = @reaction_network rType begin
        2.0, X + Y ⟾ XY                   #Ignores mass kinetics. This will have reaction rate corresponding to 2.0.
        2.0X, X + Y --> XY                 #Reaction rate needs not be constant. This will have reaction rate corresponding to 2.0*[X]*[X]*[Y].
        XY+log(X)^2, X + Y --> XY          #Reaction rate accepts quite complicated expressions (user defined functions must first be registered using the @reaction_func macro).
        hill(XY,2,2,2), X + Y --> XY       #Reaction inis activated by XY according to a hill function. hill(x,v,K,N).
        mm(XY,2,2), X + Y --> XY           #Reaction inis activated by XY according to a michaelis menten function. mm(x,v,K).
    end

    ### Multipple Reactions on a Single Line ###
    rn = @reaction_network rType begin
        (2.0,1.0), X + Y ↔ XY              #Identical to reactions (2.0, X + Y --> XY) and (1.0, XY --> X + Y).
        2.0, (X,Y) --> 0                   #This corresponds to both X and Y degrading at rate 2.0.
        (2.0, 1.0), (X,Y) --> 0            #This corresponds to X and Y degrading at rates 2.0 and 1.0, respectively.
        2.0, (X1,Y1) --> (X2,Y2)           #X1 and Y1 becomes X2 and Y2, respectively, at rate 2.0.
    end

    ### Adding Parameters ###
    kB = 2.0; kD = 1.0
    p = [kB, kD]
    p = []
    rn = @reaction_network type begin
        (kB, kD), X + Y ↔ XY            #Lets you define parameters outside on network. Parameters can be changed without recalling the network.
    end kB, kD

    ### Defining New Functions ###
    @reaction_func my_hill_repression(x, v, k, n) = v*k^n/(k^n+x^n)     #Creates and adds a new function that the @reaction_network macro can see.
    r = @reaction_network MyReactionType begin
        my_hill_repression(x, v_x, k_x, n_x), 0 --> x                       #After it has been added in @reaction_func the function can be used when defining new reaction networks.
    end v_x k_x n_x

    ### Simulating Reaction Networks ###
    probODE = ODEProblem(rn, args...; kwargs...)        #Using multiple dispatch the reaction network can be used as input to create ODE, SDE and Jump problems.
    probSDE = SDEProblem(rn, args...; kwargs...)
    probJump = JumpProblem(prob,aggregator::Direct,rn)
"""

"""
    @reaction_network

Generates a subtype of an `AbstractReactionNetwork` that encodes a chemical
reaction network, and complete ODE, SDE and jump representations of the system.
See the [Chemical Reaction Model
docs](http://docs.juliadiffeq.org/latest/models/biological.html) for details on
parameters to the macro.
"""
macro reaction_network(name, ex::Expr, p...)
    coordinate(name, MacroTools.striplines(ex), p, :no___noise___scaling)
end

#Macro to create a reaction network model. Multiple dispatch is used to allow for SDE noise scalling.
macro reaction_network(name, scale_noise, ex::Expr, p...)
    in(scale_noise, p) || (p = (p..., scale_noise))
    coordinate(name, MacroTools.striplines(ex), p, scale_noise)
end

#If no type name is given, creates a network with a default name.
macro reaction_network(ex::Expr, p...)
    coordinate(:reaction_network, MacroTools.striplines(ex), p, :no___noise___scaling)
end

################# query-based macros:

"""
    @min_reaction_network

Generates a subtype of an `AbstractReactionNetwork` that only encodes a chemical
reaction network. Use [`addodes!`](@ref), [`addsdes!`](@ref) or
[`addjumps!`](@ref) to complete the network for specific problem types. It accepts
the same arguments as [`@reaction_network`](@ref).
"""
macro min_reaction_network(name, ex::Expr, p...)
    min_coordinate(name, MacroTools.striplines(ex), p, :no___noise___scaling)
end

#Macro to create a reaction network model. Multiple dispatch is used to allow for SDE noise scalling.
macro min_reaction_network(name, scale_noise, ex::Expr, p...)
    in(scale_noise, p) || (p = (p..., scale_noise))
    min_coordinate(name, MacroTools.striplines(ex), p, scale_noise)
end

#If no type name is given, creates a network with a default name.
macro min_reaction_network(ex::Expr, p...)
    min_coordinate(:min_reaction_network, MacroTools.striplines(ex), p, :no___noise___scaling)
end

"""
    @empty_reaction_network networktype

Generates a subtype of an `AbstractReactionNetwork` that encodes an empty
chemical reaction network. `networktype` is an optional parameter that specifies
the type of the generated network. Use [`addspecies!`](@ref), [`addparam!`](@ref)
and [`addreaction!`](@ref) to extend the network.  Use [`addodes!`](@ref),
[`addsdes!`](@ref) or [`addjumps!`](@ref) to complete the network for specific
problem types.
"""
macro empty_reaction_network(name::Symbol=:min_reaction_network)
    min_coordinate(name, MacroTools.striplines(:(begin end)), (), :no___noise___scaling)
end

#################

#Declare various arrow types symbols used for the empty set (also 0).
empty_set = Set{Symbol}([:∅])
fwd_arrows = Set{Symbol}([:>, :→, :↣, :↦, :⇾, :⟶, :⟼, :⥟, :⥟, :⇀, :⇁, :⇒, :⟾])
bwd_arrows = Set{Symbol}([:<, :←, :↢, :↤, :⇽, :⟵, :⟻, :⥚, :⥞, :↼, :↽, :⇐, :⟽])
double_arrows = Set{Symbol}([:↔, :⟷, :⇄, :⇆, :⇔, :⟺])
no_mass_arrows = Set{Symbol}([:⇐, :⟽, :⇒, :⟾, :⇔, :⟺])      #Using this arrows will disable the program from multiplying reaction rates with the substrate concentrations. Gives user full control of reaction rates.

funcdict = Dict{Symbol, Function}()                             #Stores user defined functions.

#Coordination function, actually does all the work of the macro.
function coordinate(name, ex::Expr, p, scale_noise)

    # minimal reaction network components
    (reactions, reactants, parameters, syms, params) = get_minnetwork(ex, p)

    # expressions for ODEs
    (f_expr, f, f_rhs, symjac, jac, paramjac, f_symfuncs, jac_prototype) = genode_exprs(reactions, reactants, parameters, syms)
    odefun = :(ODEFunction(f; jac=$jac, jac_prototype=$jac_prototype, paramjac=$paramjac, syms=$syms))

    # expressions for SDEs
    (g_expr, g, g_funcs, p_matrix) = gensde_exprs(reactions, reactants, parameters, scale_noise)
    sdefun = :(SDEFunction(f, g; jac=$jac, jac_prototype=$jac_prototype, paramjac=$paramjac, syms=$syms))

    # expressions for jumps
    (jump_rate_expr, jump_affect_expr, jumps) = get_jumps(reactions, reactants, parameters)

    # expressions for regular jumps
    regular_jumps = get_regularjumps(reactions, reactants, parameters)

    # Build the type
    exprs = Vector{Expr}(undef,0)
    typeex,constructorex = maketype(DiffEqBase.AbstractReactionNetwork, name, f, f_rhs, f_symfuncs, g, g_funcs, jumps, regular_jumps, Meta.quot(jump_rate_expr), Meta.quot(jump_affect_expr), p_matrix, syms, scale_noise; params=params, reactions=reactions, jac=jac, paramjac=paramjac, symjac=symjac, syms_to_ints=reactants, params_to_ints=parameters, odefun=odefun, sdefun=sdefun)
    push!(exprs,typeex)
    push!(exprs,constructorex)

    # add type functions
    append!(exprs, gentypefun_exprs(name))

    # return as one expression block
    expr_arr_to_block(exprs)
end

# min_reaction_network coordination function, actually does all the work of the macro.
function min_coordinate(name, ex::Expr, p, scale_noise)

    # minimal reaction network components
    (reactions, reactants, parameters, syms, params) = get_minnetwork(ex, p)

    # Build the type
    exprs = Vector{Expr}(undef,0)
    typeex,constructorex = maketype(DiffEqBase.AbstractReactionNetwork, name, nothing, nothing,
                                    nothing, nothing, nothing, nothing, nothing, nothing,
                                    nothing, nothing, syms, scale_noise; params=params,
                                    reactions=reactions, symjac=nothing, syms_to_ints=reactants,
                                    params_to_ints=parameters)
    push!(exprs,typeex)
    push!(exprs,constructorex)

    # add type functions
    append!(exprs, gentypefun_exprs(name, gen_inplace=false, gen_outofplace=false))

    # return as one expression block
    expr_arr_to_block(exprs)
end

# ODE expressions
function genode_exprs(reactions, reactants, parameters, syms; build_jac=true,
                                                              build_paramjac=true,
                                                              sparse_jac=false,
                                                              build_symfuncs=true,
                                                              zeroout_jac=false)
    f_expr     = get_f(reactions, reactants)
    f          = make_func(f_expr, reactants, parameters)
    f_rhs      = ExprValues[element.args[2] for element in f_expr]
    f_symfuncs = build_symfuncs ? hcat([SymEngine.Basic(f) for f in f_rhs]) : nothing
    symjac, jac, jac_protoype, paramjac = get_jacs(f_rhs, syms, reactions, reactants, parameters; build_jac=build_jac, sparse_jac=sparse_jac, build_paramjac=build_paramjac, zeroout_jac=zeroout_jac)

    (f_expr,f,f_rhs,symjac,jac,paramjac,f_symfuncs, jac_protoype)
end

# SDE expressions
function gensde_exprs(reactions, reactants, parameters, scale_noise)
    g_expr   = get_g(reactions, reactants, scale_noise)
    g        = make_func(g_expr, reactants, parameters)
    g_funcs  = ExprValues[element.args[2] for element in g_expr]
    p_matrix = zeros(length(reactants), length(reactions))

    (g_expr,g,g_funcs,p_matrix)
end

# generate the minimal network components
function get_minnetwork(ex::Expr, p)
    reactions = get_reactions(ex)           ::Vector{ReactionStruct}
    reactants = get_reactants(reactions)    ::OrderedDict{Symbol,Int}
    parameters = get_parameters(p)          ::OrderedDict{Symbol,Int}

    syms = collect(keys(reactants))
    params = collect(keys(parameters))
    (in(:t,union(syms,params))) && error("t is reserved for the time variable and may neither be used as a reactant nor a parameter")

    update_reaction_info(reactions,syms)

    (reactions,reactants,parameters,syms,params)
end

#Generates a vector containing a number of reaction structures, each containing the infromation about one reaction.
function get_reactions(ex::Expr, reactions = Vector{ReactionStruct}(undef,0))
    for line in ex.args
        (line.head != :tuple) && (continue)
        (rate,r_line) = line.args
        if r_line.head  == :-->
            r_line = Expr(:call,:→,r_line.args[1],r_line.args[2])
        end

        arrow = r_line.args[1]  ::Symbol
        if in(arrow,double_arrows)
            (typeof(rate) == Expr && rate.head == :tuple) || error("Error: Must provide a tuple of reaction rates when declaring a bi-directional reaction.")
            push_reactions(reactions, r_line.args[2], r_line.args[3], rate.args[1], !in(arrow,no_mass_arrows))
            push_reactions(reactions, r_line.args[3], r_line.args[2], rate.args[2], !in(arrow,no_mass_arrows))
        elseif in(arrow,fwd_arrows)
            push_reactions(reactions, r_line.args[2], r_line.args[3], rate, !in(arrow,no_mass_arrows))
        elseif in(arrow,bwd_arrows)
            push_reactions(reactions, r_line.args[3], r_line.args[2], rate, !in(arrow,no_mass_arrows))
        else
            throw("malformed reaction")
        end
    end
    return reactions
end

#Structure containing information about one reactant in one reaction.
struct ReactantStruct
    reactant::Symbol
    stoichiometry::Int
end

#Structure containing information about one Reaction. Contain all its substrates and products as well as its rate. Contains an specialized constructor.
struct ReactionStruct
    substrates::Vector{ReactantStruct}
    products::Vector{ReactantStruct}
    netstoich::Vector{ReactantStruct}
    rate_org::ExprValues
    rate_DE::ExprValues
    rate_SSA::ExprValues
    dependants::Vector{Symbol}
    is_pure_mass_action::Bool

    function ReactionStruct(s::Vector{ReactantStruct}, p::Vector{ReactantStruct}, 
                            ns::Vector{ReactantStruct},
                            ro::ExprValues, rde::ExprValues, rssa::ExprValues,
                            dep::Vector{Symbol}, isma::Bool)
        new(s,p,ns,ro,rde,rssa,dep,isma)
    end

    function ReactionStruct(sub_line::ExprValues, prod_line::ExprValues, rate::ExprValues, use_mass_kin::Bool)
        sub = add_reactants!(sub_line,1,Vector{ReactantStruct}(undef,0))
        prod = add_reactants!(prod_line,1,Vector{ReactantStruct}(undef,0))
        ns = netstoich(sub, prod)

        rate_DE = isempty(sub) ? rate : mass_rate_DE(sub, use_mass_kin, rate)
        rate_SSA =  isempty(sub) ? rate : mass_rate_SSA(sub, use_mass_kin, rate)
        new(sub, prod, ns, rate, rate_DE, rate_SSA, [], use_mass_kin)
    end
    function ReactionStruct(r::ReactionStruct, syms::Vector{Symbol})
        deps = unique!(recursive_content(r.rate_DE,syms,Vector{Symbol}()))
        is_ma = r.is_pure_mass_action && (length(recursive_content(r.rate_org,syms,Vector{Symbol}()))==0)
        new(r.substrates, r.products, r.netstoich, r.rate_org, r.rate_DE, r.rate_SSA, deps, is_ma)
    end
end

# calculates the net stoichiometry of a reaction
function netstoich(substrates::Vector{ReactantStruct}, products::Vector{ReactantStruct})
    nsdict = Dict{Symbol,Int}(s.reactant => -s.stoichiometry for s in substrates)

    for prod in products
        rsym = prod.reactant
        rcoef = prod.stoichiometry
        @inbounds nsdict[rsym] = haskey(nsdict, rsym) ? nsdict[rsym] + rcoef : rcoef
    end

    net_stoich = Vector{ReactantStruct}()
    for (specsym,coef) in nsdict
        (coef != zero(coef)) && push!(net_stoich, ReactantStruct(specsym,coef))
    end

    net_stoich
end

#Calculates the rate used by ODEs and SDEs. If we want to use masskinetics we have to include substrate concentration, taking higher order terms into account.
function mass_rate_DE(substrates::Vector{ReactantStruct}, use_mass_kin::Bool, old_rate::ExprValues)
    rate = Expr(:call, :*, old_rate)
    coef = one(typeof(substrates[1].stoichiometry))
    if use_mass_kin
        for sub in substrates
            stoich = sub.stoichiometry
            coef *= factorial(stoich)
            push!(rate.args, stoich > one(stoich) ? Expr(:call, :^, sub.reactant, stoich) : sub.reactant)
        end
        (coef > one(coef)) && (rate.args[2] = :($old_rate / $coef))
    end
    return rate
end

#Calculates the rate used by SSAs. If we want to use masskinetics we have to include substrate concentration, taking higher order terms into account.
function mass_rate_SSA(substrates::Vector{ReactantStruct}, use_mass_kin::Bool, old_rate::ExprValues)
    rate = Expr(:call, :*, old_rate)
    use_mass_kin && foreach(sub -> push!(rate.args, :(binomial($(sub.reactant),$(sub.stoichiometry)))), substrates)
    return rate
end

#Takes a reaction line and creates reactions from it and pushes those to the reaction array. Used to creat multiple reactions from e.g. 1.0, (X,Y) --> 0.
function push_reactions(reactions::Vector{ReactionStruct}, sub_line::ExprValues, prod_line::ExprValues, rate::ExprValues, use_mass_kin::Bool)
    lengs = [tup_leng(sub_line), tup_leng(prod_line), tup_leng(rate)]
    (count(lengs.==1) + count(lengs.==maximum(lengs)) < 3) && (throw("malformed reaction"))
    for i = 1:maximum(lengs)
        push!(reactions, ReactionStruct(get_tup_arg(sub_line,i), get_tup_arg(prod_line,i), get_tup_arg(rate,i), use_mass_kin))
    end
end

#Recursive function that loops through the reactants in an reaction line and finds the reactants and their stoichiometry. Recursion makes it able to handle e.g. 2(X+Y+3(Z+XY)) (probably one will not need it though).
function add_reactants!(ex::ExprValues, mult::Int, reactants::Vector{ReactantStruct})
    if typeof(ex)!=Expr
        (ex == 0 || in(ex,empty_set)) && (return reactants)
        if in(ex, getfield.(reactants,:reactant))
            idx = findall(x -> x==ex ,getfield.(reactants,:reactant))[1]
            reactants[idx] = ReactantStruct(ex,mult+reactants[idx].stoichiometry)
        else
            push!(reactants, ReactantStruct(ex,mult))
        end
    elseif ex.args[1] == :*
        add_reactants!(ex.args[3],mult*ex.args[2],reactants)
    elseif ex.args[1] == :+
        for i = 2:length(ex.args)
            add_reactants!(ex.args[i],mult,reactants)
        end
    else
        throw("malformed reaction")
    end
    return reactants
end

#For each reaction, sets its dependencies and whenever it is a pure mass action reaction.
function update_reaction_info(reactions::Vector{ReactionStruct},syms::Vector{Symbol})
    for i = 1:length(reactions)
        reactions[i] = ReactionStruct(reactions[i],syms)
    end
end

#From the vector with all reactions, generates a dictionary with all reactants. Each reactant will point to a number so that X --> means X will be replaced with u[1] in the equations.
function get_reactants(reactions::Vector{ReactionStruct})
    reactants = OrderedDict{Symbol,Int}()
    r_count = 0    ::Int
    for reaction in reactions
        for sub in reaction.substrates
            (!haskey(reactants,sub.reactant)) && (reactants[sub.reactant] = r_count += 1)
        end
        for prod in reaction.products
            (!haskey(reactants,prod.reactant)) && (reactants[prod.reactant] = r_count += 1)
        end
    end
    return reactants
end

#Generates a dictionary with all parameters.
function get_parameters(p)
    parameters = OrderedDict{Symbol,Int}()
    p_count = 0    ::Int
    for parameter in p
        (!haskey(parameters,parameter)) && (parameters[parameter] = p_count += 1)
    end
    return parameters
end

#Produces an array of expressions. Each entry corresponds to a line in the function f, which constitutes the deterministic part of the system. The Expressions can be used for debugging, making LaTex code, or creating the real f function for simulating the network.
function get_f(reactions::Vector{ReactionStruct}, reactants::OrderedDict{Symbol,Int})
    f = Vector{Expr}(undef,length(reactants))
    @inbounds for i = 1:length(f)
        f[i] = :(internal_var___du[$i] = $(Expr(:call, :+)))
    end
    @inbounds for reaction in reactions
        rl = recursive_clean!(reaction.rate_DE)        
        @inbounds for r in reaction.netstoich
            sidx = reactants[r.reactant]
            scoef = r.stoichiometry
            ex = recursive_clean!(:($scoef * $(deepcopy(reaction.rate_DE))))
            !(ex isa Number && iszero(ex)) && push!(f[sidx].args[2].args, ex)
        end
    end

    @inbounds for line in f
        if length(line.args[2].args) == 1
            @assert line.args[2].args[1] == :+
            line.args[2] = 0
        else
            line.args[2] = clean_subtractions(line.args[2])
        end
    end

    return f
end

#Produces an array of expressions. Each entry corresponds to a line in the function g, which constitutes the stochastic part of the system. Uses the Guillespie Approach for creating Langevin equations.  The Expressions can be used for debugging, making LaTex code, or creating the real f function for simulating the network.
function get_g(reactions::Vector{ReactionStruct}, reactants::OrderedDict{Symbol,Int}, scale_noise::Symbol)
    numspec = length(reactants)
    numrxs = length(reactions)
    g = Vector{Expr}(undef,numspec*numrxs)
    idx = 0
    @inbounds for k = 1:numrxs        
        # initialize to zero for all reactions
        @inbounds for (ssym,sidx) in reactants
            g[idx += 1] = :(internal_var___du[$sidx,$k] = 0)
        end

        ns = reactions[k].netstoich        
        @inbounds for r in ns
            sidx = reactants[r.reactant]
            scoef = r.stoichiometry
            ex = recursive_clean!( :($scoef * $scale_noise * sqrt(abs($(deepcopy(reactions[k].rate_DE))))) )
            g[sidx + (k-1)*numspec] = :(internal_var___du[$sidx,$k] = $ex)
        end
    end

    return g
end

#Makes the various jacobian elements required.
function get_jacs(f_rhs::Vector{ExprValues}, syms::Vector{Symbol}, reactions::Vector{ReactionStruct}, 
                                                                   reactants::OrderedDict{Symbol,Int}, 
                                                                   parameters::OrderedDict{Symbol,Int};
                                                                   build_jac = true,                                                                   
                                                                   build_paramjac = true,
                                                                   sparse_jac = false,
                                                                   zeroout_jac = false)
    # jacobian logic                                                               
    if build_jac
        if sparse_jac
            jac_prototype, jac = calculate_sparse_jac(reactions, reactants, parameters)
            symjac = nothing
        else
            symjac = calculate_symjac(reactions, reactants) 
            jac = calculate_jac(deepcopy(symjac), reactants, parameters, zeroout_jac) 
            jac_prototype = nothing
        end
    else
        jac_prototype = jac = symjac = nothing
    end
    
    paramjac = build_paramjac ? calculate_paramjac(deepcopy(f_rhs), reactants, parameters) : nothing

    if isnothing(symjac)
        return (symjac, jac, jac_prototype, paramjac) 
    else
        return (Expr(:quote, symjac), jac, jac_prototype, paramjac)
    end
end

# assumes subs is non-empty!
function massaction_ratelaw_deriv(spec::Symbol, scoef::Int, rate::ExprValues, subs::Vector{ReactantStruct})
    # only case where we should get a constant derivative
    ((length(subs) == 1) && (subs[1].stoichiometry == one(scoef))) && (return rate)

    # otherwise have non-constant derivative
    rl = Expr(:call, :*, rate)
    denom = one(scoef)
    for sub in subs
        name = sub.reactant
        stoich = sub.stoichiometry        
        if name == spec
            expon = stoich - one(stoich)
            denom *= factorial(expon)
            if expon == one(expon)
                push!(rl.args, name)
            elseif expon > one(expon)
                push!(rl.args, Expr(:call, :^, name, expon))
            end
        else
            denom *= factorial(stoich)
            push!(rl.args, stoich > one(stoich) ? Expr(:call, :^, name, stoich) : name)
        end
    end
    (denom != one(denom)) && (rl.args[2] = :($(rl.args[2])/$denom))

    rl
end

# given an expression representing a running sum, a ratelaw like expression, 
# and a stoichiometric coefficient, update the running sum with the coef
# times the rate law. specialized for mass action reactions
@inline function addstoich_to_exprsum(ex::ExprValues, rl::ExprValues, scoef)
    if (ex isa Number) && iszero(ex)
        if scoef < zero(scoef)
            newex = (scoef == -one(scoef)) ? :(-$rl) : :($scoef * $rl)
        else
            newex = (scoef == one(scoef)) ? rl : :($scoef * $rl)
        end
    else                        
        if scoef < zero(scoef)
            scoef = -scoef
            newex = (scoef == one(scoef)) ? rl : :($scoef * $rl)
            newex = :($ex - $newex)
        else
            newex = (scoef == one(scoef)) ? rl : :($scoef * $rl)
            newex = :($ex + $newex)
        end
    end    
    newex
end

# given an expression representing a running sum, a ratelaw like expression, 
# and a stoichiometric coefficient, update the running sum with the coef
# times the rate law. Specialized for SymEngine ratelaws.
@inline function addstoich_to_exprsum(ex::ExprValues, rl::SymEngine.Basic, scoef)
    if (ex isa Number) && iszero(ex)
        newex = Meta.parse(string(scoef*rl))
    else
        if scoef < zero(scoef)
            newex = :($ex - $(Meta.parse(string(-scoef*rl))))
        else
            newex = :($ex + $(Meta.parse(string(scoef*rl))))
        end
    end
    newex
end

# get the correct index type for a dict mapping (i,j) to ExprValues, 
# or for a Matrix of ExprValues
@inline getkeytype(d::Dict)                       = keytype(d)
@inline getkeytype(d::Matrix)                     = CartesianIndex{2}
@inline resolveindex(T::Type,i,j)                 = T(i,j)
@inline resolveindex(T::Type{Tuple{Int,Int}},i,j) = (i,j)
@inline checkforkey(d::Dict, key)                 = haskey(d, key)
@inline checkforkey(d::Matrix, key)               = true

# generate Jacobian. Using preceding functions it supports jacexprs as
# a Dict mapping (i,j) => ExprValues or as a dense matrix of ExprValues
function jac_as_exprvalues!(jacexprs, reactions, reactants) 
    internal_vars = [Symbol(:internal_variable___,var) for var in keys(reactants)]
    ivtosym = Dict(zip(internal_vars, keys(reactants)))
    symtoiv = Dict(zip(keys(reactants), internal_vars))

    @inbounds for rx in reactions         
        if rx.is_pure_mass_action 
            for sub in rx.substrates
                spec  = sub.reactant
                scoef = sub.stoichiometry
                @inbounds j = reactants[spec]

                # derivative with respect to this species
                dratelaw = massaction_ratelaw_deriv(spec, scoef, rx.rate_org, rx.substrates)

                # determine stoichiometric coefficent to multiply by 
                @inbounds for ns in rx.netstoich
                    key = resolveindex(getkeytype(jacexprs), reactants[ns.reactant], j)
                    if checkforkey(jacexprs, key) 
                        jacexprs[key] = addstoich_to_exprsum(jacexprs[key], deepcopy(dratelaw), ns.stoichiometry)
                    else
                        jacexprs[key] = addstoich_to_exprsum(0, deepcopy(dratelaw), ns.stoichiometry)
                    end
                end
            end            
        else
            ratelaw = SymEngine.Basic(recursive_replace!(recursive_clean!(deepcopy(rx.rate_DE)), symtoiv))
            @inbounds for dep in rx.dependants
                #dratelaw = diff(ratelaw, dep)
                dratelaw = diff(ratelaw, symtoiv[dep])
                j = reactants[dep]

                @inbounds for ns in rx.netstoich
                    key = resolveindex(getkeytype(jacexprs), reactants[ns.reactant], j)
                    if checkforkey(jacexprs, key) 
                        jacexprs[key] = addstoich_to_exprsum(jacexprs[key], dratelaw, ns.stoichiometry)
                    else
                        jacexprs[key] = addstoich_to_exprsum(0, dratelaw, ns.stoichiometry)
                    end
                    jacexprs[key] = recursive_replace!(jacexprs[key], ivtosym)
                end
            end
        end
    end
    nothing
end

# generate a dense matrix of ExprValues representing the Jacobian
function calculate_symjac(reactions, reactants)
    nspecs   = length(reactants) 
    jacexprs = Matrix{ExprValues}(undef, nspecs, nspecs)
    fill!(jacexprs, 0)
    jac_as_exprvalues!(jacexprs, reactions, reactants)
    jacexprs
end

# Generate the Jacobian evaluation function as a single expression.
# zeroout_jac kwarg controls whether to first fill the matrix with zeros and then only fill in non-zero entries
# vs looping through and putting expressions for entries that are always zero.
# This is needed for sufficiently large dense matrices to allow compilation. (e.g. 1000x1000)
function calculate_jac(symjac::Matrix{ExprValues}, reactants::OrderedDict{Symbol,Int}, parameters::OrderedDict{Symbol,Int}, zeroout_jac=false)
    func_body = Expr(:block)
    zeroout_jac && push!(func_body.args, :(fill!(internal___var___J, zero(eltype(internal___var___J)))))
    @inbounds for j = 1:size(symjac)[2], i = 1:size(symjac)[1]        
        if symjac[i,j] isa Number
            ex = symjac[i,j]
            (!zeroout_jac || !iszero(ex)) && push!(func_body.args, :(internal___var___J[$i,$j] = $(ex)))
        else
            ex = recursive_replace!(symjac[i,j],(reactants,:internal___var___u), (parameters, :internal___var___p))
            splitplus!(ex)
            push!(func_body.args, :(internal___var___J[$i,$j] = $(ex)))
        end
    end
    push!(func_body.args,:(return nothing))

    return :((internal___var___J,internal___var___u,internal___var___p,t) -> @inbounds $func_body)
end

# convert Dict mapping (i,j) => val to sparse matrix
function dict_to_sparsemat(d, m, n; vals=nothing)
    V   = isnothing(vals) ? collect(values(d)) : vals
    len = length(d)
    I   = Vector{Int}(undef,len)
    J   = Vector{Int}(undef,len)
    for (i,k) in enumerate(keys(d))
        I[i] = k[1]
        J[i] = k[2]
    end    
    return sparse(I,J,V,m,n)
end

# create a sparse Jacobian
function calculate_sparse_jac(reactions, reactants, parameters)

    # get the elements and structure of sparse matrix
    jacexprs = Dict{Tuple{Int,Int},ExprValues}()
    jac_as_exprvalues!(jacexprs, reactions, reactants)
    jac_prototype = dict_to_sparsemat(jacexprs, length(reactants), length(reactants), vals=ones(length(jacexprs)))

    # build the function
    jfun = Expr(:block)
    nspecs = length(reactants)
    push!(jfun.args, :(internal___var___Jvals = nonzeros(internal___var___J)))
    rows = rowvals(jac_prototype)
    @inbounds for j = 1:nspecs
        @inbounds for k in nzrange(jac_prototype,j)
            i = rows[k]
            ex = (jacexprs[(i,j)] isa Number) ? jacexprs[(i,j)] : recursive_replace!(jacexprs[(i,j)],(reactants,:internal___var___u), (parameters,:internal___var___p))
            push!(jfun.args, :(internal___var___Jvals[$k] = $ex))
        end
    end
    push!(jfun.args, :(return nothing))
    jfunex = :((internal___var___J,internal___var___u,internal___var___p,t) -> @inbounds $jfun)
    
    return jac_prototype, jfunex
end

#Makes the Jacobian, with respect to parameter values. 
function calculate_paramjac(f_rhs::Vector{ExprValues}, reactants::OrderedDict{Symbol,Int}, parameters::OrderedDict{Symbol,Int})
    func_body = Expr(:block)
    for j = 1:length(parameters), i = 1:length(reactants)    
        paramjac_entry = Meta.parse(string(diff(SymEngine.Basic(f_rhs[i]), parameters.keys[j])))
        push!(func_body.args, :(internal___var___pJ[$i,$j] = $(recursive_replace!(paramjac_entry,(reactants,:internal___var___u), (parameters, :internal___var___p)))))
    end
    push!(func_body.args,:(return internal___var___pJ))
    return :((internal___var___pJ,internal___var___u,internal___var___p,t) -> @inbounds $func_body) 
end


### Jump related functions ###

#Creates expressions for jump affects and rates. Also creates and array with MassAction, ConstantRate and VariableRate Jumps.
function get_jumps(reactions::Vector{ReactionStruct}, reactants::OrderedDict{Symbol,Int}, parameters::OrderedDict{Symbol,Int}; minimal_jumps=false)
    rates = Vector{ExprValues}()
    affects = Vector{Vector{Expr}}()
    jumps = Expr(:tuple)
    idx = 0
    for reaction in reactions
        hastime = recursive_contains(:t,reaction.rate_SSA)
        if !minimal_jumps || !reaction.is_pure_mass_action || hastime
            idx += 1
            push!(rates, recursive_clean!(deepcopy(reaction.rate_SSA)))
            push!(affects, Expr[:(integrator.u[$(reactants[r.reactant])] += $(r.stoichiometry)) for r in reaction.netstoich])
            syntax_rate = recursive_replace!(deepcopy(rates[idx]), (reactants, :internal_var___u), (parameters, :internal_var___p))
            hastime ? push!(jumps.args,Expr(:call,:VariableRateJump)) : push!(jumps.args,Expr(:call,:ConstantRateJump))
            push!(jumps.args[idx].args, :((internal_var___u,internal_var___p,t) -> @inbounds $syntax_rate))
            push!(jumps.args[idx].args, :(integrator -> @inbounds $(expr_arr_to_block(deepcopy(affects[idx])))))
        end
    end
    return (Tuple(rates),Tuple(affects),jumps)
end

# generates regular jumps representation
function get_regularjumps(reactions::Vector{ReactionStruct}, reactants::OrderedDict{Symbol,Int}, parameters::OrderedDict{Symbol,Int}; minimal_jumps=false)
    reg_rates = Expr(:block)
    reg_c = Expr(:block)
    @inbounds for (idx,reaction) in enumerate(reactions)
        rate = recursive_clean!(deepcopy(reaction.rate_SSA))
        syntax_rate = recursive_replace!(rate, (reactants,:internal_var___u), (parameters, :internal_var___p))
        push!(reg_rates.args, :(@inbounds internal_var___out[$idx]= $syntax_rate))
        foreach(r -> push!(reg_c.args,:(@inbounds internal_var___dc[$(reactants[r.reactant]),$idx]=$(r.stoichiometry))), reaction.netstoich)
    end
    reg_jumps = :(RegularJump((internal_var___out,internal_var___u,internal_var___p,t)->$reg_rates,(internal_var___dc,internal_var___u,internal_var___p,t,internal_var___mark)->$reg_c,zeros($(length(reactants)),$(length(reactions)));constant_c=true))
    return reg_jumps
end


### Pre Defined Functions that can be inserted into the reaction rates ###

#Hill function made avaiable
hill_name = Set{Symbol}([:hill, :Hill, :h, :H, :HILL])
function hill(expr::Expr)
    return :($(expr.args[3])*($(expr.args[2])^$(expr.args[5]))/($(expr.args[4])^$(expr.args[5])+$(expr.args[2])^$(expr.args[5])))
end
hillR_name = Set{Symbol}([:hill_repressor, :hillr, :hillR, :HillR, :hR, :hR, :Hr, :HR, :HILLR])
function hillR(expr::Expr)
    return :($(expr.args[3])*($(expr.args[4])^$(expr.args[5]))/($(expr.args[4])^$(expr.args[5])+$(expr.args[2])^$(expr.args[5])))
end

#Michaelis menten function made avaiable.
mm_name = Set{Symbol}([:MM, :mm, :Mm, :mM, :M, :m])
function mm(expr::Expr)
    return :($(expr.args[3])*$(expr.args[2])/($(expr.args[4])+$(expr.args[2])))
end
#Michaelis menten function made avaiable.
mmR_name = Set{Symbol}([:mm_repressor, :MMR, :mmr, :mmR, :MmR, :mMr, :MR, :mr, :Mr, :mR])
function mmR(expr::Expr)
    return :($(expr.args[3])*$(expr.args[4])/($(expr.args[4])+$(expr.args[2])))
end

#Allows the user to define new function and enable the @reaction_network macro to see them.
macro reaction_func(expr)
    name = expr.args[1].args[1]
    args = expr.args[1].args[2:end]
    maths = expr.args[2].args[2]

    funcdict[name]  = x -> replace_names(maths, args, x)
end

