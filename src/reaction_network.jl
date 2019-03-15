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
    (f_expr, f, f_rhs, symjac, jac, paramjac, f_symfuncs) = genode_exprs(reactions, reactants, parameters, syms)
    odefun = :(ODEFunction(f; jac=$jac, jac_prototype=nothing, paramjac=$paramjac, syms=$syms))

    # expressions for SDEs
    (g_expr, g, g_funcs, p_matrix) = gensde_exprs(reactions, reactants, parameters, scale_noise)
    sdefun = :(SDEFunction(f, g; jac=$jac, jac_prototype=nothing, paramjac=$paramjac, syms=$syms))

    # expressions for jumps
    (jump_rate_expr, jump_affect_expr, jumps, regular_jumps) = get_jumps(reactions, reactants, parameters)

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

# SDE expressions
function gensde_exprs(reactions, reactants, parameters, scale_noise)
    g_expr   = get_g(reactions, reactants, scale_noise)
    g        = make_func(g_expr, reactants, parameters)
    g_funcs  = ExprValues[element.args[2] for element in g_expr]
    p_matrix = zeros(length(reactants), length(reactions))

    (g_expr,g,g_funcs,p_matrix)
end

# ODE expressions
function genode_exprs(reactions, reactants, parameters, syms; build_jac=true,
                                                              build_symfuncs=true)
    f_expr                = get_f(reactions, reactants)
    f                     = make_func(f_expr, reactants, parameters)
    f_rhs                 = ExprValues[element.args[2] for element in f_expr]
    symjac, jac, paramjac = build_jac ? get_jacs(f_rhs, syms, reactants, parameters) : (nothing,nothing,nothing)
    f_symfuncs            = build_symfuncs ? hcat([SymEngine.Basic(f) for f in f_rhs]) : nothing

    (f_expr,f,f_rhs,symjac,jac,paramjac,f_symfuncs)
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
    rate_org::ExprValues
    rate_DE::ExprValues
    rate_SSA::ExprValues
    dependants::Vector{Symbol}
    is_pure_mass_action::Bool

    function ReactionStruct(s::Vector{ReactantStruct}, p::Vector{ReactantStruct},
                            ro::ExprValues, rde::ExprValues, rssa::ExprValues,
                            dep::Vector{Symbol}, isma::Bool)
        new(s,p,ro,rde,rssa,dep,isma)
    end

    function ReactionStruct(sub_line::ExprValues, prod_line::ExprValues, rate::ExprValues, use_mass_kin::Bool)
        sub = add_reactants!(sub_line,1,Vector{ReactantStruct}(undef,0))
        prod = add_reactants!(prod_line,1,Vector{ReactantStruct}(undef,0))

        rate_DE = mass_rate_DE(sub, use_mass_kin, rate)
        rate_SSA =  mass_rate_SSA(sub, use_mass_kin, rate)
        new(sub, prod, rate, rate_DE, rate_SSA, [], use_mass_kin)
    end
    function ReactionStruct(r::ReactionStruct, syms::Vector{Symbol})
        deps = unique!(recursive_content(r.rate_DE,syms,Vector{Symbol}()))
        is_ma = r.is_pure_mass_action && (length(recursive_content(r.rate_org,syms,Vector{Symbol}()))==0)
        new(r.substrates, r.products, r.rate_org, r.rate_DE, r.rate_SSA, deps, is_ma)
    end
end

#Calculates the rate used by ODEs and SDEs. If we want to use masskinetics we have to include substrate concentration, taking higher order terms into account.
function mass_rate_DE(substrates::Vector{ReactantStruct}, use_mass_kin::Bool, old_rate::ExprValues)
    rate = Expr(:call, :*, old_rate)
    use_mass_kin && foreach(sub -> push!(rate.args,:($(Expr(:call, :^, sub.reactant, sub.stoichiometry))/$(factorial(sub.stoichiometry)))), substrates)
    return rate
end

#Calculates the rate used by SSAs. If we want to use masskinetics we have to include substrate concentration, taking higher order terms into account.
function mass_rate_SSA(substrates::Vector{ReactantStruct}, use_mass_kin::Bool, old_rate::ExprValues)
    rate = Expr(:call, :*, old_rate)
    use_mass_kin && foreach(sub -> push!(rate.args, :(binomial($(sub.reactant),$(sub.stoichiometry)))), substrates)
    return rate
end

#Returns the length of a expression tuple, or 1 if it is not an expression tuple (probably a  Symbol/Numerical).
function tup_leng(ex::ExprValues)
    (typeof(ex)==Expr && ex.head == :tuple) && (return length(ex.args))
    return 1
end

#Gets the i'th element in a expression tuple, or return the input itself if it is not an expression tuple (probably a  Symbol/Numerical).
function get_tup_arg(ex::ExprValues,i::Int)
    (tup_leng(ex) == 1) && (return ex)
    return ex.args[i]
end

#Takes a reaction line and creates reactions from it and pushes those to the reaction array. Used to creat multiple reactions from e.g. 1.0, (X,Y) --> 0.
function push_reactions(reactions::Vector{ReactionStruct}, sub_line::ExprValues, prod_line::ExprValues, rate::ExprValues, use_mass_kin::Bool)
    lengs = [tup_leng(sub_line), tup_leng(prod_line), tup_leng(rate)]
    (count(lengs.==1) + count(lengs.==maximum(lengs)) < 3) && (throw("malformed reaction"))
    for i = 1:maximum(lengs)
        push!(reactions, ReactionStruct(get_tup_arg(sub_line,i), get_tup_arg(prod_line,i), get_tup_arg(rate,i), use_mass_kin))
    end
end

#Recursive function that loops through the reactants in an reaction line and finds the reactants and their stochiometry. Recursion makes it able to handle e.g. 2(X+Y+3(Z+XY)) (probably one will not need it though).
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
    for i = 1:length(f)
        f[i] = :(internal_var___du[$i] = $(Expr(:call, :+)))
    end
    for reaction in deepcopy(reactions)
        for reactant in union(getfield.(reaction.products, :reactant),getfield.(reaction.substrates, :reactant))
            ex = recursive_clean!(:($(get_stoch_diff(reaction,reactant)) * $(reaction.rate_DE)))
            !(ex isa Number && iszero(ex)) && push!(f[reactants[reactant]].args[2].args, ex)
        end
    end

    for line in f
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
    g = Vector{Expr}(undef,length(reactions)*length(reactants))
    idx = 0
    for reactant in keys(reactants), i = 1:length(reactions)
        g[idx += 1] = recursive_clean!(:(internal_var___du[$(reactants[reactant]),$i] = $scale_noise * $(get_stoch_diff(reactions[i],reactant)) * sqrt(abs($(deepcopy(reactions[i].rate_DE))))))
    end
    return g
end

#Computes how much the stoichiometry in a single reactant changes for a reaction. Only really interesting if the reactant is both a product and substrate.
function get_stoch_diff(reaction::ReactionStruct, reactant::Symbol)
    stoch = 0
    for prod in reaction.products
        (reactant == prod.reactant) && (stoch += prod.stoichiometry)
    end
    for sub in reaction.substrates
        (reactant == sub.reactant) && (stoch -= sub.stoichiometry)
    end
    return stoch
end

function splitplus!(ex)
  dosplit = ex.head == :(=) && ex.args[2] isa Expr && ex.args[2].head == :call && ex.args[2].args[1] == :(+)
  if dosplit
    summands = ex.args[2].args[2:end]
    ex.args[2] = foldl((x,y)->(:(($x + $y))), summands)
  end
  dosplit
end

#Creates an expression which can be evaluated to an actual function. Input is an array of expression were each entry is a line in the function. Uses the array of expressions generated in either get_f or get_g.
function make_func(func_expr::Vector{Expr},reactants::OrderedDict{Symbol,Int}, parameters::OrderedDict{Symbol,Int})
    system = Expr(:block)
    for func_line in deepcopy(func_expr)
        ex = recursive_replace!(func_line, (reactants,:internal_var___u), (parameters, :internal_var___p))
        splitplus!(ex)
        push!(system.args,ex)
    end
    push!(system.args, :(nothing))
    return :((internal_var___du,internal_var___u,internal_var___p,t) -> @inbounds $system)
end

#Creates expressions for jump affects and rates. Also creates and array with MassAction, ConstantRate and VariableRate Jumps.
function get_jumps(reactions::Vector{ReactionStruct}, reactants::OrderedDict{Symbol,Int}, parameters::OrderedDict{Symbol,Int}; minimal_jumps=false)
    rates = Vector{ExprValues}(undef,length(reactions))
    affects = Vector{Vector{Expr}}(undef,length(reactions))
    jumps = Expr(:tuple)
    reg_rates = Expr(:block)
    reg_c = Expr(:block)
    idx = 0
    for reaction in deepcopy(reactions)
        rates[idx += 1] = recursive_clean!(reaction.rate_SSA)
        affects[idx] = Vector{Expr}(undef,0)
        reactant_set = union(getfield.(reaction.products, :reactant),getfield.(reaction.substrates, :reactant))
        foreach(r -> push!(affects[idx], :(integrator.u[$(reactants[r])] += $(get_stoch_diff(reaction,r)))), reactant_set)
        syntax_rate = recursive_replace!(deepcopy(rates[idx]), (reactants,:internal_var___u), (parameters, :internal_var___p))

        if minimal_jumps && reaction.is_pure_mass_action
            recursive_contains(:t,rates[idx]) && push!(jumps.args,Expr(:call,:VariableRateJump))
        #    ma_sub_stoch = :(reactant_stoich = [[]])
        #    ma_stoch_change = :(reactant_stoich = [[]])
        #    foreach(sub -> push!(ma_sub_stoch.args[2].args[1].args),:($(reactants[sub.reactant])=>$(sub.stoichiometry)),reaction.substrates)
        #    foreach(reactant -> push!(ma_stoch_change.args[2].args[1].args),:($(reactants[reactant.reactant])=>$(get_stoch_diff(reaction,reactant))),reaction.substrates)
        #    push!(jumps.args,:(MassActionJump($(reaction.rate_org),$(ma_sub_stoch),$(ma_stoch_change))))
        else
            recursive_contains(:t,rates[idx]) ? push!(jumps.args,Expr(:call,:VariableRateJump)) : push!(jumps.args,Expr(:call,:ConstantRateJump))
            push!(jumps.args[idx].args, :((internal_var___u,internal_var___p,t) -> @inbounds $syntax_rate))
            push!(jumps.args[idx].args, :(integrator -> @inbounds $(expr_arr_to_block(deepcopy(affects[idx])))))
        end
        push!(reg_rates.args,:(@inbounds internal_var___out[$idx]= $syntax_rate))
        foreach(r -> push!(reg_c.args,:(@inbounds internal_var___dc[$(reactants[r]),$idx]=$(get_stoch_diff(reaction,r)))), reactant_set)
    end
    reg_jumps = :(RegularJump((internal_var___out,internal_var___u,internal_var___p,t)->$reg_rates,(internal_var___dc,internal_var___u,internal_var___p,t,internal_var___mark)->$reg_c,zeros($(length(reactants)),$(length(reactions)));constant_c=true))
    return (Tuple(rates),Tuple(affects),jumps,reg_jumps)
end

"""
clean_subtractions(ex::Expr)
Replace additions of negative terms with subtractions.
This is a fairly stupid function which is designed for a specific problem
with reaction networks. It is neither recursive nor very general.
Return :: cleaned out expression

From Latexify.jl with permission:
[see](https://github.com/JuliaDiffEq/DiffEqBiological.jl/issues/89#issuecomment-462147882)
"""
function clean_subtractions(ex::Expr)
    ex.args[1] != :+ && return ex

    term = ex.args[2]

    ### Sort out the first term
    if term isa Expr && length(term.args) >= 3 && term.args[1:2] == [:*, -1]
        result = :(- *($(term.args[3:end]...)))
    else
        result = :($term)
    end

    ### Sort out the other terms
    for term in ex.args[3:end]
        if term isa Expr && length(term.args) >= 3 && term.args[1:2] == [:*, -1]
            result = :($result - *($(term.args[3:end]...)))
        else
            result = :($result + $term)
        end
    end
    return result
end

#Recursively traverses an expression and removes things like X^1, 1*X. Will not actually have any affect on the expression when used as a function, but will make it much easier to look at it for debugging, as well as if it is transformed to LaTeX code.
function recursive_clean!(expr::ExprValues)
    (expr == :no___noise___scaling) && (return 1)
    (typeof(expr)!=Expr) && (return expr)
    for i = 1:length(expr.args)
        expr.args[i] = recursive_clean!(expr.args[i])
    end
    (expr.args[1] == :^) && (expr.args[3] == 1) && (return expr.args[2])
    if expr.args[1] == :*
        in(0,expr.args) && (return 0)
        i = 1
        while (i = i + 1) <= length(expr.args)
             if (typeof(expr.args[i]) == Expr) && (expr.args[i].head == :call) && (expr.args[i].args[1] == :*)
                 for arg in expr.args[i].args
                     (arg != :*) && push!(expr.args, arg)
                 end
             end
        end
        for i = length(expr.args):-1:2
            (typeof(expr.args[i]) == Expr) && (expr.args[i].head == :call) && (expr.args[i].args[1] == :*) && deleteat!(expr.args,i)
            (expr.args[i] == 1) && deleteat!(expr.args,i)
        end
        (length(expr.args) == 2) && (return expr.args[2])                   # We have a multiplication of only one thing, return only that thing.
        (length(expr.args) == 1) && (return 1)                              #We have only * and no real argumenys.
        (length(expr.args) == 3) && (expr.args[2] == -1) && return :(-$(expr.args[3]))
        (length(expr.args) == 3) && (expr.args[3] == -1) && return :(-$(expr.args[2]))
    end
    if expr.head == :call
        (expr.args[1] == :/) && (expr.args[3] == 1) && (return expr.args[2])
        haskey(funcdict, expr.args[1]) && return funcdict[expr.args[1]](expr.args[2:end])
        in(expr.args[1],hill_name) && return hill(expr)
        in(expr.args[1],hillR_name) && return hillR(expr)
        in(expr.args[1],mm_name) && return mm(expr)
        in(expr.args[1],mmR_name) && return mmR(expr)
        (expr.args[1] == :binomial) && (expr.args[3] == 1) && return expr.args[2]
        #@isdefined($(expr.args[1])) || error("Function $(expr.args[1]) not defined.")
    end
    return expr
end

#Recursively traverses an expression and replace instances of variables and parmaters with things that the DifferentialEquations packakes simulation algorithms can understand. E.g. X --> u[1], kB1 --> p[1] etc.
function recursive_replace!(expr::ExprValues, replace_requests::Tuple{OrderedDict{Symbol,Int},Symbol}...)
    if typeof(expr) == Symbol
        for rr in replace_requests
            (haskey(rr[1],expr)) && (return :($(rr[2])[$(rr[1][expr])]))
        end
    elseif typeof(expr) == Expr
        for i = 1:length(expr.args)
            expr.args[i] = recursive_replace!(expr.args[i], replace_requests...)
        end
    end
    return expr
end

#Recursively traverses an expression and replaces a symbol with another.
function recursive_replace!(expr::ExprValues, replace_requests::Dict{Symbol,Symbol})
    if typeof(expr) == Symbol
        haskey(replace_requests,expr) && return replace_requests[expr]
    elseif typeof(expr) == Expr
        for i = 1:length(expr.args)
            expr.args[i] = recursive_replace!(expr.args[i], replace_requests)
        end
    end
    return expr
end

#Recursive Contains, checks whenever an expression contains a certain symbol.
function recursive_contains(s,ex)
    (typeof(ex)!=Expr) && (return s==ex)
    for arg in ex.args
        recursive_contains(s,arg) && (return true)
    end
    return false
end

#Parses an expression, and returns a set with all symbols in the expression, which is also a part of the provided vector with symbols (syms).
function recursive_content(ex,syms::Vector{Symbol},content::Vector{Symbol})
    if ex isa Symbol
        in(ex,syms) && push!(content,ex)
    elseif ex isa Expr
        foreach(arg -> recursive_content(arg,syms,content), ex.args)
    end
    return content
end

function recursive_content(ex,symsmap::OrderedDict{Symbol,Int},content::Vector{Symbol})
    if ex isa Symbol        
        haskey(symsmap,ex) && push!(content,ex)
    elseif ex isa Expr  
        foreach(arg -> recursive_content(arg,symsmap,content), ex.args)
    end
    return content
end

#Makes the various jacobian elements required.
function get_jacs(f_rhs::Vector{ExprValues}, syms::Vector{Symbol}, reactants::OrderedDict{Symbol,Int}, parameters::OrderedDict{Symbol,Int})
    symjac = calculate_symjac(deepcopy(f_rhs), syms)
    jac = calculate_jac(deepcopy(symjac), reactants, parameters)
    paramjac = calculate_paramjac(deepcopy(f_rhs), reactants, parameters)
    return (Expr(:quote, symjac), jac, paramjac)
end

#Makes the Symbolic Jacobian.
function calculate_symjac(f_rhs::Vector{ExprValues}, syms)
    n = length(syms); internal_vars = [Symbol(:internal_variable___,var) for var in syms]
    symfuncs = [SymEngine.Basic(recursive_replace!(f,Dict(zip(syms,internal_vars)))) for f in f_rhs]
    jacexprs = Matrix{ExprValues}(undef, n, n)
    for j = 1:n, i = 1:n
        symjacij = diff(symfuncs[i],internal_vars[j])
        jacexprs[i,j] = :($(recursive_replace!(Meta.parse(string(symjacij)),Dict(zip(internal_vars,syms)))))
    end
    jacexprs
end

#Makes the Jacobian.
function calculate_jac(symjac::Matrix{ExprValues}, reactants::OrderedDict{Symbol,Int}, parameters::OrderedDict{Symbol,Int})
    func_body = Expr(:block)
    for j = 1:size(symjac)[2], i = 1:size(symjac)[1]
        push!(func_body.args, :(internal___var___J[$i,$j] = $(recursive_replace!(symjac[i,j],(reactants,:internal___var___u), (parameters, :internal___var___p)))))
    end
    push!(func_body.args,:(return internal___var___J))
    return :((internal___var___J,internal___var___u,internal___var___p,t) -> @inbounds $func_body)
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

#Turns an array of expressions to a expression block with corresponding expressions.
function expr_arr_to_block(exprs)
  block = :(begin end)
  foreach(expr -> push!(block.args, expr), exprs)
  return block
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

#This will be called whenever a function stored in funcdict is called.
function replace_names(expr, old_names, new_names)
    mapping = Dict(zip(old_names, new_names))
    MacroTools.postwalk( x -> x in old_names ? x= mapping[x] : x, expr)
end
