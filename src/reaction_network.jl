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
macro reaction_network(name, ex::Expr, p...)
    coordinate(name, MacroTools.striplines(ex), p, :no___noise___scaling)
end

#Macro to create a reaction network model. Multiple dispatch is used to allow for SDE noise scalling.
macro reaction_network(name, scale_noise, ex::Expr, p...)
    in(scale_noise, p) || (p = (p..., scale_noise))
    coordinate(name, MacroTools.striplines(ex), p, scale_noise)
end

#Used to give a warning if someone uses the old macro.
macro reaction_network(ex::Expr)
    error("The Reaction Network DSL have been deprecated in favor for a new one. With only slight modification old code can be made to work with the new DSL. In addition the new one provides lots of additional functionality. Please view the documentation for more information.")
end

#Declare various arrow types symbols used for the empty set (also 0).
empty_set = Set{Symbol}([:∅])
fwd_arrows = Set{Symbol}([:>, :→, :↣, :↦, :⇾, :⟶, :⟼, :⥟, :⥟, :⇀, :⇁, :⇒, :⟾])
bwd_arrows = Set{Symbol}([:<, :←, :↢, :↤, :⇽, :⟵, :⟻, :⥚, :⥞, :↼, :↽, :⇐, :⟽])
double_arrows = Set{Symbol}([:↔, :⟷, :⇄, :⇆, :⇔, :⟺])
no_mass_arrows = Set{Symbol}([:⇐, :⟽, :⇒, :⟾, :⇔, :⟺])      #Using this arrows will disable the program from multiplying reaction rates with the substrate concentrations. Gives user full control of reaction rates.

funcdict = Dict{Symbol, Function}()                             #Stores user defined functions.

#Coordination function, actually does all the work of the macro.
function coordinate(name, ex::Expr, p, scale_noise)
    reactions = get_reactions(ex)           ::Vector{ReactionStruct}
    reactants = get_reactants(reactions)    ::OrderedDict{Symbol,Int}
    parameters = get_parameters(p)          ::OrderedDict{Symbol,Int}

    syms = collect(keys(reactants))
    params = collect(keys(parameters))
    (in(:t,union(syms,params))) && error("t is reserved for the time variable and may neither be used as a reactant nor a parameter")

    update_reaction_info(reactions,syms)

    f_expr = get_f(reactions, reactants)
    f = make_func(f_expr, reactants, parameters)

    g_expr = get_g(reactions, reactants, scale_noise)
    g = make_func(g_expr, reactants, parameters)
    p_matrix = zeros(length(reactants), length(reactions))

    (jump_rate_expr, jump_affect_expr, jumps, regular_jumps) = get_jumps(reactions, reactants, parameters)

    f_rhs = [element.args[2] for element in f_expr]
    #symjac = Expr(:quote, calculate_jac(deepcopy(f_rhs), syms))
    f_symfuncs = hcat([SymEngine.Basic(f) for f in f_rhs])

    # Build the type
    exprs = Vector{Expr}(undef,0)

    ## only get the right-hand-side of the equations.
    f_funcs = [element.args[2] for element in f_expr]
    g_funcs = [element.args[2] for element in g_expr]

    typeex,constructorex = maketype(name, f, f_funcs, f_symfuncs, g, g_funcs, jumps, regular_jumps, Meta.quot(jump_rate_expr), Meta.quot(jump_affect_expr), p_matrix, syms; params=params, reactions=reactions)#, symjac=symjac)

    push!(exprs,typeex)
    push!(exprs,constructorex)

    ## Overload the type so that it can act as a function.
    overloadex = :(((f::$name))(du, u, p, t::Number) = f.f(du, u, p, t)) |> esc
    push!(exprs,overloadex)

    ## Add a method which allocates the `du` and returns it instead of being inplace
    overloadex = :(((f::$name))(u,p,t::Number) = (du=similar(u); f(du,u,p,t); du)) |> esc
    push!(exprs,overloadex)

    # export type constructor
    def_const_ex = :(($name)()) |> esc
    push!(exprs,def_const_ex)

    expr_arr_to_block(exprs)
end

#Generates a vector containing a number of reaction structures, each containing the infromation about one reaction.
function get_reactions(ex::Expr)
    reactions = Vector{ReactionStruct}(undef,0)
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
    rate_org::Any
    rate_DE::Any
    rate_SSA::Any
    dependants::Vector{Symbol}
    is_pure_mass_action::Bool

    function ReactionStruct(sub_line::Any, prod_line::Any, rate::Any, use_mass_kin::Bool)
        sub = add_reactants!(sub_line,1,Vector{ReactantStruct}(undef,0))
        prod = add_reactants!(prod_line,1,Vector{ReactantStruct}(undef,0))

        rate_DE = mass_rate_DE(sub, use_mass_kin, rate)
        rate_SSA =  mass_rate_SSA(sub, use_mass_kin, rate)
        new(sub, prod, rate, rate_DE, rate_SSA, [], use_mass_kin)
    end
    function ReactionStruct(r::ReactionStruct, syms::Vector{Symbol})
        deps = recursive_content(r.rate_DE,syms,Vector{Symbol}())
        is_ma = r.is_pure_mass_action && (length(recursive_content(r.rate_org,syms,Vector{Symbol}()))==0)
        new(r.substrates, r.products, r.rate_org, r.rate_DE, r.rate_SSA, deps, is_ma)
    end
end

#Calculates the rate used by ODEs and SDEs. If we want to use masskinetics we have to include substrate concentration, taking higher order terms into account.
function mass_rate_DE(substrates::Vector{ReactantStruct}, use_mass_kin::Bool, old_rate::Any)
    rate = Expr(:call, :*, old_rate)
    use_mass_kin && foreach(sub -> push!(rate.args,:($(Expr(:call, :^, sub.reactant, sub.stoichiometry))/$(factorial(sub.stoichiometry)))), substrates)
    return rate
end

#Calculates the rate used by SSAs. If we want to use masskinetics we have to include substrate concentration, taking higher order terms into account.
function mass_rate_SSA(substrates::Vector{ReactantStruct}, use_mass_kin::Bool, old_rate::Any)
    rate = Expr(:call, :*, old_rate)
    use_mass_kin && foreach(sub -> push!(rate.args, :(binomial($(sub.reactant),$(sub.stoichiometry)))), substrates)
    return rate
end

#Returns the length of a expression tuple, or 1 if it is not an expression tuple (probably a  Symbol/Numerical).
function tup_leng(ex::Any)
    (typeof(ex)==Expr && ex.head == :tuple) && (return length(ex.args))
    return 1
end

#Gets the i'th element in a expression tuple, or return the input itself if it is not an expression tuple (probably a  Symbol/Numerical).
function get_tup_arg(ex::Any,i::Int)
    (tup_leng(ex) == 1) && (return ex)
    return ex.args[i]
end

#Takes a reaction line and creates reactions from it and pushes those to the reaction array. Used to creat multiple reactions from e.g. 1.0, (X,Y) --> 0.
function push_reactions(reactions::Vector{ReactionStruct}, sub_line::Any, prod_line::Any, rate::Any, use_mass_kin::Bool)
    lengs = [tup_leng(sub_line), tup_leng(prod_line), tup_leng(rate)]
    (count(lengs.==1) + count(lengs.==maximum(lengs)) < 3) && (throw("malformed reaction"))
    for i = 1:maximum(lengs)
        push!(reactions, ReactionStruct(get_tup_arg(sub_line,i), get_tup_arg(prod_line,i), get_tup_arg(rate,i), use_mass_kin))
    end
end

#Recursive function that loops through the reactants in an reaction line and finds the reactants and their stochiometry. Recursion makes it able to handle e.g. 2(X+Y+3(Z+XY)) (probably one will not need it though).
function add_reactants!(ex::Any, mult::Int, reactants::Vector{ReactantStruct})
    if typeof(ex)!=Expr
        (ex == 0 || in(ex,empty_set)) && (return reactants)
        if in(ex, getfield.(reactants,:reactant))
            idx = find(x -> x==ex ,getfield.(reactants,:reactant))[1]
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

#For each reaction, sets its dependencies and whenever it is a pure mass action reaction.
function update_reaction_info(reactions::Vector{ReactionStruct},syms::Vector{Symbol})
    for i = 1:length(reactions)
        reactions[i] = ReactionStruct(reactions[i],syms)
    end
end

#Produces an array of expressions. Each entry corresponds to a line in the function f, which constitutes the deterministic part of the system. The Expressions can be used for debugging, making LaTex code, or creating the real f function for simulating the network.
function get_f(reactions::Vector{ReactionStruct}, reactants::OrderedDict{Symbol,Int})
    f = Vector{Expr}(undef,length(reactants))
    for i = 1:length(f)
        f[i] = :(internal_var___du[$i] = $(Expr(:call, :+)))
    end
    for reaction in deepcopy(reactions)
        for reactant in union(getfield.(reaction.products, :reactant),getfield.(reaction.substrates, :reactant))
            push!(f[reactants[reactant]].args[2].args, recursive_clean!(:($(get_stoch_diff(reaction,reactant)) * $(reaction.rate_DE))))
        end
    end
    return f
end

#Produces an array of expressions. Each entry corresponds to a line in the function g, which constitutes the stochastic part of the system. Uses the Guillespie Approach for creating Langevin equations.  The Expressions can be used for debugging, making LaTex code, or creating the real f function for simulating the network.
function get_g(reactions::Vector{ReactionStruct}, reactants::OrderedDict{Symbol,Int}, scale_noise::Symbol)
    g = Vector{Expr}(undef,length(reactions)*length(reactants))
    idx = 0
    for reactant in keys(reactants), i = 1:length(reactions)
            g[idx += 1] = recursive_clean!(:(internal_var___du[$(reactants[reactant]),$i] = $scale_noise * $(get_stoch_diff(reactions[i],reactant)) * sqrt($(deepcopy(reactions[i].rate_DE)))))
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

#Creates an expression which can be evaluated to an actual function. Input is an array of expression were each entry is a line in the function. Uses the array of expressions generated in either get_f or get_g.
function make_func(func_expr::Vector{Expr},reactants::OrderedDict{Symbol,Int}, parameters::OrderedDict{Symbol,Int})
    system = Expr(:block)
    for func_line in deepcopy(func_expr)
        push!(system.args, recursive_replace!(func_line, (reactants,:internal_var___u), (parameters, :internal_var___p)))
    end
    return :((internal_var___du,internal_var___u,internal_var___p,t) -> $system)
end

#Creates expressions for jump affects and rates. Also creates and array with MassAction, ConstantRate and VariableRate Jumps.
function get_jumps(reactions::Vector{ReactionStruct}, reactants::OrderedDict{Symbol,Int}, parameters::OrderedDict{Symbol,Int})
    rates = Vector{Any}(undef,length(reactions))
    affects = Vector{Vector{Expr}}(undef,length(reactions))
    jumps = Expr(:tuple)
    reg_rates = Expr(:block)
    reg_c = Expr(:block)
    idx = 0
    for reaction in deepcopy(reactions)
        rates[idx += 1] = recursive_clean!(reaction.rate_SSA)
        affects[idx] = Vector{Expr}(undef,0)
        reactant_set = union(getfield.(reaction.products, :reactant),getfield.(reaction.substrates, :reactant))
        foreach(r -> push!(affects[idx],:(@inbounds integrator.u[$(reactants[r])] += $(get_stoch_diff(reaction,r)))), reactant_set)
        syntax_rate = recursive_replace!(deepcopy(rates[idx]), (reactants,:internal_var___u), (parameters, :internal_var___p))
        #if reaction.is_pure_mass_action
        #    ma_sub_stoch = :(reactant_stoich = [[]])
        #    ma_stoch_change = :(reactant_stoich = [[]])
        #    foreach(sub -> push!(ma_sub_stoch.args[2].args[1].args),:($(reactants[sub.reactant])=>$(sub.stoichiometry)),reaction.substrates)
        #    foreach(reactant -> push!(ma_stoch_change.args[2].args[1].args),:($(reactants[reactant.reactant])=>$(get_stoch_diff(reaction,reactant))),reaction.substrates)
        #    push!(jumps.args,:(MassActionJump($(reaction.rate_org),$(ma_sub_stoch),$(ma_stoch_change))))
        #else
            recursive_contains(:t,rates[idx]) ? push!(jumps.args,Expr(:call,:VariableRateJump)) : push!(jumps.args,Expr(:call,:ConstantRateJump))
            push!(jumps.args[idx].args, :((internal_var___u,internal_var___p,t) -> $syntax_rate))
            push!(jumps.args[idx].args, :(integrator -> $(expr_arr_to_block(deepcopy(affects[idx])))))
        #end
        push!(reg_rates.args,:(internal_var___out[$idx]=$syntax_rate))
        foreach(r -> push!(reg_c.args,:(internal_var___dc[$(reactants[r]),$idx]=$(get_stoch_diff(reaction,r)))), reactant_set)
    end
    reg_jumps = :(RegularJump((internal_var___out,internal_var___u,internal_var___p,t)->$reg_rates,(internal_var___dc,internal_var___u,internal_var___p,t,internal_var___mark)->$reg_c,zeros($(length(reactants)),$(length(reactions)));constant_c=true))
    return (Tuple(rates),Tuple(affects),jumps,reg_jumps)
end

#Recursively traverses an expression and removes things like X^1, 1*X. Will not actually have any affect on the expression when used as a function, but will make it much easier to look at it for debugging, as well as if it is transformed to LaTeX code.
function recursive_clean!(expr::Any)
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
        in(expr.args[1],mm_name) && return mm(expr)
        (expr.args[1] == :binomial) && (expr.args[3] == 1) && return expr.args[2]
        #@isdefined($(expr.args[1])) || error("Function $(expr.args[1]) not defined.")
    end
    return expr
end

#Recursively traverses an expression and replace instances of variables and parmaters with things that the DifferentialEquations packakes simulation algorithms can understand. E.g. X --> u[1], kB1 --> p[1] etc.
function recursive_replace!(expr::Any, replace_requests::Tuple{OrderedDict{Symbol,Int},Symbol}...)
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
function recursive_replace!(expr::Any, replace_requests::Dict{Symbol,Symbol})
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
    if typeof(ex)!=Expr
        in(ex,syms) && push!(content,ex)
    else
        foreach(arg -> recursive_content(arg,syms,content), ex.args)
    end
    return content
end

#Makes the Jacobian.
function calculate_jac(f_expr::Vector{Expr}, syms)
    n = length(syms); internal_vars = [Symbol(:internal_variable___,var) for var in syms]
    symjac = Matrix{SymEngine.Basic}(undef, n, n);
    symfuncs = [SymEngine.Basic(recursive_replace!(f,Dict(zip(syms,internal_vars)))) for f in f_expr]
    for i = 1:n, j = 1:n
        symjac[i,j] = diff(symfuncs[i],internal_vars[j])
    end
    @show symjac
    map!(symentry -> SymEngine.Basic(recursive_replace!(parse(string(symentry)),Dict(zip(internal_vars,syms)))),symjac)
    return symjac
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

#Michaelis menten function made avaiable.
mm_name = Set{Symbol}([:MM, :mm, :Mm, :mM, :M, :m])
function mm(expr::Expr)
    return :($(expr.args[3])*$(expr.args[2])/($(expr.args[4])+$(expr.args[2])))
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
