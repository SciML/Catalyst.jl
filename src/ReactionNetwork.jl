"""
Macro that inputs an expression corresponding to a reaction network and output a Reaction Network Structure that can be used as input to generation of SDE and ODE and Jump problems.
Most arrows accepted (both right, left and bi drectional arrows).
Note that while --> is a correct arrow, neither <-- nor <--> works.
Using non-filled arrows (⇐, ⟽, ⇒, ⟾, ⇔, ⟺) will disable mass kinetics and lets you cutomize reaction rates yourself.
Use 0 or ∅ for degradation/creation to/from nothing.
Example systems:
    rn = @reaction_network_new rType begin  #Creates a reaction network of type rType.
        2.0, X + Y --> XY                  #This will have reaction rate corresponding to 2.0*[X][Y]
        2.0, XY ← X + Y                    #Identical to 2.0, X + Y --> XY
        (2.0,1.0), XY ← X + Y              #Identical to reactions (2.0, X + Y --> XY) and (1.0, XY --> X + Y).
        2.0, X + Y ⟾ XY                   #Ignores mass kinetics. This will have reaction rate corresponding to 2.0*[X].
        hill(XY,2,2,2), X + Y --> XY       #Reaction inis activated by XY according to a hill function. hill(x,v,K,N).
        mm(XY,2,2), X + Y --> XY           #Reaction inis activated by XY according to a michaelis menten function. mm(x,v,K).
        2.0, (X,Y) --> 0                   #This corresponds to both X and Y degrading at rate 2.0.
        (2.0, 1.0), (X,Y) --> 0            #This corresponds to X and Y degrading at rates 2.0 and 1.0, respectively.
        2.0, (X1,Y1) --> (X2,Y2)           #X1 and Y1 becomes X2 and Y2, respectively, at rate 2.0.
    end

    kB = 2.0; kD = 1.0
    p = [kB, kD]
    p = []
    rn = @reaction_network_new type begin
        (kB, kD), X + Y --> XY            #Lets you define parameters outside on network. Parameters can be changed without recalling the network.
    end kB, kD
"""
#Macro to create a reaction network model.
macro reaction_network_new(name, ex::Expr, p...)
    coordinate(name, ex, p)
end

#Declare various arrow types symbols used for the empty set (also 0).
empty_set = Set{Symbol}([:∅])
fwd_arrows = Set{Symbol}([:>, :→, :↣, :↦, :⇾, :⟶, :⟼, :⥟, :⥟, :⇀, :⇁, :⇒, :⟾])
bwd_arrows = Set{Symbol}([:<, :←, :↢, :↤, :⇽, :⟵, :⟻, :⥚, :⥞, :↼, :↽, :⇐, :⟽])
double_arrows = Set{Symbol}([:↔, :⟷, :⥎, :⥐, :⇄, :⇆, :⇋, :⇌, :⇔, :⟺])
no_mass_arrows = Set{Symbol}([:⇐, :⟽, :⇒, :⟾, :⇔, :⟺])      #Using this arrows will disable the program from multiplying reaction rates with the substrate concentrations. Gives user full control of reaction rates.
disallowed_reactants = Set{Symbol}([:du, :u, :p, :t])           #These are not allowed since they are used in "return :((du,u,p,t) -> $system)", if a variable these gets replaced with e.g. u[1], which is bad. To Do: Make change so that these can be used.

#Coordination function, actually does all the work of the macro.
function coordinate(name, ex::Expr, p)
    reactions = get_reactions(ex)           ::Vector{ReactionStruct}
    reactants = get_reactants(reactions)    ::OrderedDict{Symbol,Int64}
    parameters = get_parameters(p)          ::OrderedDict{Symbol,Int64}

    syms = collect(keys(reactants))
    params = collect(keys(parameters))

    f_expr = get_f(reactions, reactants)
    f = make_func(f_expr, reactants, parameters)

    g_expr = get_g(reactions, reactants)
    g = make_func(g_expr, reactants, parameters)
    p_matrix = zeros(length(reactants), length(reactions))

    (jump_rate_expr, jump_affect_expr) = get_jump_expr(reactions, reactants)
    jumps = get_jumps(jump_rate_expr, jump_affect_expr,reactants,reactants)

    f_rhs = [element.args[2] for element in f_expr]
    symjac = Expr(:quote, calculate_jac(f_rhs, syms))
    f_symfuncs = hcat([SymEngine.Basic(f) for f in f_rhs])

    # Build the type
    exprs = Vector{Expr}(0)

    ## only get the right-hand-side of the equations.
    f_funcs = [element.args[2] for element in f_expr]
    g_funcs = [element.args[2] for element in g_expr]

    typeex,constructorex = maketype(name, f, f_funcs, f_symfuncs, g, g_funcs, jumps, Meta.quot(jump_rate_expr), Meta.quot(jump_affect_expr), p_matrix, syms; params=params, symjac=symjac)
    push!(exprs,typeex)
    push!(exprs,constructorex)

    # export type constructor
    def_const_ex = :(($name)()) |> esc
    push!(exprs,def_const_ex)

    expr_arr_to_block(exprs)
end

#Generates a vector containing a number of reaction structures, each containing the infromation about one reaction.
function get_reactions(ex::Expr)
    reactions = Vector{ReactionStruct}(0)      ::Vector{ReactionStruct}
    for line in ex.args
        (line.head != :tuple) && (continue)
        (rate,r_line) = line.args

        if r_line.head  == :-->
            r_line = Expr(:call,:→,r_line.args[1],r_line.args[2])
        end

        arrow = r_line.args[1]  ::Symbol
        if in(arrow,double_arrows)
            push_reactions(reactions::Vector{ReactionStruct}, r_line.args[2], r_line.args[3], rate.args[1], !in(arrow,no_mass_arrows))
            push_reactions(reactions::Vector{ReactionStruct}, r_line.args[3], r_line.args[2], rate.args[2], !in(arrow,no_mass_arrows))
        elseif in(arrow,fwd_arrows)
            push_reactions(reactions::Vector{ReactionStruct}, r_line.args[2], r_line.args[3], rate, !in(arrow,no_mass_arrows))
        elseif in(arrow,bwd_arrows)
            push_reactions(reactions::Vector{ReactionStruct}, r_line.args[3], r_line.args[2], rate, !in(arrow,no_mass_arrows))
        else
            throw("malformed reaction")
        end
    end
    return reactions
end

#Structure containing information about one reactant in one reaction.
struct ReactantStruct
    reactant::Symbol
    stoichiometry::Int64
end

#Structure containing information about one Reaction. Contain all its substrates and products as well as its rate. Contains an specialized constructor.
struct ReactionStruct
    substrates::Vector{ReactantStruct}
    products::Vector{ReactantStruct}
    rate::Any
    function ReactionStruct(sub_line::Any, prod_line::Any, rate::Any, use_mass_kin::Bool)
        sub = add_reactants!(sub_line,1,Vector{ReactantStruct}(0))
        prod = add_reactants!(prod_line,1,Vector{ReactantStruct}(0))
        use_mass_kin && (rate = mass_rate(sub,rate))
        new(sub,prod,rate)
    end
end

#If we want to use mass kinetics, modifies rate accordingly. Called in ReactionStruct constructor if use_mass_kin is true.
function mass_rate(substrates::Vector{ReactantStruct},old_rate::Any)
    rate = Expr(:call, :*, old_rate)
    for sub in substrates
        push!(rate.args,Expr(:call, :^, sub.reactant, sub.stoichiometry))
    end
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
function add_reactants!(ex::Any, mult::Int64, reactants::Vector{ReactantStruct})
    if typeof(ex)!=Expr
        (ex == 0 || in(ex,empty_set)) && (return reactants)
        in(ex,disallowed_reactants) && throw("Can not use reactant names: u, du, p, t. These are used in function arguments.")
        push!(reactants, ReactantStruct(ex,mult))
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
    reactants = OrderedDict{Symbol,Int64}()
    r_count = 0    ::Int64
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
    parameters = OrderedDict{Symbol,Int64}()
    p_count = 0    ::Int64
    for parameter in p
        (!haskey(parameters,parameter)) && (parameters[parameter] = p_count += 1)
    end
    return parameters
end

#From the reactions and reactants generates f, the functions describing the deterministic time evolution of the system.
function get_f(reactions::Vector{ReactionStruct}, reactants::OrderedDict{Symbol,Int64})
    f = Vector{Expr}(length(reactants))
    for i = 1:length(f)
        f[i] = :(du[$i] = $(Expr(:call, :+)))
    end

    #Loops through all reactions. For all products and substrates loops ads their rate of change to the corresponding line in the system (off differential equations).
    for reaction in deepcopy(reactions)
        for prod in reaction.products
            push!(f[reactants[prod.reactant]].args[2].args, recursive_clean!(:($(reaction.rate) * $(prod.stoichiometry))))
        end
        for sub in reaction.substrates
            push!(f[reactants[sub.reactant]].args[2].args, recursive_clean!(:(-$(reaction.rate) * $(sub.stoichiometry))))
        end
    end
    return f
end

function get_g(reactions::Vector{ReactionStruct}, reactants::OrderedDict{Symbol,Int64})
    g = Vector{Expr}(length(reactions)*length(reactants))
    idx = 0
    for reactant in keys(reactants), i = 1:length(reactions)
            g[idx += 1] = recursive_clean!(:(du[$(reactants[reactant]),$i] = $(get_stoch_diff(reactions[i],reactant)) * sqrt($(deepcopy(reactions[i].rate)))))
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

function make_func(func_expr::Vector{Expr},reactants::OrderedDict{Symbol,Int64},parameters::OrderedDict{Symbol,Int64})
    system = Expr(:block)
    for func_line in deepcopy(func_expr)
        push!(system.args, recursive_replace!(func_line, (reactants,:u), (parameters, :p)))
    end
    return :((du,u,p,t) -> $system)
end

function get_jump_expr(reactions::Vector{ReactionStruct}, reactants::OrderedDict{Symbol,Int64})
    rates = Vector{Any}(length(reactions))
    affects = Vector{Vector{Expr}}(length(reactions))
    idx = 0
    for reaction in deepcopy(reactions)
        rates[idx += 1] = recursive_clean!(reaction.rate)
        affects[idx] = Vector{Expr}(0)
        foreach(prod -> push!(affects[idx],:(integrator.u[$(reactants[prod.reactant])] += $(prod.stoichiometry))), reaction.products)
        foreach(sub -> push!(affects[idx],:(integrator.u[$(reactants[sub.reactant])] -= $(sub.stoichiometry))), reaction.substrates)
    end
    return (Tuple(rates),Tuple(affects))
end

function get_jumps(rates::Tuple, affects::Tuple,reactants::OrderedDict{Symbol,Int64},parameters::OrderedDict{Symbol,Int64})
    jumps = Expr(:tuple)
    for i = 1:length(rates)
        push!(jumps.args,Expr(:call,:ConstantRateJump))
    end
    for i = 1:length(rates)
        push!(jumps.args[i].args, :((u,p,t) -> $(recursive_replace!(deepcopy(rates[i]), (reactants,:u), (parameters, :p)))))
        push!(jumps.args[i].args, :(integrator -> $(expr_arr_to_block(deepcopy(affects[i])))))
    end
    return jumps
end

function recursive_clean!(expr::Any)
    (typeof(expr)!=Expr) && (return expr)
    for i = 1:length(expr.args)
        expr.args[i] = recursive_clean!(expr.args[i])
    end
    (expr.args[1] == :^) && (expr.args[3] == 1) && (return expr.args[2])
    if expr.args[1] == :*
        for i = length(expr.args):-1:2
            (expr.args[i] == 1) && deleteat!(expr.args,i)                   #Removes all multiplications by 1.
        end
        (length(expr.args) == 2) && (return expr.args[2])                   # We have a multiplication of only one thing, return only that thing.
        (length(expr.args) == 1) && (return 1)                              #We have only * and no real argumenys.
    end
    if expr.head == :call
        in(expr.args[1],hill_name) && return hill(expr)
        (in(expr.args[1],mm_name)) && (return mm(expr))
    end
    return expr
end

function recursive_replace!(expr::Any, replace_requests::Tuple{OrderedDict{Symbol,Int64},Symbol}...)
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

#Makes the Jacobian
function calculate_jac(f_expr::Vector{Expr}, syms)
    symjac = Matrix{SymEngine.Basic}( length(syms), length(syms))
    symfuncs = [SymEngine.Basic(f) for f in f_expr]
    for i in eachindex(f_expr)
        for j in eachindex(syms)
          symjac[i,j] = diff(symfuncs[i],syms[j])
        end
    end
    symjac
end

#Turns an array of expressions to a expression block with corresponding expressions.
function expr_arr_to_block(exprs)
  block = :(begin end)
  foreach(expr -> push!(block.args, expr), exprs)
  block
end

#hill function made avaiable
hill_name = Set{Symbol}([:hill, :Hill, :h, :H, :HILL])
function hill(expr::Expr)
    return :($(expr.args[3])*($(expr.args[2])^$(expr.args[5]))/($(expr.args[4])^$(expr.args[5])+$(expr.args[2])^$(expr.args[5])))
end

#michaelis menten function made avaiable.
mm_name = Set{Symbol}([:MM, :mm, :Mm, :mM, :M, :m])
function mm(expr::Expr)
    return :($(expr.args[3])*$(expr.args[2])/($(expr.args[4])+$(expr.args[2])))
end
