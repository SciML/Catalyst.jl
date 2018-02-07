#To Do:
# - Allow for several functions in one line, e.g:
#     kDeg, (X,Y,Z) --> ∅
#     (kDegX, kDegY, kDegZ), (X,Y,Z) --> ∅
# - Modify variable replacement so that p, t, u and du can be used as input variable names or reaction rates.
#     Low prio.
# - Add support for converting to polynomials and use algebraic methods (fincing equilibirumpoints etc.)
#     Long term.

#Declare various arrow types symbols used for the empty set (also 0).
empty_set = Set{Symbol}([:∅])
fwd_arrows = Set{Symbol}([:>, :→, :↣, :↦, :⇾, :⟶, :⟼, :⥟, :⥟, :⇀, :⇁, :⇒, :⟾])
bwd_arrows = Set{Symbol}([:<, :←, :↢, :↤, :⇽, :⟵, :⟻, :⥚, :⥞, :↼, :↽, :⇐, :⟽])
double_arrows = Set{Symbol}([:↔, :⟷, :⥎, :⥐, :⇄, :⇆, :⇋, :⇌, :⇔, :⟺])
no_mass_arrows = Set{Symbol}([:⇐, :⟽, :⇒, :⟾, :⇔, :⟺])      #Using this arrows will disable the program from multiplying reaction rates with the substrate concentrations. GIves user full control of reaction rates.
disallowed_reactants = Set{Symbol}([:du, :u, :p, :t])           #These are not allowed since they are used in "return :((du,u,p,t) -> $system)", if a variable these gets replaced with e.g. u[1], which is bad.

"""
Macro that inputs an expression corresponding to a reaction entwork and output a Reaction Netqork Structure that can be used as input to generation of SDE and ODE and Jump problems.
Must arrows accepted (both right, left and bi drectional arrows).
Using arrows no fileld arrows (⇐, ⟽, ⇒, ⟾, ⇔, ⟺) will disable mass kinetics and lets you cutomize reaction rates yourself.
Example system:
    2.0, X + Y --> XY       #This will have reaction rate corresponding to 2.0*[X][Y]
    2.0X, X + Y ⟾ XY       #This will have reaction rate corresponding to 2.0*[X]
    (hill(X,2,2,2),kD), X + Y ⟷ XY    #Reaction in forward direction is activated by X according to a hill function. Reaction in backward direction have a rate according to constant kD, declared elsewere in your program.
Note that while --> is a correct arrow, neither <-- nor <--> works.
"""
#Macro to create a reaction network model. Generates expressions for the various things you might want. Last line is executed and constructions a ReactionNetwork structure containing the infromation,
macro reaction_network_new(ex::Expr, p...)
    reactions = get_reactions(ex)           ::Vector{ReactionStruct}
    reactants = get_reactants(reactions)    ::Dict{Symbol,Int64}
    parameters = get_parameters(p)          ::Dict{Symbol,Int64}

    f = recursive_equify!(get_f(reactions, reactants), reactants, parameters)                                                ::Expr   #For ODEs
    f_expr = Expr(:quote,recursive_equify!(get_f(reactions, reactants), Dict{Symbol,Int64}(), Dict{Symbol,Int64}()))         ::Expr   #For ODEs
    g = recursive_equify!(get_g(reactions, reactants), reactants, parameters)                                                ::Expr   #For SDEs
    g_expr = Expr(:quote,recursive_equify!(get_g(reactions, reactants), Dict{Symbol,Int64}(), Dict{Symbol,Int64}()))         ::Expr   #For SDEs
    jumps = recursive_equify!(get_jumps(reactions, reactants), reactants, parameters)                                        ::Expr   #For Gillespie Simulations
    jumps_expr = Expr(:quote,recursive_equify!(get_jumps(reactions, reactants), Dict{Symbol,Int64}(), Dict{Symbol,Int64}())) ::Expr   #For Gillespie Simulations
    return :(ReactionNetwork($f, $f_expr, $g, $g_expr, $jumps, $jumps_expr, zeros($(length(reactants)),$(length(reactions)))))
end

#Generates a dictionary with all parameters.
function get_parameters(p)
    parameters = Dict{Symbol,Int64}()
    p_count = 0    ::Int64
    for parameter in p
        (!haskey(parameters,parameter)) && (parameters[parameter] = p_count += 1)
    end
    return parameters
end

#Generates a vector containing a number of reaction structures, each containing the infromation about one reaction.
function get_reactions(ex::Expr)
    reactions = Vector{ReactionStruct}(0)      ::Vector{ReactionStruct}     #The reactions are saved here.
    for line in ex.args
        (line.head != :tuple) && (continue)
        (rate,r_line) = line.args

        #Allows --> to be used in addition to normal arrows. This gives a different expression so have to be thrown around a little. Arrows  <--> and <-- do not yield correct expressions and cannot be used.
        if r_line.head  == :-->
            r_line = Expr(:call,:→,r_line.args[1],r_line.args[2])
        end

        #Checks what type of arrow we have (what direction) and generates reactions accordingly.
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

#Structure containing information about one Reaction. Contain all its substrates and products as well as its rate.
struct ReactionStruct
    substrates::Vector{ReactantStruct}
    products::Vector{ReactantStruct}
    rate::Any
    #Construction. Genearates a reaction from one line of expression.
    #use_mass_kin = true will multiple the reaction rate by the concentration of the substrates (one usually want this).
    #direction says which direction the reaction goes in (X --> Y is true, X <-- Y is false. For X <--> Y this will be called in each direction).
    function ReactionStruct(sub_line::Any, prod_line::Any, rate::Any, use_mass_kin::Bool)
        sub = add_reactants!(sub_line,1,Vector{ReactantStruct}(0))
        prod = add_reactants!(prod_line,1,Vector{ReactantStruct}(0))
        use_mass_kin && (rate = mass_rate(sub,rate))
        new(sub,prod,rate)
    end
end

#If we want to use mass kinetics, fixes that.
function mass_rate(substrates::Vector{ReactantStruct},old_rate::Any)
    rate = Expr(:call, :*, old_rate)
    for sub in substrates
        push!(rate.args,Expr(:call, :^, sub.reactant, sub.stoichiometry))
    end
    return rate
end

function tup_leng(ex::Any)
    (typeof(ex)==Expr && ex.head == :tuple) && (return length(ex.args))
    return 1
end

function get_tup_arg(ex::Any,i::Int)
    (tup_leng(ex) == 1) && (return ex)
    return ex.args[i]
end

function push_reactions(reactions::Vector{ReactionStruct}, sub_line::Any, prod_line::Any, rate::Any, use_mass_kin::Bool)
    lengs = [tup_leng(sub_line), tup_leng(prod_line), tup_leng(rate)]
    (count(lengs.==1) + count(lengs.==maximum(lengs)) < 3) && (throw("malformed reaction"))
    for i = 1:maximum(lengs)
        push!(reactions, ReactionStruct(get_tup_arg(sub_line,i), get_tup_arg(prod_line,i), get_tup_arg(rate,i), use_mass_kin))
    end
end

#Recursive function that loops through the reactants in an reaction line and finds the reactants and their stochiometry. Recursion makes it able to handle e.g. 2(X+Y+3(Z+XY)) (probably one will not need it though).
function add_reactants!(ex::Any, mult::Int64, reactants::Vector{ReactantStruct})
    #We have a symbol (or 0), we have found a reactant or possibly the empty set.
    if typeof(ex)!=Expr
        (ex == 0 || in(ex,empty_set)) && (return reactants)
        in(ex,disallowed_reactants) && throw("Can not use reactant names: u, du, p, t. These are used in function arguments.")
        push!(reactants, ReactantStruct(ex,mult))
    elseif ex.args[1] == :*         #We have found something on the form 2X, gets the stochiometry and recal. The recal will be on a sybol. Alterantive is that we have e.g. 2(X+3Y) for this reason stochiometry is stored in recall.
        add_reactants!(ex.args[3],mult*ex.args[2],reactants)
    elseif ex.args[1] == :+            # We have a sum of reactants, loops through all reactants.
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
    reactants = Dict{Symbol,Int64}()
    r_count = 0    ::Int64
    #For all reactions, checks all products and substrates. Add them to the dictionary (if not already in it) and updates countr of number of reactant types.
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

#From the reactions and reactants generates f, the functions describing the deterministic time evolution of the system.
function get_f(reactions::Vector{ReactionStruct}, reactants::Dict{Symbol,Int64})
    #Generates the system base.
    system = Expr(:block)

    #Ensures every line start something like du[2] = ...
    #Here ... is a sum of some terms (to be added to expression in the next step).
    for i = 1:length(reactants)
        line = :(du[$i] = $(Expr(:call, :+)))
        push!(system.args,line)
    end

    #Loops through all reactions. For all products and substrates loops ads their rate of change to the corresponding line in the system (off differential equations).
    for reaction in reactions
        for prod in reaction.products
            push!(system.args[reactants[prod.reactant]].args[2].args, :($(reaction.rate) * $(prod.stoichiometry)))
        end
        for sub in reaction.substrates
            push!(system.args[reactants[sub.reactant]].args[2].args, :(-$(reaction.rate) * $(sub.stoichiometry)))
        end
    end
    return :((du,u,p,t) -> $system)
end

#From the reactions and reactants generates g, the functions describing the noise of the system for Gillespie SDE simulations (noise multiplies by sqrt of reaction rate and the stochiometric change due to a certain reaction).
function get_g(reactions::Vector{ReactionStruct}, reactants::Dict{Symbol,Int64})
    system = Expr(:block)               #Creates an empty system to put the equations in.
    for reactant in keys(reactants)     #For every reactan.
        for i = 1:length(reactions)     #For every reaction.
            line = :(du[$(reactants[reactant]),$i] = $(get_stoch_diff(reactions[i],reactant)) * sqrt($(reactions[i].rate)))     #Get and inserts an entry corresponding to the noise rate for that reactant in that reaction (0 of reactant not part of reaction).
            push!(system.args,line)
        end
    end
    return :((du,u,p,t) -> $system)
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

#Generates a tuple of constant rate jumps to be used for Gillespie simulations.
function get_jumps(reactions::Vector{ReactionStruct}, reactants::Dict{Symbol,Int64})
    return_tuple = :(())            #Expression for a tuple, contains all the ConatantRateJumps we create.
    for reaction in reactions
        system = Expr(:block)       #System corresponding to the affect functions changes in the integrator.
        #Loops through all substrates and products in the reaction and gets their stoichiometry change.
        for prod in reaction.products
            push!(system.args,:(integrator.u[$(reactants[prod.reactant])] += $(prod.stoichiometry)))
        end
        for sub in reaction.substrates
            push!(system.args,:(integrator.u[$(reactants[sub.reactant])] -= $(sub.stoichiometry)))
        end
        push!(return_tuple.args, :(ConstantRateJump((u,p,t) -> $(reaction.rate),integrator -> $system)))    #Adding a reaction to the return tuple. Contains the rate function and the affect function.
    end
    return return_tuple
end

#Recursive function that replaces all X, Y etc. with du[1], du[2] etc. Also removes stuff like 1*u[1]^1, this to make it easier to understand the equations if you print them.
function recursive_equify!(expr::Any, reactants::Dict{Symbol,Int64}, parameters::Dict{Symbol,Int64})
    if typeof(expr) == Symbol           #If we have a symbol, check if it is one of the reactants. If so replace it accordingly.
        (haskey(reactants,expr)) && (return :(u[$(reactants[expr])]))
        (haskey(parameters,expr)) && (return :(p[$(parameters[expr])]))
    elseif typeof(expr) == Expr         #If we have an expression, do recursion on its parts.
        for i = 1:length(expr.args)
            expr.args[i] = recursive_equify!(expr.args[i], reactants, parameters)
        end
        (expr.args[1] == :^) && (expr.args[3] == 1) && (return expr.args[2])    #If we have to the power of 1, skip that.
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
    end
    return expr
end

#The output structure, contains all the interesting information in the system.
struct ReactionNetwork
    f::Function
    f_expr::Expr
    g::Function
    g_expr::Expr
    jumps::Tuple{ConstantRateJump,Vararg{ConstantRateJump}}
    jumps_expr::Expr
    p_matrix::Array{Float64,2}
end


#hill function made avaiable
hill_name = Set{Symbol}([:hill, :Hill, :h, :H, :HILL])
function hill(expr::Expr)
    return :($(expr.args[4])*($(expr.args[2])^$(expr.args[3]))/($(expr.args[5])^$(expr.args[3])+$(expr.args[2])^$(expr.args[3])))
end

#michaelis menten function made avaiable.
mm_name = Set{Symbol}([:MM, :mm, :Mm, :mM, :M, :m])
function mm(expr::Expr)
    return :($(expr.args[3])*$(expr.args[2])/($(expr.args[4])+$(expr.args[2])))
end
