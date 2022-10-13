### Fetch Required Packages ###

# Very WIP (wanted to do with {S,T} and stuff, but got various errors, and wasn't needed right now).
### Spatial Reaction Structure. ###
### Describing a spatial reaction that involves species from two neighbouring compartments. ###
struct SpatialReaction
    """The rate function (excluding mass action terms)."""
    rate::Any
    """Reaction substrates."""
    substrates::Tuple{Vector,Vector}
    """Reaction products."""
    products::Tuple{Vector,Vector}
    """
    `false` (default) if `rate` should be multiplied by mass action terms to give the rate law.
    `true` if `rate` represents the full reaction rate law.
    """
    only_use_rate::Bool

    function SpatialReaction(rate,substrates,products;only_use_rate=false)
        new(rate,substrates,products,only_use_rate)
    end
    
end

### Lattice Reaction Network Structure ###
# Couples:
# A network of compartments (a meta graph that can contain some additional infro for each compartment).
# A reaction network (that is simulated within each compartment).
# A set of spatial reactions (denoting interaction between comaprtments).
struct LatticeReactionSystem # <: MT.AbstractTimeDependentSystem # Adding this part messes up show, disabling me from creating LRSs
    """The spatial reactions defined between individual nodes."""
    rs::ReactionSystem
    """The spatial reactions defined between individual nodes."""
    spatial_reactions::Vector{SpatialReaction}
    """The graph on which the lattice is defined."""
    lattice::DiGraph
    """Dependent (state) variables representing amount of each species. Must not contain the
    independent variable."""
    
    function LatticeReactionSystem(rs, spatial_reactions, lattice::DiGraph)
        return new(rs, spatial_reactions, lattice)
    end
    function LatticeReactionSystem(rs, spatial_reactions, lattice::SimpleGraph)
        return new(rs, spatial_reactions, DiGraph(lattice))
    end
end


# Function for creating a single ReactionSystem structure from a LatticeReactionSystem input.
function make_spatial_reaction_system(lrs::LatticeReactionSystem)
    @unpack rs,spatial_reactions,lattice = lrs

    # Extends an empty base network, with one copy of the input network for each node on the graph.    
    @named rs_base = ReactionSystem(rs.iv)
    rs_combined = flatten(compose(rs_base, map(c -> make_subsystem(rs,c),vertices(lattice))))
    
    connection_systems = Vector{ReactionSystem}()
    for sr in spatial_reactions, (ei,e) in enumerate(edges(lattice))
        substrates_src = map(sub -> ParentScope(sym_to_var(sub[1],rs_combined,e.src)),sr.substrates[1])
        substrates_dst = map(sub -> ParentScope(sym_to_var(sub[1],rs_combined,e.dst)),sr.substrates[2])
        products_src = map(prod -> ParentScope(sym_to_var(prod[1],rs_combined,e.src)),sr.products[1])
        products_dst = map(prod -> ParentScope(sym_to_var(prod[1],rs_combined,e.dst)),sr.products[2])
        subs_stoich_src = map(sub -> sub[2],sr.substrates[1])
        subs_stoich_dst = map(sub -> sub[2],sr.substrates[2])
        prods_stoich_src = map(prod -> prod[2],sr.products[1])
        prods_stoich_dst = map(prod -> prod[2],sr.products[2])
        rs_internal = ReactionSystem([Reaction(sr.rate,[substrates_src;substrates_dst],[products_src;products_dst],[subs_stoich_src;subs_stoich_dst],[prods_stoich_src;prods_stoich_dst])],rs_base.iv;name=Symbol(:connection___,ei))
        push!(connection_systems,rs_internal)
    end
    return flatten(compose(rs_combined, connection_systems))
end

# Helper functions.
make_subsystem(rs::ReactionSystem,comp::Int64) = @set! rs.name = Symbol("container___",comp)
sym_to_var(sym,rs,comp) = get_var_to_name(rs)[Symbol("container___",comp,:₊,sym)]

### ODEProblem ###
# Creates an ODEProblem from a LatticeReactionSystem.
# Currently parameter values are identical across the grid.
# Default species concentrations are set in the input to the problem. 
# If a specific default value is given in the input value in the graph, that is used isntead.
function DiffEqBase.ODEProblem(lrs::LatticeReactionSystem, u0, tspan,
                               p = DiffEqBase.NullParameters(), args...;
                               check_length = false, name = nameof(lrs.rs),
                               combinatoric_ratelaws = get_combinatoric_ratelaws(lrs.rs),
                               include_zero_odes = true, remove_conserved = false,
                               checks = false, kwargs...)
                               
    # Creates the expanded reaction system
    rs_expanded = make_spatial_reaction_system(lrs)

    # Creates spatial u0 and p vectors
    container_syms = keys(get_var_to_name(lrs.rs))
    nCells = length(vertices(lrs.lattice))
    nConnections = length(edges(lrs.lattice))
    u0_expanded = vcat(map(u0_i -> expand_sym_val_pair(u0_i,container_syms,nCells,nConnections), u0)...)
    p_expanded = vcat(map(p_i -> expand_sym_val_pair(p_i,container_syms,nCells,nConnections), p)...)
   
    u0map = symmap_to_varmap(rs_expanded,u0_expanded)
    pmap = symmap_to_varmap(rs_expanded,p_expanded)

    osys = convert(ODESystem, rs_expanded; name, combinatoric_ratelaws, include_zero_odes, checks,
                   remove_conserved)
    return ODEProblem(osys, u0map, tspan, pmap, args...; check_length, kwargs...)
end

# Helper function.
# Converts e.g. :x => 1.0 to [container___1₊X => 1.0, container___2₊X => 1.0 ...]
function expand_sym_val_pair(pair, container_syms, nCells, nConnections)
    base,n = in(pair[1],container_syms) ? (:container___,nCells) : (:connection___,nConnections)
    syms = map(i -> Symbol(base,i,:₊,pair[1]), 1:n)
    vals = (pair[2] isa Vector) ? pair[2] : fill(pair[2],n)
    return Pair.(syms,vals)
end


### SDEProblem ###
# Basically a copy of ODEProblem
function DiffEqBase.SDEProblem(lrs::LatticeReactionSystem, u0, tspan,
    p = DiffEqBase.NullParameters(), args...;
    check_length = false, name = nameof(rs),
    combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
    include_zero_odes = true, remove_conserved = false,
    checks = false, kwargs...)
                               
    # Creates the expanded reaction system
    rs_expanded = make_spatial_reaction_system(lrs)

    # Creates spatial u0 and p vectors
    container_syms = keys(get_var_to_name(lrs.rs))
    nCells = length(vertices(lrs.lattice))
    nConnections = length(edges(lrs.lattice))
    u0_expanded = vcat(map(u0_i -> expand_sym_val_pair(u0_i,container_syms,nCells,nConnections), u0)...)
    p_expanded = vcat(map(p_i -> expand_sym_val_pair(p_i,container_syms,nCells,nConnections), p)...)
   
    u0map = symmap_to_varmap(rs_expanded,u0_expanded)
    pmap = symmap_to_varmap(rs_expanded,p_expanded)

    sde_sys = convert(SDESystem, rs_expanded; name, combinatoric_ratelaws, include_zero_odes, checks,remove_conserved)
    p_matrix = zeros(length(get_states(sde_sys)), numreactions(lrs.rs))
    return SDEProblem(sde_sys, u0map, tspan, pmap, args...; check_length, noise_rate_prototype = p_matrix, kwargs...)
end

### Discrete and Jump Problems ###
# DiscreteProblem from AbstractReactionNetwork
function DiffEqBase.DiscreteProblem(lrs::LatticeReactionSystem, u0, tspan::Tuple,
    p = DiffEqBase.NullParameters(), args...;
    name = nameof(rs),
    combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
    checks = false, kwargs...)
                               
    # Creates the expanded reaction system
    rs_expanded = make_spatial_reaction_system(lrs)

    # Creates spatial u0 and p vectors
    container_syms = keys(get_var_to_name(lrs.rs))
    nCells = length(vertices(lrs.lattice))
    nConnections = length(edges(lrs.lattice))
    u0_expanded = vcat(map(u0_i -> expand_sym_val_pair(u0_i,container_syms,nCells,nConnections), u0)...)
    p_expanded = vcat(map(p_i -> expand_sym_val_pair(p_i,container_syms,nCells,nConnections), p)...)
   
    u0map = symmap_to_varmap(rs_expanded,u0_expanded)
    pmap = symmap_to_varmap(rs_expanded,p_expanded)
 
    jsys = convert(JumpSystem, rs_expanded; name, combinatoric_ratelaws, checks)
    return DiscreteProblem(jsys, u0map, tspan, pmap, args...; kwargs...)
end

# JumpProblem from AbstractReactionNetwork
function JumpProcesses.JumpProblem(lrs::LatticeReactionSystem, prob, aggregator, args...;
   name = nameof(rs),
   combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
   checks = false, kwargs...)

    jsys = convert(JumpSystem, make_spatial_reaction_system(lrs); name, combinatoric_ratelaws, checks)
    return JumpProblem(jsys, prob, aggregator, args...; kwargs...)
end

### Plotting ###

# Should at some point overload "plot()" command.
# Plots the output from a spatial simulation as an animation over time.
function plot_spatial_sol(sol,lattice,var;v_max=maximum(maximum.(sol.u)),samp_freq=1,nodesize=0.3,linewidth=1,kwargs...)
    nComps = length(vertices(lattice));
    trajectories = map(vals -> vals[(var-1)*nComps+1:var*nComps] ./ v_max, sol.u)[1:samp_freq:end]

    # Has to do this to fix x and y positions of all nodes, else the graph moves in every fram (=very bad).
    base_plot = graphplot(lattice); 
    xs = base_plot.series_list[end][:x]
    ys = base_plot.series_list[end][:y]

    # Function that plots the graph for a certain set of solution values.
    function plot_frame(vals)
        graphplot(lattice,x=xs,y=ys,
        nodeshape=:circle, nodesize=nodesize,
        axis_buffer=0.6,
        curves=false,
        color=:black,
        nodecolor=map(val -> RGB{Float64}(val,val,0.0), vals),
        linewidth=linewidth,kwargs...)
    end

    # Creates an animation.
    return  @animate for vals in trajectories
        plot_frame(vals)
    end
end