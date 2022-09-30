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
    """The reaction system expanded onto a graph."""
    rs::ReactionSystem
    """The original (non-spatial) reaction system that was expanded onto a graph."""
    rs_base::ReactionSystem
    """The spatial reactions defined between individual nodes."""
    spatial_reactions::Vector{SpatialReaction}
    """The graph on which the lattice is defined."""
    lattice::MetaGraph
    """Dependent (state) variables representing amount of each species. Must not contain the
    independent variable."""
    
    function LatticeReactionSystem(rs_base, spatial_reactions, lattice)
        rs = make_spatial_reaction_system(rs_base, spatial_reactions, lattice)
        return new(rs, rs_base, spatial_reactions, lattice)
    end
end


# Function for creating a single ReactionSystem structure from a LatticeReactionSystem input.
function make_spatial_reaction_system(rs::ReactionSystem,srs::Vector{SpatialReaction},lattice::MetaGraph)
    # Extends an empty base network, with one copy of the input network for each node on the graph.
    @named rs_base = ReactionSystem(rs.iv)
    rs_combined = flatten(compose(rs_base, map(c -> make_subsystem(rs,c,lattice),vertices(lattice))))

    # Loops through all connections, and spatial reactions, and add them (in both directions) for each connection.
    spatial_reactions = Vector{Reaction}()
    for sr in srs, e in edges(lattice)
        # Could create a loop to generate all 8 expressions, but code is easier to read this way.
        subs1_r1 = map(sub -> sym_to_var(sub,rs_combined,e.src,lattice), first.(sr.substrates[1]))
        subs2_r1 = map(sub -> sym_to_var(sub,rs_combined,e.dst,lattice), first.(sr.substrates[2]))
        subs1_r2 = map(sub -> sym_to_var(sub,rs_combined,e.dst,lattice), first.(sr.substrates[1]))
        subs2_r2 = map(sub -> sym_to_var(sub,rs_combined,e.src,lattice), first.(sr.substrates[2]))
        prods1_r1 = map(sub -> sym_to_var(sub,rs_combined,e.src,lattice), first.(sr.products[1]))
        prods2_r1 = map(sub -> sym_to_var(sub,rs_combined,e.dst,lattice), first.(sr.products[2]))
        prods1_r2 = map(sub -> sym_to_var(sub,rs_combined,e.dst,lattice), first.(sr.products[1]))
        prods2_r2 = map(sub -> sym_to_var(sub,rs_combined,e.src,lattice), first.(sr.products[2]))

        push!(spatial_reactions,Reaction(sr.rate,[subs1_r1;subs2_r1],[prods1_r1;prods2_r1],[last.(sr.substrates[1]);last.(sr.substrates[2])],[last.(sr.products[1]);last.(sr.products[2])]))
        push!(spatial_reactions,Reaction(sr.rate,[subs1_r2;subs2_r2],[prods1_r2;prods2_r2],[last.(sr.substrates[1]);last.(sr.substrates[2])],[last.(sr.products[1]);last.(sr.products[2])]))
    end

    # Adds the spatial reactions to the initial network.
    return extend(rs_combined,ReactionSystem(spatial_reactions,rs_combined.iv;name=rs_combined.name))
end

# Helper functions.
get_name(comp::Int64,lattice::MetaGraph) = has_prop(lattice,comp,:name) ? get_prop(lattice,comp,:name) : Symbol("container___",comp)
make_subsystem(rs::ReactionSystem,comp::Int64,lattice::MetaGraph) = @set! rs.name = get_name(comp,lattice)
sym_to_var(sym,rs,comp,lattice) = rs.var_to_name[Symbol(get_name(comp,lattice),:₊,sym)]

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
    
    # Converts the [:X => ...] to similar, but for each compartment variable.
    state_names = Symbol.(getfield.(states(lrs.rs_base),:f)) # Not sure if there's a better way to do this.
    ps_names = Symbol.(parameters(lrs.rs_base)) # Not sure if there's a better way to do this.
    u0_full = vcat(map(comp -> map(species -> sym_to_var(species,lrs.rs,comp,lrs.lattice) => get_val(comp,species,lrs.lattice,u0), state_names), vertices(lrs.lattice))...)
    p_full = vcat(map(comp -> map(species -> sym_to_var(species,lrs.rs,comp,lrs.lattice) => get_val(comp,species,lrs.lattice,p), ps_names), vertices(lrs.lattice))...)

    u0map = symmap_to_varmap(lrs.rs, u0_full)
    pmap_lattice = symmap_to_varmap(lrs.rs, p_full)
    pmap_global = symmap_to_varmap(lrs.rs,filter(p_pair->!in(p_pair[1],ps_names), p))

    osys = convert(ODESystem, lrs.rs; name, combinatoric_ratelaws, include_zero_odes, checks,
                   remove_conserved)
    return ODEProblem(osys, u0map, tspan, [pmap_lattice;pmap_global], args...; check_length, kwargs...)
end

# Helper function.
# Gets the value of a species (or parameter) for a given cell and ODEProblem input.
get_val(comp,sym,lattice,vals) = has_prop(lattice,comp,sym) ? get_prop(lattice,comp,sym) : Dict(vals)[sym] 


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