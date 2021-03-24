

# This function takes a graph and a set of indices
# and returns an iterator over assignments
# i -> e[1], j -> e[2], k -> e[3], ...
# for all n-edges of the graph, where n is the
# number of indices. The assignments should be
# of type Dict{Symbol,Symbol}.
#
# For an ordinary graph, the 1-edges are the vertices
# and the 2-edges are the edges in the usual sense.
#
# Important: as the order of the indices is not
# defined, this function should assign elements
# to all permutations of the edges, ie. for
# i => 4, j => 5 it should also return i => 5, j => 4
function graph_iterator end

struct BasicGraph{T}
    verts::AbstractVector{T}
    edges::AbstractVector{Pair{T,T}}
end
    
function graph_iterator(graph::BasicGraph{T}, idcs::Set{Symbol}) where {T}
    if length(idcs) == 1
        i = first(idcs)
        return ( Dict(i => vert) for vert in graph.verts )
    elseif length(idcs) == 2
        i, j = collect(idcs)
        edges_r = ( edge.second => edge.first for edge in graph.edges )
        edges_sym = union(graph.edges, edges_r)
        
        return ( Dict(i => edge.first, j => edge.second) for edge in edges_sym )
    end
    
    return nothing
end

struct PeriodicGrid{N} 
    size::NTuple{N, Int64}
end

PeriodicGrid(sizes...) = PeriodicGrid(sizes)

function get_neighbour(coords, dim::Int64, amt::Int64, size::NTuple{N,Int64})::Vector{Int64} where N
    ret = collect(coords)
    ret[dim] += amt
    ret[dim] = mod(ret[dim] - 1, size[dim]) + 1    # transform to lie in 1:size[dim]
    ret
end

coord_to_str(coords) = "(" * join(coords, ",") * ")"
grid_coords(grid::PeriodicGrid{N}) where N = Iterators.product([ 1:grid.size[i] for i in 1:N ]...)
    
function get_edge_iter(i::Symbol, j::Symbol, dim::Int64, grid::PeriodicGrid)
    return Iterators.flatten(( Dict(i => coord_to_str(coords), 
                                      j => coord_to_str(get_neighbour(coords, dim, sgn, grid.size))) 
                               for coords in grid_coords(grid) ) 
                               for sgn in [-1,+1])
end

function graph_iterator(grid::PeriodicGrid, idcs::Set{Symbol})
    if length(idcs) == 1
        idx = first(idcs)
        return ( Dict(idx => coord_to_str(coords)) for coords in grid_coords(grid) ) 
    elseif length(idcs) == 2
        i, j = collect(idcs)
        return Iterators.flatten(get_edge_iter(i, j, dim, graph) for dim in 1:N)
    end
    
    return nothing
end

# Convert a list of expressions of the form
# G[i], H[j,k] into a dictionary 
# of the form Dict(i => G, j => H, k => H)
function get_idx_dict(exprs)::Dict{Symbol,Any}
    ret = Dict{Symbol,Any}()
    
    for idx_map in exprs
        (isa(idx_map, Expr) && idx_map.head == :ref) || throw("malformed index assignment: $(idx_map)")
        graph_name = idx_map.args[1]
        graph = Base.MainInclude.eval(graph_name)
        idcs = idx_map.args[2:end]
        
        for idx in idcs
            idx isa Symbol || throw("index not a symbol: $(idx)")
            idx in keys(ret) && throw("multiple index assignment: $(idx)")
            ret[idx] = graph
        end
    end
    
    return ret
end

# Iterate over all valid combinations of indices for a spatial reaction
function iterate_idcs(idcs::Set{Symbol}, idx_dict::AbstractDict)
    isempty(idcs) && return (Dict{Symbol,Any}(),)
    
    iterators = []
    for graph in unique(values(idx_dict))
        idcs_graph = filter(idx -> idx_dict[idx] == graph, idcs)
        isempty(idcs_graph) && continue
        iterator = graph_iterator(graph, idcs_graph)
        iterator === nothing && throw("graph iteration not implemented for $(length(idcs_graph)) indices")
        push!(iterators, iterator)
    end
    
    return (merge(dicts...) for dicts in Iterators.product(iterators...))
end

# Get indices occurring in a reaction
function get_idcs(ex::ExprValues)::Set{Symbol}
    if ex isa Expr
        if ex.head == :ref
            return Set(ex.args[2:end])
        end
        
        return mapreduce(get_idcs, union, ex.args)
    end
    
    return Set{Symbol}()
end
        
# Expand spatial reactants of the form X[i] with the given assignment of
# vertices to indices
function expand_reactant_idcs(ex::ExprValues, idx_map::AbstractDict)
    if ex isa Expr
        if ex.head == :ref
            name = ex.args[1]
            name isa Symbol || throw("only species can be indexed")
                    
            length(ex.args) == 2 || throw("species $(ex.args[1]) can only have one index")
            
            idx = ex.args[2]
            idx_mapped = idx_map[idx]
            
            return Symbol(name, "_", idx_mapped)
        end
        
        return Expr(ex.head, map(arg -> expand_reactant_idcs(arg, idx_map), ex.args)...)
    end
    
    ex
end

# Expand reaction rates of the form a[i,j] with the given assignment of
# vertices to indices
function expand_rate_idcs(ex::ExprValues, idx_map::AbstractDict)
    if ex isa Expr
        if ex.head == :ref
            name = ex.args[1]
            name isa Symbol || throw("only constant parameters can be indexed")
            
            idcs = ex.args[2:end]
            idcs_mapped = map(x -> idx_map[x], idcs)
            if length(idcs) == 1
                return Symbol(name, "_", idcs_mapped[1])
            else
                return Symbol(name, "_(", join(idcs_mapped, ","), ")")
            end
        end
        
        return Expr(ex.head, map(arg -> expand_rate_idcs(arg, idx_map), ex.args)...)
    end
    
    ex
end
        
function expand_spatial_reactions!(reactions::AbstractVector{Expr}, rate::ExprValues, sub_line::ExprValues, prod_line::ExprValues, arrow::Symbol, 
                                   idx_dict::AbstractDict)
    idcs_sub = get_idcs(sub_line)
    idcs_prod = get_idcs(prod_line)
    idcs = union(idcs_sub, idcs_prod)
    
    issubset(idcs, keys(idx_dict)) || throw("invalid index in reaction")
    
    idcs_rates = get_idcs(rate)
    issubset(idcs_rates, idcs) || throw("invalid index in reaction rate")
    
    for idx_map in iterate_idcs(idcs, idx_dict)
        rate_mapped = expand_rate_idcs(rate, idx_map)
        sub_line_mapped = expand_reactant_idcs(sub_line, idx_map)
        prod_line_mapped = expand_reactant_idcs(prod_line, idx_map)
        
        r_line = Expr(:call, arrow, sub_line_mapped, prod_line_mapped)

        push!(reactions, Expr(:tuple, rate_mapped, r_line))
    end
end
        
function get_spatial_reactions(ex::Expr, idx_dict::AbstractDict)
    reactions = Expr[]
    
    for line in ex.args
        (line isa Expr && line.head == :tuple) || continue
        
        (rate, r_line) = line.args
        (r_line.head  == :-->) && (r_line = Expr(:call,:→,r_line.args[1],r_line.args[2]))

        expand_spatial_reactions!(reactions, rate, r_line.args[2], r_line.args[3], r_line.args[1], idx_dict)
    end
    
    Expr(:block, reactions...)
end
        
macro spatial_reaction_network(ex::Expr, graphs...)
    idx_dict = get_idx_dict(graphs)
    esc(get_spatial_reactions(MacroTools.striplines(ex), idx_dict))
end